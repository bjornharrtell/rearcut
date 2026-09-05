//! `rearcut` is a Rust port of [mapbox/earcut](https://github.com/mapbox/earcut) (and
//! [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp)), a fast and small polygon
//! triangulation library based on ear slicing with z-order curve hashing for O(n) average
//! case performance.
//!
//! The algorithm triangulates simple polygons (which may be non-convex and contain holes)
//! given as one or more flat rings of 2D coordinates.
//!
//! # Example
//!
//! ```
//! use rearcut::earcut;
//!
//! // a quadrilateral: (10,0) (0,50) (60,60) (70,10)
//! let data = [10.0, 0.0, 0.0, 50.0, 60.0, 60.0, 70.0, 10.0];
//! let triangles: Vec<u32> = earcut(&data, &[], 2);
//! assert_eq!(triangles.len(), 6); // two triangles
//! assert!(rearcut::deviation(&data, &[], 2, &triangles) < 1e-9);
//! ```
//!
//! Polygons with holes are supported via `hole_indices`, which lists the starting vertex
//! index of each hole ring within `data`:
//!
//! ```
//! use rearcut::earcut;
//!
//! let data = [
//!     0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0, // outer ring
//!     20.0, 20.0, 80.0, 20.0, 80.0, 80.0, 20.0, 80.0, // hole
//! ];
//! let triangles: Vec<u32> = earcut(&data, &[4], 2);
//! assert_eq!(triangles.len(), 8 * 3);
//! ```

use std::fmt::Debug;

/// Sentinel used for arena-relative "null" links (there's no node at this index).
const NULL: u32 = u32::MAX;

/// A trait implemented for the integer types usable as output triangle indices.
///
/// Implemented for the common unsigned integer types. Panics (via `try_from`) if the
/// number of input vertices exceeds what the chosen index type can represent.
pub trait EarcutIndex: Copy + Debug {
    fn from_usize(v: usize) -> Self;
    fn to_usize(self) -> usize;
}

macro_rules! impl_earcut_index {
    ($($t:ty),*) => {
        $(
            impl EarcutIndex for $t {
                #[inline]
                fn from_usize(v: usize) -> Self {
                    <$t>::try_from(v).expect("vertex index does not fit in the chosen index type")
                }

                #[inline]
                fn to_usize(self) -> usize {
                    self as usize
                }
            }
        )*
    };
}

impl_earcut_index!(u16, u32, u64, usize);

/// top bit of `Node::i_steiner` flags a steiner point; the rest holds the vertex index.
const STEINER_BIT: u32 = 1 << 31;
const INDEX_MASK: u32 = !STEINER_BIT;

// `Node` is kept as small as possible (40 bytes, 8-byte aligned) so large rings fit more
// nodes per cache line: the steiner flag is packed into the top bit of the vertex index
// (mirrors upstream's `i_steiner` trick) instead of a separate `bool` field, and the vertex
// index itself is `u32` (an input with more than u32::MAX vertices would already exhaust
// memory for its coordinate array well before this becomes a limitation).
#[derive(Clone, Copy)]
struct Node {
    /// original vertex index in the input coordinates array (lower 31 bits) + steiner flag
    /// (bit 31)
    i_steiner: u32,
    /// z-order curve value
    z: i32,
    x: f64,
    y: f64,
    /// previous and next vertex nodes in a polygon ring
    prev: u32,
    next: u32,
    /// previous and next nodes in z-order (NULL if none)
    prev_z: u32,
    next_z: u32,
}

impl Node {
    #[inline]
    fn new(i: u32, x: f64, y: f64) -> Self {
        Node {
            i_steiner: i,
            x,
            y,
            prev: NULL,
            next: NULL,
            z: 0,
            prev_z: NULL,
            next_z: NULL,
        }
    }

    #[inline(always)]
    fn i(&self) -> u32 {
        self.i_steiner & INDEX_MASK
    }

    #[inline(always)]
    fn is_steiner(&self) -> bool {
        self.i_steiner & STEINER_BIT != 0
    }

    #[inline(always)]
    fn set_steiner(&mut self) {
        self.i_steiner |= STEINER_BIT;
    }
}

/// edges per block in the hole-bridge block-bbox index (see `Arena::build_block_index`)
const BLOCK_SIZE: i32 = 16;
const BLOCK_INDEX_MIN_NODES: usize = 256;

/// Arena of nodes backing the doubly linked list(s) used by the algorithm.
struct Arena {
    nodes: Vec<Node>,
    /// set by `filter_points` whenever it removes at least one node; read by
    /// `earcut_linked`'s stall handler to decide whether another clip pass is worth
    /// attempting before the costlier stages (mirrors the module-level flag in earcut.js).
    filtered_out: bool,

    /// Block-bbox index for `find_hole_bridge` (ported from earcut.hpp issue #183 fix): one
    /// `[min_x, min_y, max_x, max_y]` bbox per `BLOCK_SIZE` consecutive ring edges, so the
    /// leftward-ray scan can skip whole blocks in O(1) instead of walking the whole merged
    /// ring. Grown append-only: the outer ring seeds it, then each merged hole appends its
    /// own blocks. Buffers are reused/grown across calls.
    ///
    /// Each block owns an explicit list of node handles (`block_nodes[block_start[b]..
    /// block_start[b + 1]]`) rather than a ring range, because merging a hole splices a whole
    /// new ring run into the middle of an already-indexed block's range — so a range-based
    /// block would have to re-walk those nodes (which its own blocks already cover) on every
    /// later query. With explicit lists each block stays bounded by `BLOCK_SIZE`.
    ///
    /// `filter_points` only drops collinear/coincident points, so a stale bbox stays a
    /// conservative superset of its live edges (never a false skip); the scan skips dead
    /// nodes (`p.prev.next != p`). Blocks are scanned in append (not ring) order, so the
    /// chosen bridge can differ from an un-indexed scan — a different but equally valid
    /// result.
    block_bbox: Vec<f64>,
    block_nodes: Vec<u32>,
    block_start: Vec<u32>,
    num_blocks: usize,
    /// true only while `eliminate_holes` merges holes, so `remove_node` keeps the block
    /// index live (via `grow_block`)
    index_active: bool,

    /// scratch buffers for `index_curve`/`sort_by_z`, reused across (possibly recursive,
    /// via `split_earcut`) calls to avoid a heap allocation each time.
    z_order_buf: Vec<(i32, u32)>,
    z_order_scratch: Vec<(i32, u32)>,
}

impl Arena {
    fn new(capacity: usize) -> Self {
        Arena {
            nodes: Vec::with_capacity(capacity),
            filtered_out: false,
            block_bbox: Vec::new(),
            block_nodes: Vec::new(),
            block_start: Vec::new(),
            num_blocks: 0,
            index_active: false,
            z_order_buf: Vec::new(),
            z_order_scratch: Vec::new(),
        }
    }

    /// Clear all per-triangulation state while keeping every buffer's allocation, so a
    /// reused `Arena` (see `Earcut`) never reallocates once its buffers have grown to fit
    /// the largest input seen so far.
    fn reset(&mut self, capacity_hint: usize) {
        self.nodes.clear();
        self.nodes.reserve(capacity_hint);
        self.filtered_out = false;
        self.num_blocks = 0;
        self.index_active = false;
    }

    /// `Node` is kept as small as possible (40 bytes, 8-byte aligned) so large rings fit more
    /// nodes per cache line. Node "handles" (`prev`/`next`/`prev_z`/`next_z`, and every `u32`
    /// returned by `create_node`/`insert_node`) are **byte offsets** into `Arena::nodes`'s
    /// backing storage, not element indices: `create_node` computes `index * NODE_SIZE` once
    /// when a node is created, and `Arena::get`/`get_mut` dereference via `byte_add` (a plain
    /// pointer add) instead of slice indexing, which would redo `idx * NODE_SIZE` and a
    /// bounds check on every single access. Both sit on the load-to-load dependency chain of
    /// a ring walk, so on small polygons — where that chain *is* the runtime — plain indexing
    /// measures ~25% slower. This
    /// mirrors georust/earcut's `NodeOffset` design.
    const NODE_SIZE: usize = std::mem::size_of::<Node>();

    #[inline]
    fn create_node(&mut self, i: u32, x: f64, y: f64) -> u32 {
        let index = self.nodes.len();
        self.nodes.push(Node::new(i, x, y));
        // guard against a byte offset colliding with the `NULL` sentinel or overflowing `u32`;
        // in practice this would require ~100M+ vertices, already impractical for other reasons
        assert!(
            (index + 1) * Self::NODE_SIZE < NULL as usize,
            "too many nodes for a u32 byte-offset arena"
        );
        (index * Self::NODE_SIZE) as u32
    }

    #[inline(always)]
    fn get(&self, off: u32) -> &Node {
        debug_assert!(
            (off as usize).is_multiple_of(Self::NODE_SIZE),
            "misaligned node offset"
        );
        debug_assert!(
            (off as usize / Self::NODE_SIZE) < self.nodes.len(),
            "arena offset out of bounds"
        );
        // Safety: every offset handled by this arena is either `NULL` (checked by callers
        // before dereferencing) or was returned by `create_node`/`insert_node` as
        // `index * NODE_SIZE` for some `index < self.nodes.len()` at the time of creation;
        // `nodes` only ever grows (nodes are unlinked from the doubly linked list, not
        // deallocated), so any non-`NULL` offset passed here always lands on a live element.
        unsafe { &*self.nodes.as_ptr().byte_add(off as usize) }
    }

    #[inline(always)]
    fn get_mut(&mut self, off: u32) -> &mut Node {
        debug_assert!(
            (off as usize).is_multiple_of(Self::NODE_SIZE),
            "misaligned node offset"
        );
        debug_assert!(
            (off as usize / Self::NODE_SIZE) < self.nodes.len(),
            "arena offset out of bounds"
        );
        // Safety: see `get`.
        unsafe { &mut *self.nodes.as_mut_ptr().byte_add(off as usize) }
    }

    /// create a node and optionally link it with previous one (in a circular doubly linked list)
    fn insert_node(&mut self, i: u32, x: f64, y: f64, last: u32) -> u32 {
        let p = self.create_node(i, x, y);

        if last == NULL {
            self.get_mut(p).prev = p;
            self.get_mut(p).next = p;
        } else {
            let last_next = self.get(last).next;
            self.get_mut(p).next = last_next;
            self.get_mut(p).prev = last;
            self.get_mut(last_next).prev = p;
            self.get_mut(last).next = p;
        }
        p
    }

    fn remove_node(&mut self, p: u32) {
        let (prev, next, prev_z, next_z) = {
            let n = self.get(p);
            (n.prev, n.next, n.prev_z, n.next_z)
        };
        self.get_mut(next).prev = prev;
        self.get_mut(prev).next = next;

        if prev_z != NULL {
            self.get_mut(prev_z).next_z = next_z;
        }
        if next_z != NULL {
            self.get_mut(next_z).prev_z = prev_z;
        }

        // keep the hole-bridge index's block bboxes covering the healed prev->next edge
        if self.index_active {
            self.grow_block(prev, next);
        }
    }

    /// Block-bbox index buffers: size once from the input upper bound and reuse across calls.
    fn build_block_index(&mut self, max_nodes: usize, num_holes: usize) {
        // upper bound: every input node indexed once, +2 bridge nodes per hole, plus a
        // partial trailing block per appended segment (outer ring + one per hole)
        let block_size = BLOCK_SIZE as usize;
        let max_blocks = (max_nodes + 2 * num_holes).div_ceil(block_size) + num_holes + 2;
        if self.block_bbox.len() < max_blocks * 4 {
            self.block_bbox.resize(max_blocks * 4, 0.0);
        }
        self.num_blocks = 0;
        self.block_nodes.clear();
        self.block_start.clear();
        self.block_start.push(0);
    }

    /// index the ring run `head..stop` (exclusive) as `ceil(len / BLOCK_SIZE)` blocks; `head
    /// == stop` means the whole ring. Each block records the nodes it owns and a bbox
    /// covering both endpoints of every edge those nodes start.
    fn index_segment(&mut self, head: u32, stop: u32) {
        let mut p = head;
        loop {
            let b = self.num_blocks;
            self.num_blocks += 1;
            let mut b_min_x = f64::INFINITY;
            let mut b_min_y = f64::INFINITY;
            let mut b_max_x = f64::NEG_INFINITY;
            let mut b_max_y = f64::NEG_INFINITY;
            let mut k: i32 = 0;
            loop {
                let c = self.get(p).next; // edge p->c; bbox must bound both endpoints
                // reuse z as the owning block during eliminate_holes (see grow_block)
                self.get_mut(p).z = b as i32;
                self.block_nodes.push(p);
                let (px, py) = (self.get(p).x, self.get(p).y);
                let (cx, cy) = (self.get(c).x, self.get(c).y);
                b_min_x = b_min_x.min(px).min(cx);
                b_max_x = b_max_x.max(px).max(cx);
                b_min_y = b_min_y.min(py).min(cy);
                b_max_y = b_max_y.max(py).max(cy);
                p = c;
                k += 1;
                if k >= BLOCK_SIZE || p == stop {
                    break;
                }
            }
            self.block_start.push(self.block_nodes.len() as u32);
            let g = b * 4;
            self.block_bbox[g] = b_min_x;
            self.block_bbox[g + 1] = b_min_y;
            self.block_bbox[g + 2] = b_max_x;
            self.block_bbox[g + 3] = b_max_y;
            if p == stop {
                break;
            }
        }
    }

    /// when `filter_points` heals an edge head->tail (removing the collinear node between
    /// them), the healed edge can extend past head's frozen block bbox if its old far
    /// endpoint lived in another block; grow head's block bbox to cover tail so the
    /// leftward-ray prune can't false-skip it.
    fn grow_block(&mut self, head: u32, tail: u32) {
        let g = self.get(head).z as usize * 4;
        let (tx, ty) = (self.get(tail).x, self.get(tail).y);
        if tx < self.block_bbox[g] {
            self.block_bbox[g] = tx;
        }
        if ty < self.block_bbox[g + 1] {
            self.block_bbox[g + 1] = ty;
        }
        if tx > self.block_bbox[g + 2] {
            self.block_bbox[g + 2] = tx;
        }
        if ty > self.block_bbox[g + 3] {
            self.block_bbox[g + 3] = ty;
        }
    }

    /// the block's head node can be removed by `filter_points` during merges; advance it to
    /// the next live node so the walk doesn't start on (and immediately terminate at) a dead
    /// node. For the single full-ring seed block (`head == stop`) the same forward advance
    /// keeps them equal, so the loop still laps the whole ring instead of collapsing to an
    /// empty walk.
    #[inline]
    fn block_range(&self, b: usize) -> std::ops::Range<usize> {
        self.block_start[b] as usize..self.block_start[b + 1] as usize
    }

    /// a node's coordinates
    #[inline]
    fn xy(&self, p: u32) -> (f64, f64) {
        let n = self.get(p);
        (n.x, n.y)
    }

    /// a node is live while its ring neighbours still point back at it; `remove_node` leaves
    /// removed nodes' own links intact, so stale index entries are detected this way
    #[inline]
    fn is_live(&self, p: u32, prev: u32) -> bool {
        self.get(prev).next == p
    }

    #[inline]
    fn equals(&self, p1: u32, p2: u32) -> bool {
        let a = self.get(p1);
        let b = self.get(p2);
        a.x == b.x && a.y == b.y
    }

    /// `area`, for a middle vertex whose coordinates the caller already has
    #[inline]
    fn area_at(&self, p: u32, qx: f64, qy: f64, r: u32) -> f64 {
        let p = self.get(p);
        let r = self.get(r);
        area_xy(p.x, p.y, qx, qy, r.x, r.y)
    }

    /// signed area of a triangle
    #[inline]
    fn area(&self, p: u32, q: u32, r: u32) -> f64 {
        let p = self.get(p);
        let q = self.get(q);
        let r = self.get(r);
        area_xy(p.x, p.y, q.x, q.y, r.x, r.y)
    }
}

impl Default for Arena {
    fn default() -> Self {
        Arena::new(0)
    }
}

/// Triangulate a polygon given as a flat array of vertex coordinates.
///
/// `data` is a flat array of vertex coordinates (`dim` numbers per vertex, x/y first).
/// `hole_indices` lists the vertex index (not coordinate index) at which each hole ring
/// starts in `data`; the first ring (indices `0..hole_indices[0]`, or all of `data` if
/// there are no holes) is the outer ring.
///
/// Returns the triangulation as flat triplets of vertex indices into `data`.
///
/// This allocates a fresh internal arena on every call; for repeated triangulations,
/// [`Earcut`] reuses its arena's allocations across calls instead.
pub fn earcut<N: EarcutIndex>(data: &[f64], hole_indices: &[usize], dim: usize) -> Vec<N> {
    let mut triangles = Vec::new();
    earcut_into(data, hole_indices, dim, &mut triangles);
    triangles
}

/// Same as [`earcut`], but appends into a caller-provided `Vec` (which is cleared first).
/// Useful to reuse the output allocation across repeated triangulations; the internal arena
/// is still freshly allocated each call — use [`Earcut`] to reuse both.
pub fn earcut_into<N: EarcutIndex>(
    data: &[f64],
    hole_indices: &[usize],
    dim: usize,
    triangles: &mut Vec<N>,
) {
    let mut arena = Arena::new(0);
    earcut_impl(&mut arena, data, hole_indices, dim, triangles);
}

/// A reusable triangulator: keeps its internal arena's allocations (nodes, hole-bridge
/// block index, z-order sort scratch buffers) across calls, so repeated triangulations
/// (even of different, unrelated polygons) only allocate while growing to the largest
/// input seen so far, instead of once per call like the free [`earcut`]/[`earcut_into`]
/// functions.
#[derive(Default)]
pub struct Earcut {
    arena: Arena,
}

impl Earcut {
    /// Create a new, empty reusable triangulator.
    pub fn new() -> Self {
        Self {
            arena: Arena::new(0),
        }
    }

    /// Triangulate `data`/`hole_indices`/`dim` (see [`earcut`]), appending into
    /// `triangles` (cleared first), reusing this `Earcut`'s internal buffers.
    pub fn earcut_into<N: EarcutIndex>(
        &mut self,
        data: &[f64],
        hole_indices: &[usize],
        dim: usize,
        triangles: &mut Vec<N>,
    ) {
        earcut_impl(&mut self.arena, data, hole_indices, dim, triangles);
    }
}

/// Shared triangulation body driving `arena`; used by both the one-shot free functions
/// (fresh `Arena` each call) and `Earcut` (persistent `Arena` reused across calls).
fn earcut_impl<N: EarcutIndex>(
    arena: &mut Arena,
    data: &[f64],
    hole_indices: &[usize],
    dim: usize,
    triangles: &mut Vec<N>,
) {
    triangles.clear();

    if data.is_empty() || dim == 0 {
        return;
    }

    // A simple polygon with `n` vertices produces at most `n - 2` triangles;
    // reserving up front avoids repeated reallocation as triangles are pushed.
    let n_estimate = data.len() / dim;
    triangles.reserve(n_estimate.saturating_sub(2) * 3);

    let has_holes = !hole_indices.is_empty();
    let outer_len = if has_holes {
        hole_indices[0] * dim
    } else {
        data.len()
    };

    arena.reset(n_estimate + 8);

    let outer_node = linked_list(arena, data, 0, outer_len, dim, true);
    let mut outer_node = match outer_node {
        Some(n) => n,
        None => return,
    };

    if arena.get(outer_node).next == arena.get(outer_node).prev {
        return;
    }

    let mut min_x = 0.0;
    let mut min_y = 0.0;
    let mut inv_size = 0.0;

    if has_holes {
        outer_node = eliminate_holes(arena, data, hole_indices, outer_node, dim);
    }

    if data.len() > 80 * dim {
        min_x = data[0];
        min_y = data[1];
        let mut max_x = min_x;
        let mut max_y = min_y;

        let mut i = dim;
        while i < outer_len {
            let x = data[i];
            let y = data[i + 1];
            if x < min_x {
                min_x = x;
            }
            if y < min_y {
                min_y = y;
            }
            if x > max_x {
                max_x = x;
            }
            if y > max_y {
                max_y = y;
            }
            i += dim;
        }

        inv_size = f64::max(max_x - min_x, max_y - min_y);
        inv_size = if inv_size != 0.0 {
            32767.0 / inv_size
        } else {
            0.0
        };
    }

    earcut_linked(arena, outer_node, triangles, min_x, min_y, inv_size);
}

/// create a circular doubly linked list from polygon points in the specified winding order
fn linked_list(
    arena: &mut Arena,
    data: &[f64],
    start: usize,
    end: usize,
    dim: usize,
    clockwise: bool,
) -> Option<u32> {
    let mut last = NULL;

    if clockwise == (signed_area(data, start, end, dim) > 0.0) {
        let mut i = start;
        while i < end {
            last = arena.insert_node((i / dim) as u32, data[i], data[i + 1], last);
            i += dim;
        }
    } else {
        let mut i = end;
        while i > start {
            i -= dim;
            last = arena.insert_node((i / dim) as u32, data[i], data[i + 1], last);
        }
    }

    if last != NULL {
        let next = arena.get(last).next;
        if arena.equals(last, next) {
            arena.remove_node(last);
            last = arena.get(last).next;
        }
    }

    if last == NULL { None } else { Some(last) }
}

/// Remove collinear or coincident points; removability depends only on a node's immediate
/// neighbors, so we sweep forward and re-check the predecessor after each removal. With
/// `end == None` we sweep the whole ring, lapping until nothing is removable (the fixpoint
/// the clipper needs). With an explicit `end` we heal only the dirty window around a
/// bridge/diagonal cut, stopping at `end` rather than lapping.
fn filter_points(arena: &mut Arena, start: u32, end: Option<u32>) -> u32 {
    let full = end.is_none();
    let mut end = end.unwrap_or(start);

    let mut p = start;
    let mut again;
    loop {
        again = false;
        let (px, py, p_next, p_prev, steiner) = {
            let n = arena.get(p);
            (n.x, n.y, n.next, n.prev, n.is_steiner())
        };
        let degenerate = p != p_next && !steiner && {
            // coincident with the next point, or collinear with both neighbours
            let (nx, ny) = arena.xy(p_next);
            (px == nx && py == ny) || arena.area_at(p_prev, px, py, p_next) == 0.0
        };
        if degenerate {
            if full || p == end {
                end = p_prev;
            }
            arena.filtered_out = true;
            arena.remove_node(p);
            p = p_prev;
            again = true;
        } else if full || p != end {
            p = p_next;
            again = !full;
        }
        if !(again || p != end) {
            break;
        }
    }

    end
}

/// main ear slicing loop which triangulates a polygon (given as a linked list)
fn earcut_linked<N: EarcutIndex>(
    arena: &mut Arena,
    mut ear: u32,
    triangles: &mut Vec<N>,
    min_x: f64,
    min_y: f64,
    inv_size: f64,
) {
    // interlink polygon nodes in z-order
    if inv_size != 0.0 {
        index_curve(arena, ear, min_x, min_y, inv_size);
    }

    let mut stop = ear;
    let mut cured = false;

    // iterate through ears, slicing them one by one
    loop {
        let (prev, next, prev_z, next_z, ex, ey) = {
            let n = arena.get(ear);
            (n.prev, n.next, n.prev_z, n.next_z, n.x, n.y)
        };
        if prev == next {
            break;
        }

        if arena.area_at(prev, ex, ey, next) < 0.0
            && (if inv_size != 0.0 {
                is_ear_hashed(
                    arena, prev, next, prev_z, next_z, ex, ey, min_x, min_y, inv_size,
                )
            } else {
                is_ear(arena, ear)
            })
        {
            let pi = arena.get(prev).i() as usize;
            let ei = arena.get(ear).i() as usize;
            let ni = arena.get(next).i() as usize;
            triangles.push(N::from_usize(pi));
            triangles.push(N::from_usize(ei));
            triangles.push(N::from_usize(ni));

            arena.remove_node(ear);

            ear = next;
            stop = next;
            continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if ear == stop {
            // try filtering collinear/coincident points and slicing again — repeat as long
            // as filtering actually removes nodes, since each removal can expose new ears
            arena.filtered_out = false;
            ear = filter_points(arena, ear, None);
            if arena.filtered_out {
                stop = ear;
                continue;
            }

            // filtering is exhausted: cure small local self-intersections once, then retry
            if !cured {
                ear = cure_local_intersections(arena, ear, triangles);
                stop = ear;
                cured = true;
                continue;
            }

            // as a last resort, try splitting the remaining polygon into two
            split_earcut(arena, ear, triangles, min_x, min_y, inv_size);
            break;
        }
    }
}

/// check whether a polygon node forms a valid ear with adjacent nodes
fn is_ear(arena: &Arena, ear: u32) -> bool {
    let a = arena.get(ear).prev;
    let b = ear;
    let c = arena.get(ear).next;

    let (ax, ay) = arena.xy(a);
    let (bx, by) = arena.xy(b);
    let (cx, cy) = arena.xy(c);

    let x0 = ax.min(bx).min(cx);
    let y0 = ay.min(by).min(cy);
    let x1 = ax.max(bx).max(cx);
    let y1 = ay.max(by).max(cy);

    let mut p = arena.get(c).next;
    while p != a {
        let (px, py, p_prev, p_next) = {
            let n = arena.get(p);
            (n.x, n.y, n.prev, n.next)
        };
        if px >= x0
            && px <= x1
            && py >= y0
            && py <= y1
            && !(ax == px && ay == py)
            && point_in_triangle(ax, ay, bx, by, cx, cy, px, py)
            && arena.area_at(p_prev, px, py, p_next) >= 0.0
        {
            return false;
        }
        p = p_next;
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn is_ear_hashed(
    arena: &Arena,
    a: u32,
    c: u32,
    prev_z: u32,
    next_z: u32,
    bx: f64,
    by: f64,
    min_x: f64,
    min_y: f64,
    inv_size: f64,
) -> bool {
    let (ax, ay) = (arena.get(a).x, arena.get(a).y);
    let (cx, cy) = (arena.get(c).x, arena.get(c).y);

    let x0 = ax.min(bx).min(cx);
    let y0 = ay.min(by).min(cy);
    let x1 = ax.max(bx).max(cx);
    let y1 = ay.max(by).max(cy);

    let min_z = z_order(x0, y0, min_x, min_y, inv_size);
    let max_z = z_order(x1, y1, min_x, min_y, inv_size);

    let mut p = prev_z;
    while p != NULL {
        let pn = arena.get(p);
        if pn.z < min_z {
            break;
        }
        let (px, py) = (pn.x, pn.y);
        if p != c
            && px >= x0
            && px <= x1
            && py >= y0
            && py <= y1
            && !(ax == px && ay == py)
            && point_in_triangle(ax, ay, bx, by, cx, cy, px, py)
            && arena.area_at(pn.prev, px, py, pn.next) >= 0.0
        {
            return false;
        }
        p = pn.prev_z;
    }

    let mut n = next_z;
    while n != NULL {
        let nn = arena.get(n);
        if nn.z > max_z {
            break;
        }
        let (nx, ny) = (nn.x, nn.y);
        if n != c
            && nx >= x0
            && nx <= x1
            && ny >= y0
            && ny <= y1
            && !(ax == nx && ay == ny)
            && point_in_triangle(ax, ay, bx, by, cx, cy, nx, ny)
            && arena.area_at(nn.prev, nx, ny, nn.next) >= 0.0
        {
            return false;
        }
        n = nn.next_z;
    }

    true
}

/// go through all polygon nodes and cure small local self-intersections
fn cure_local_intersections<N: EarcutIndex>(
    arena: &mut Arena,
    start: u32,
    triangles: &mut Vec<N>,
) -> u32 {
    let mut p = start;
    let mut start = start;
    let mut cured = false;
    loop {
        let a = arena.get(p).prev;
        let b = arena.get(arena.get(p).next).next;

        if intersects(arena, a, p, arena.get(p).next, b, false)
            && locally_inside(arena, a, b)
            && locally_inside(arena, b, a)
        {
            triangles.push(N::from_usize(arena.get(a).i() as usize));
            triangles.push(N::from_usize(arena.get(p).i() as usize));
            triangles.push(N::from_usize(arena.get(b).i() as usize));

            let p_next = arena.get(p).next;
            arena.remove_node(p);
            arena.remove_node(p_next);

            p = b;
            start = b;
            cured = true;
        }
        p = arena.get(p).next;
        if p == start {
            break;
        }
    }

    if cured {
        filter_points(arena, p, None)
    } else {
        p
    }
}

/// try splitting polygon into two and triangulate them independently
fn split_earcut<N: EarcutIndex>(
    arena: &mut Arena,
    start: u32,
    triangles: &mut Vec<N>,
    min_x: f64,
    min_y: f64,
    inv_size: f64,
) {
    let mut a = start;
    loop {
        let mut b = arena.get(arena.get(a).next).next;
        while b != arena.get(a).prev {
            if arena.get(a).i() != arena.get(b).i() && is_valid_diagonal(arena, a, b) {
                let c = split_polygon(arena, a, b);

                let a_next = arena.get(a).next;
                let new_a = filter_points(arena, a, Some(a_next));
                let c_next = arena.get(c).next;
                let new_c = filter_points(arena, c, Some(c_next));

                earcut_linked(arena, new_a, triangles, min_x, min_y, inv_size);
                earcut_linked(arena, new_c, triangles, min_x, min_y, inv_size);
                return;
            }
            b = arena.get(b).next;
        }
        a = arena.get(a).next;
        if a == start {
            break;
        }
    }
}

/// link every hole into the outer loop, producing a single-ring polygon without holes
fn eliminate_holes(
    arena: &mut Arena,
    data: &[f64],
    hole_indices: &[usize],
    mut outer_node: u32,
    dim: usize,
) -> u32 {
    let mut queue: Vec<u32> = Vec::with_capacity(hole_indices.len());

    let len = hole_indices.len();
    for i in 0..len {
        let start = hole_indices[i] * dim;
        let end = if i < len - 1 {
            hole_indices[i + 1] * dim
        } else {
            data.len()
        };
        if let Some(list) = linked_list(arena, data, start, end, dim, false) {
            if list == arena.get(list).next {
                arena.get_mut(list).set_steiner();
            }
            queue.push(get_leftmost(arena, list));
        }
    }

    queue.sort_unstable_by(|&a, &b| compare_xy_slope(arena, a, b));

    // The block index amortizes its setup cost on larger rings. Small polygons are faster with
    // the original linear bridge scan because they do not have enough edges for bbox pruning to
    // pay for building and checking the index.
    let use_block_index = data.len() / dim >= BLOCK_INDEX_MIN_NODES;
    if use_block_index {
        arena.build_block_index(data.len() / dim, queue.len());
        arena.index_segment(outer_node, outer_node);
        arena.index_active = true;
    }

    // process holes from left to right
    for hole in queue {
        outer_node = eliminate_hole(arena, hole, outer_node);
    }
    arena.index_active = false;

    // collapse collinear/coincident points across the whole merged ring once before clipping
    filter_points(arena, outer_node, None)
}

fn compare_xy_slope(arena: &Arena, a: u32, b: u32) -> std::cmp::Ordering {
    let (an, bn) = (arena.get(a), arena.get(b));
    if an.x != bn.x {
        return an.x.partial_cmp(&bn.x).unwrap();
    }
    if an.y != bn.y {
        return an.y.partial_cmp(&bn.y).unwrap();
    }
    let a_next = arena.get(an.next);
    let b_next = arena.get(bn.next);
    let slope_a = (a_next.y - an.y) / (a_next.x - an.x);
    let slope_b = (b_next.y - bn.y) / (b_next.x - bn.x);
    slope_a
        .partial_cmp(&slope_b)
        .unwrap_or(std::cmp::Ordering::Equal)
}

/// find a bridge between vertices that connects hole with an outer ring and link it
fn eliminate_hole(arena: &mut Arena, hole: u32, outer_node: u32) -> u32 {
    let bridge = find_hole_bridge(arena, hole, outer_node);
    let bridge = match bridge {
        Some(b) => b,
        None => return outer_node,
    };

    let bridge_reverse = split_polygon(arena, bridge, hole);

    if arena.index_active {
        // index the merged-in segment before filtering: in ring order the splice runs
        // bridge -> hole -> bridge_reverse -> bridge2 -> (bridge's old next), covering the
        // hole's edges and both new slit edges. filter_points below only drops
        // collinear/coincident points, so these bboxes stay valid (conservative) supersets.
        let bridge2 = arena.get(bridge_reverse).next;
        let bridge2_next = arena.get(bridge2).next;
        arena.index_segment(bridge, bridge2_next);
    }

    // heal collinear/coincident points around the two new slit edges
    let bridge_reverse_next = arena.get(bridge_reverse).next;
    filter_points(arena, bridge_reverse, Some(bridge_reverse_next));
    let bridge_next = arena.get(bridge).next;
    filter_points(arena, bridge, Some(bridge_next))
}

/// David Eberly's algorithm for finding a bridge between hole and outer polygon, using the
/// block-bbox index built by `eliminate_holes` to skip whole runs of ring edges that can't
/// possibly beat the current best crossing.
fn find_hole_bridge(arena: &mut Arena, hole: u32, outer_node: u32) -> Option<u32> {
    if arena.num_blocks == 0 {
        return find_hole_bridge_linear(arena, hole, outer_node);
    }

    let p = outer_node;
    let (hx, hy) = arena.xy(hole);
    let mut qx = f64::NEG_INFINITY;
    let mut m: Option<u32> = None;

    // find a segment intersected by a ray from the hole's leftmost vertex to the left;
    // segment's endpoint with lesser x will be potential connection vertex, unless they
    // intersect at a vertex, then choose the vertex
    if arena.equals(hole, p) {
        return Some(p);
    }

    // scan blocks; skip any whose bbox can't hold a crossing that beats qx and lies left of hx
    let mut b = 0;
    let mut g = 0;
    while b < arena.num_blocks {
        let (b_min_x, b_min_y, b_max_x, b_max_y) = (
            arena.block_bbox[g],
            arena.block_bbox[g + 1],
            arena.block_bbox[g + 2],
            arena.block_bbox[g + 3],
        );
        if hy < b_min_y || hy > b_max_y || b_min_x > hx || b_max_x <= qx {
            b += 1;
            g += 4;
            continue;
        }

        for i in arena.block_range(b) {
            let p = arena.block_nodes[i];
            let (px, py, p_prev, p_next) = {
                let n = arena.get(p);
                (n.x, n.y, n.prev, n.next)
            };
            // skip nodes removed by filter_points (stale in the index)
            if arena.is_live(p, p_prev) {
                let (nx_, ny_) = arena.xy(p_next);
                if hx == nx_ && hy == ny_ {
                    return Some(p_next);
                } else if hy <= py && hy >= ny_ && ny_ != py {
                    let x = px + (hy - py) * (nx_ - px) / (ny_ - py);
                    if x <= hx && x > qx {
                        qx = x;
                        m = Some(if px < nx_ { p } else { p_next });
                        if x == hx {
                            return m; // hole touches outer segment; pick leftmost endpoint
                        }
                    }
                }
            }
        }
        b += 1;
        g += 4;
    }

    let m0 = m?;

    // look for points inside the triangle of hole vertex, segment intersection and endpoint;
    // if there are no points found, we have a valid connection; otherwise choose the vertex
    // of the minimum angle with the ray as connection vertex
    let (mx, my) = arena.xy(m0);
    let t_min_y = hy.min(my); // the triangle's y span; x span is [mx, hx]
    let t_max_y = hy.max(my);
    let mut tan_min = f64::INFINITY;
    let mut m = m0;

    // scan the same blocks; skip any whose bbox can't overlap the triangle's
    // [mx,hx]x[t_min_y,t_max_y] box
    let mut b = 0;
    let mut g = 0;
    while b < arena.num_blocks {
        let (b_min_x, b_min_y, b_max_x, b_max_y) = (
            arena.block_bbox[g],
            arena.block_bbox[g + 1],
            arena.block_bbox[g + 2],
            arena.block_bbox[g + 3],
        );
        if b_max_x < mx || b_min_x > hx || b_max_y < t_min_y || b_min_y > t_max_y {
            b += 1;
            g += 4;
            continue;
        }

        for i in arena.block_range(b) {
            let p = arena.block_nodes[i];
            let (px, py, p_prev, p_next) = {
                let n = arena.get(p);
                (n.x, n.y, n.prev, n.next)
            };
            if arena.is_live(p, p_prev) // skip nodes removed by filter_points
                && hx >= px
                && px >= mx
                && hx != px
                && point_in_triangle(
                    if hy < my { hx } else { qx },
                    hy,
                    mx,
                    my,
                    if hy < my { qx } else { hx },
                    hy,
                    px,
                    py,
                )
            {
                let tan = (hy - py).abs() / (hx - px); // tangential

                // if hole point sits on p's horizontal edge (T-junction touch): the bridge
                // runs along that edge — locally_inside rejects it as collinear, but it's valid
                let mx_ = arena.get(m).x;
                if (locally_inside(arena, p, hole) || {
                    let n = arena.get(p_next);
                    py == hy && n.y == hy && n.x > hx
                }) && (tan < tan_min
                    || (tan == tan_min
                        && (px > mx_ || (px == mx_ && sector_contains_sector(arena, m, p)))))
                {
                    m = p;
                    tan_min = tan;
                }
            }
        }
        b += 1;
        g += 4;
    }

    Some(m)
}

/// Linear bridge search for small polygons, where building and scanning the block index costs
/// more than the bbox checks can save.
fn find_hole_bridge_linear(arena: &Arena, hole: u32, outer_node: u32) -> Option<u32> {
    let (hx, hy) = arena.xy(hole);
    let mut qx = f64::NEG_INFINITY;
    let mut m: Option<u32> = None;

    if arena.equals(hole, outer_node) {
        return Some(outer_node);
    }

    let mut p = outer_node;
    loop {
        let node = arena.get(p);
        let p_next = node.next;
        let (px, py) = (node.x, node.y);
        let next = arena.get(p_next);
        let (nx, ny) = (next.x, next.y);
        if hx == nx && hy == ny {
            return Some(p_next);
        }
        if hy <= py && hy >= ny && ny != py {
            let x = px + (hy - py) * (nx - px) / (ny - py);
            if x <= hx && x > qx {
                qx = x;
                m = Some(if px < nx { p } else { p_next });
                if x == hx {
                    return m;
                }
            }
        }
        p = p_next;
        if p == outer_node {
            break;
        }
    }

    let m0 = m?;
    let (mx, my) = arena.xy(m0);
    let mut tan_min = f64::INFINITY;
    let mut m = m0;

    let mut p = outer_node;
    loop {
        let node = arena.get(p);
        let p_next = node.next;
        let (px, py) = (node.x, node.y);
        if hx >= px
            && px >= mx
            && hx != px
            && point_in_triangle(
                if hy < my { hx } else { qx },
                hy,
                mx,
                my,
                if hy < my { qx } else { hx },
                hy,
                px,
                py,
            )
        {
            let tan = (hy - py).abs() / (hx - px);
            let next = arena.get(p_next);
            let current_mx = arena.get(m).x;
            if (locally_inside(arena, p, hole) || (py == hy && next.y == hy && next.x > hx))
                && (tan < tan_min
                    || (tan == tan_min
                        && (px > current_mx
                            || (px == current_mx && sector_contains_sector(arena, m, p)))))
            {
                m = p;
                tan_min = tan;
            }
        }
        p = p_next;
        if p == outer_node {
            break;
        }
    }

    Some(m)
}

/// whether sector in vertex m contains sector in vertex p in the same coordinates
fn sector_contains_sector(arena: &Arena, m: u32, p: u32) -> bool {
    arena.area(arena.get(m).prev, m, arena.get(p).prev) < 0.0
        && arena.area(arena.get(p).next, m, arena.get(m).next) < 0.0
}

/// interlink polygon nodes in z-order: collect into a vec, sort by z, relink
fn index_curve(arena: &mut Arena, start: u32, min_x: f64, min_y: f64, inv_size: f64) {
    arena.z_order_buf.clear();
    let mut p = start;
    loop {
        let (x, y) = (arena.get(p).x, arena.get(p).y);
        let z = z_order(x, y, min_x, min_y, inv_size);
        arena.get_mut(p).z = z;
        arena.z_order_buf.push((z, p));
        p = arena.get(p).next;
        if p == start {
            break;
        }
    }

    // reuse the arena's scratch Vecs across calls (avoids a heap allocation per
    // `index_curve` invocation, which recurses through `split_earcut`)
    let mut order = std::mem::take(&mut arena.z_order_buf);
    let mut scratch = std::mem::take(&mut arena.z_order_scratch);
    sort_by_z(&mut order, &mut scratch);

    let mut prev: u32 = NULL;
    for &(_, node) in order.iter() {
        arena.get_mut(node).prev_z = prev;
        if prev != NULL {
            arena.get_mut(prev).next_z = node;
        }
        prev = node;
    }
    arena.get_mut(prev).next_z = NULL;

    arena.z_order_buf = order;
    arena.z_order_scratch = scratch;
}

/// Sort a z-order queue. Small queues use `sort_unstable_by_key` (introsort); large ones
/// (>= `RADIX_MIN_LEN`) use four stable byte-wise LSD radix passes, which is O(n) instead of
/// O(n log n) and a clear win once n is large enough to amortize the counting passes (this
/// is where `index_curve` spends most of its time on big star/hole-grid inputs).
fn sort_by_z(order: &mut [(i32, u32)], scratch: &mut Vec<(i32, u32)>) {
    const RADIX_MIN_LEN: usize = 512;
    let n = order.len();
    if n < RADIX_MIN_LEN {
        order.sort_unstable_by_key(|&(z, _)| z);
        return;
    }

    #[inline(always)]
    fn scatter(src: &[(i32, u32)], dst: &mut [(i32, u32)], shift: u32, offsets: &mut [u32; 256]) {
        for &item in src {
            let bucket = (item.0 as u32 >> shift & 0xff) as usize;
            dst[offsets[bucket] as usize] = item;
            offsets[bucket] += 1;
        }
    }

    let mut counts = [[0u32; 256]; 4];
    for &(z, _) in order.iter() {
        let key = z as u32;
        counts[0][(key & 0xff) as usize] += 1;
        counts[1][(key >> 8 & 0xff) as usize] += 1;
        counts[2][(key >> 16 & 0xff) as usize] += 1;
        counts[3][(key >> 24 & 0xff) as usize] += 1;
    }

    scratch.resize(n, order[0]);
    let mut in_order = true;
    for (pass, counts) in counts.iter().enumerate() {
        // a pass where every item lands in the same bucket can't reorder anything: skip it
        if counts.iter().any(|&count| count as usize == n) {
            continue;
        }
        let mut offsets = [0u32; 256];
        let mut sum = 0;
        for (offset, &count) in offsets.iter_mut().zip(counts.iter()) {
            *offset = sum;
            sum += count;
        }
        let shift = pass as u32 * 8;
        if in_order {
            scatter(order, scratch, shift, &mut offsets);
        } else {
            scatter(scratch, order, shift, &mut offsets);
        }
        in_order = !in_order;
    }
    if !in_order {
        order.copy_from_slice(scratch);
    }
}

/// z-order of a point given coords and inverse of the longer side of data bbox
fn z_order(x: f64, y: f64, min_x: f64, min_y: f64, inv_size: f64) -> i32 {
    let mut x = ((x - min_x) * inv_size) as i32;
    let mut y = ((y - min_y) * inv_size) as i32;

    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    x | (y << 1)
}

/// find the leftmost node of a polygon ring
fn get_leftmost(arena: &Arena, start: u32) -> u32 {
    let mut p = start;
    let mut leftmost = start;
    loop {
        let pn = arena.get(p);
        let ln = arena.get(leftmost);
        if pn.x < ln.x || (pn.x == ln.x && pn.y < ln.y) {
            leftmost = p;
        }
        p = arena.get(p).next;
        if p == start {
            break;
        }
    }
    leftmost
}

/// signed area of the triangle `(p, q, r)`, from raw coordinates
#[inline]
fn area_xy(px: f64, py: f64, qx: f64, qy: f64, rx: f64, ry: f64) -> f64 {
    (qy - py) * (rx - qx) - (qx - px) * (ry - qy)
}

/// check if a point lies within a convex triangle
#[allow(clippy::too_many_arguments)]
#[inline]
fn point_in_triangle(
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
    px: f64,
    py: f64,
) -> bool {
    (cx - px) * (ay - py) >= (ax - px) * (cy - py)
        && (ax - px) * (by - py) >= (bx - px) * (ay - py)
        && (bx - px) * (cy - py) >= (cx - px) * (by - py)
}

/// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
fn is_valid_diagonal(arena: &Arena, a: u32, b: u32) -> bool {
    let zero_length = arena.equals(a, b)
        && arena.area(arena.get(a).prev, a, arena.get(a).next) > 0.0
        && arena.area(arena.get(b).prev, b, arena.get(b).next) > 0.0;

    arena.get(arena.get(a).next).i() != arena.get(b).i()
        && (zero_length
            || (locally_inside(arena, a, b)
                && locally_inside(arena, b, a)
                && (arena.area(arena.get(a).prev, a, arena.get(b).prev) != 0.0
                    || arena.area(a, arena.get(b).prev, b) != 0.0)))
        && !intersects_polygon(arena, a, b)
        && (zero_length || middle_inside(arena, a, b))
}

/// check if two segments intersect; by default includes collinear boundary touches
fn intersects(arena: &Arena, p1: u32, q1: u32, p2: u32, q2: u32, include_boundary: bool) -> bool {
    let (p1x, p1y) = arena.xy(p1);
    let (q1x, q1y) = arena.xy(q1);
    let (p2x, p2y) = arena.xy(p2);
    let (q2x, q2y) = arena.xy(q2);

    let o1 = area_xy(p1x, p1y, q1x, q1y, p2x, p2y);
    let o2 = area_xy(p1x, p1y, q1x, q1y, q2x, q2y);
    let o3 = area_xy(p2x, p2y, q2x, q2y, p1x, p1y);
    let o4 = area_xy(p2x, p2y, q2x, q2y, q1x, q1y);

    if ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0))
        && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))
    {
        return true;
    }

    include_boundary
        && ((o1 == 0.0 && on_segment(p1x, p1y, p2x, p2y, q1x, q1y))
            || (o2 == 0.0 && on_segment(p1x, p1y, q2x, q2y, q1x, q1y))
            || (o3 == 0.0 && on_segment(p2x, p2y, p1x, p1y, q2x, q2y))
            || (o4 == 0.0 && on_segment(p2x, p2y, q1x, q1y, q2x, q2y)))
}

/// for collinear points p, q, r, check if point q lies on segment pr
#[inline]
fn on_segment(px: f64, py: f64, qx: f64, qy: f64, rx: f64, ry: f64) -> bool {
    qx <= px.max(rx) && qx >= px.min(rx) && qy <= py.max(ry) && qy >= py.min(ry)
}

/// check if a polygon diagonal intersects any polygon segments
fn intersects_polygon(arena: &Arena, a: u32, b: u32) -> bool {
    let (an, bn) = (arena.get(a), arena.get(b));
    let min_x = an.x.min(bn.x);
    let max_x = an.x.max(bn.x);
    let min_y = an.y.min(bn.y);
    let max_y = an.y.max(bn.y);

    let a_i = an.i();
    let b_i = bn.i();

    let mut p = a;
    let (mut px, mut py, mut p_i) = {
        let n = arena.get(p);
        (n.x, n.y, n.i())
    };
    loop {
        let n = arena.get(p).next;
        let (nx, ny, n_i) = {
            let v = arena.get(n);
            (v.x, v.y, v.i())
        };
        if !((px > max_x && nx > max_x)
            || (px < min_x && nx < min_x)
            || (py > max_y && ny > max_y)
            || (py < min_y && ny < min_y))
            && p_i != a_i
            && n_i != a_i
            && p_i != b_i
            && n_i != b_i
            && intersects(arena, p, n, a, b, true)
        {
            return true;
        }
        p = n;
        if p == a {
            break;
        }
        (px, py, p_i) = (nx, ny, n_i);
    }

    false
}

/// check if a polygon diagonal is locally inside the polygon
fn locally_inside(arena: &Arena, a: u32, b: u32) -> bool {
    let (a_prev, a_next) = {
        let n = arena.get(a);
        (n.prev, n.next)
    };
    let (ax, ay) = arena.xy(a);
    let (vx, vy) = arena.xy(a_prev);
    let (wx, wy) = arena.xy(a_next);
    let (bx, by) = arena.xy(b);
    if area_xy(vx, vy, ax, ay, wx, wy) < 0.0 {
        area_xy(ax, ay, bx, by, wx, wy) >= 0.0 && area_xy(ax, ay, vx, vy, bx, by) >= 0.0
    } else {
        area_xy(ax, ay, bx, by, vx, vy) < 0.0 || area_xy(ax, ay, wx, wy, bx, by) < 0.0
    }
}

/// check if the middle point of a polygon diagonal is inside the polygon
fn middle_inside(arena: &Arena, a: u32, b: u32) -> bool {
    let (an, bn) = (arena.get(a), arena.get(b));
    let px = (an.x + bn.x) / 2.0;
    let py = (an.y + bn.y) / 2.0;

    let mut p = a;
    let (mut qx, mut qy) = arena.xy(p);
    let mut inside = false;
    loop {
        let n = arena.get(p).next;
        let (nx, ny) = arena.xy(n);
        if (qy > py) != (ny > py) && px < (nx - qx) * (py - qy) / (ny - qy) + qx {
            inside = !inside;
        }
        p = n;
        if p == a {
            break;
        }
        (qx, qy) = (nx, ny);
    }

    inside
}

/// link two polygon vertices with a bridge; if the vertices belong to the same ring, it
/// splits the polygon into two; if one belongs to the outer ring and another to a hole, it
/// merges them into a single ring
fn split_polygon(arena: &mut Arena, a: u32, b: u32) -> u32 {
    let (a_i, a_x, a_y) = {
        let n = arena.get(a);
        (n.i(), n.x, n.y)
    };
    let (b_i, b_x, b_y) = {
        let n = arena.get(b);
        (n.i(), n.x, n.y)
    };

    let a2 = arena.create_node(a_i, a_x, a_y);
    let b2 = arena.create_node(b_i, b_x, b_y);

    let an = arena.get(a).next;
    let bp = arena.get(b).prev;

    arena.get_mut(a).next = b;
    arena.get_mut(b).prev = a;

    arena.get_mut(a2).next = an;
    arena.get_mut(an).prev = a2;

    arena.get_mut(b2).next = a2;
    arena.get_mut(a2).prev = b2;

    arena.get_mut(bp).next = b2;
    arena.get_mut(b2).prev = bp;

    b2
}

/// signed (shoelace) area of a flat polygon ring
fn signed_area(data: &[f64], start: usize, end: usize, dim: usize) -> f64 {
    let mut sum = 0.0;
    let mut j = end.saturating_sub(dim);
    let mut i = start;
    while i < end {
        sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
        j = i;
        i += dim;
    }
    sum
}

/// Return the relative difference between the polygon area and the area of its
/// triangulation — a value near 0 means a correct triangulation. Useful for verifying
/// output, e.g. in tests.
pub fn deviation<N: EarcutIndex>(
    data: &[f64],
    hole_indices: &[usize],
    dim: usize,
    triangles: &[N],
) -> f64 {
    let has_holes = !hole_indices.is_empty();
    let outer_len = if has_holes {
        hole_indices[0] * dim
    } else {
        data.len()
    };

    let mut polygon_area = signed_area(data, 0, outer_len, dim).abs();
    if has_holes {
        let len = hole_indices.len();
        for i in 0..len {
            let start = hole_indices[i] * dim;
            let end = if i < len - 1 {
                hole_indices[i + 1] * dim
            } else {
                data.len()
            };
            polygon_area -= signed_area(data, start, end, dim).abs();
        }
    }

    let mut triangles_area = 0.0;
    let mut i = 0;
    while i < triangles.len() {
        let a = triangles[i].to_usize() * dim;
        let b = triangles[i + 1].to_usize() * dim;
        let c = triangles[i + 2].to_usize() * dim;
        triangles_area += ((data[a] - data[c]) * (data[b + 1] - data[a + 1])
            - (data[a] - data[b]) * (data[c + 1] - data[a + 1]))
            .abs();
        i += 3;
    }

    if polygon_area == 0.0 && triangles_area == 0.0 {
        0.0
    } else {
        ((triangles_area - polygon_area) / polygon_area).abs()
    }
}

/// Turn a polygon in ring form (outer ring followed by hole rings, each a slice of `[x, y]`
/// points) into the flat form `earcut` accepts: `(vertices, hole_indices, dimensions)`.
pub fn flatten(rings: &[Vec<[f64; 2]>]) -> (Vec<f64>, Vec<usize>, usize) {
    let mut vertices = Vec::new();
    let mut holes = Vec::new();
    let dimensions = 2;
    let mut hole_index = 0;
    let mut prev_len = 0;

    for ring in rings {
        for p in ring {
            vertices.push(p[0]);
            vertices.push(p[1]);
        }
        if prev_len > 0 {
            hole_index += prev_len;
            holes.push(hole_index);
        }
        prev_len = ring.len();
    }

    (vertices, holes, dimensions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_quad() {
        let data = [10.0, 0.0, 0.0, 50.0, 60.0, 60.0, 70.0, 10.0];
        let triangles: Vec<u32> = earcut(&data, &[], 2);
        assert_eq!(triangles.len(), 6);
        assert!(deviation(&data, &[], 2, &triangles) < 1e-9);
    }

    #[test]
    fn square_with_hole() {
        let data = [
            0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0, 20.0, 20.0, 80.0, 20.0, 80.0, 80.0,
            20.0, 80.0,
        ];
        let triangles: Vec<u32> = earcut(&data, &[4], 2);
        assert_eq!(triangles.len(), 8 * 3);
        assert!(deviation(&data, &[4], 2, &triangles) < 1e-9);
    }

    #[test]
    fn empty_input() {
        let data: [f64; 0] = [];
        let triangles: Vec<u32> = earcut(&data, &[], 2);
        assert!(triangles.is_empty());
    }

    #[test]
    fn degenerate_triangle() {
        let data = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let triangles: Vec<u32> = earcut(&data, &[], 2);
        assert!(triangles.is_empty());
    }

    #[test]
    fn reusable_earcut_matches_one_shot() {
        let quad = [10.0, 0.0, 0.0, 50.0, 60.0, 60.0, 70.0, 10.0];
        let hole_square = [
            0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0, 20.0, 20.0, 80.0, 20.0, 80.0, 80.0,
            20.0, 80.0,
        ];

        let mut earcutter = Earcut::new();
        let mut out: Vec<u32> = Vec::new();

        // triangulate a few different shapes back-to-back on the same `Earcut`, and check
        // each result matches what the one-shot `earcut` function returns.
        for _ in 0..3 {
            earcutter.earcut_into(&quad, &[], 2, &mut out);
            assert_eq!(out, earcut::<u32>(&quad, &[], 2));

            earcutter.earcut_into(&hole_square, &[4], 2, &mut out);
            assert_eq!(out, earcut::<u32>(&hole_square, &[4], 2));

            earcutter.earcut_into(&[], &[], 2, &mut out);
            assert!(out.is_empty());
        }
    }
}
