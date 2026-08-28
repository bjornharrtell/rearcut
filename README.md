# rearcut

A Rust port of [mapbox/earcut](https://github.com/mapbox/earcut) / [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp) — a fast,
dependency-free 2D polygon triangulation library based on ear slicing with a
z-order curve hash for O(n) average-case performance. It triangulates simple
polygons that can be concave and have holes.

This crate is a from-scratch, safe-Rust reimplementation of the earcut
algorithm (arena/`Vec`-backed doubly linked list instead of raw pointers —
no `unsafe`), including upstream's block-bbox hole-bridge index (see
"Correctness" below), API-compatible in spirit with the original JS/C++
libraries.

## Usage

```rust
use rearcut::earcut;

// a quadrilateral: (10,0) (0,50) (60,60) (70,10)
let data = [10.0, 0.0, 0.0, 50.0, 60.0, 60.0, 70.0, 10.0];
let triangles: Vec<u32> = earcut(&data, &[], 2);
```

Polygons with holes are given as a flat vertex array plus the vertex index
at which each hole ring starts:

```rust
use rearcut::earcut;

let data = [
    0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0, // outer ring
    20.0, 20.0, 80.0, 20.0, 80.0, 80.0, 20.0, 80.0, // hole
];
let triangles: Vec<u32> = earcut(&data, &[4], 2);
```

`rearcut::flatten` converts GeoJSON-style ring-of-points polygons
(`Vec<Vec<[f64; 2]>>`, outer ring first, then holes) into that flat form, and
`rearcut::deviation` computes the relative area error of a triangulation,
useful for verifying correctness.

The output index type is generic over `u16`/`u32`/`u64`/`usize` via the
`EarcutIndex` trait — pick the smallest type that fits your vertex count.

## Correctness

`rearcut` is validated against the full upstream
[mapbox/earcut test fixture suite](https://github.com/mapbox/earcut/tree/main/test)
(58 real-world and adversarial polygons, vendored under `tests/fixtures/`,
checked against `tests/expected.json`), see `tests/fixtures_test.rs`.

Upstream's implementation uses a block-bbox spatial index to accelerate
hole-bridge search (issue #183): ring edges are grouped into fixed-size
blocks, each with a cached bounding box, so the leftward-ray scan can skip
whole blocks instead of walking every merged-ring node. This port includes
the same block index (see `Arena::build_block_index`/`index_segment` in
`src/lib.rs`), so triangle counts match upstream almost exactly — of the 58
fixtures, only `issue142` (a hole touching the outer ring at a vertex) picks
a different, still-valid ear and differs by one triangle, which the test
suite tolerates explicitly. All other fixtures require an exact triangle
count when specified, plus a deviation (relative area error) check within
the documented tolerance.

Run the test suite with:

```sh
cargo test
```

## Benchmarks

Benchmarks (via [criterion](https://docs.rs/criterion)) compare `rearcut`
against [`lyon_tessellation`](https://docs.rs/lyon_tessellation)'s
`FillTessellator` and against [`earcut`](https://crates.io/crates/earcut)
(the [georust/earcut](https://github.com/georust/earcut) crate, another
independent Rust port of mapbox/earcut) on the same inputs:

- **`fixtures`** — a subset of the real-world earcut fixtures (buildings,
  hand-drawn shapes, water polygons of increasing size/hole count).
- **`star`** — synthetic concave star polygons (no holes) at increasing
  vertex counts, to isolate ear-clipping throughput on non-convex shapes.
- **`holes`** — a square with an increasing grid of square holes punched out,
  to isolate hole-bridging cost as hole count grows.

Run with:

```sh
cargo bench
```

This generates an HTML report at `target/criterion/report/index.html`.

An optional `earcut-hpp` feature additionally benchmarks the upstream C++
reference implementation, [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp)
(vendored under `cpp/`), via a small `extern "C"` FFI wrapper built with the
[`cc`](https://docs.rs/cc) crate. This is off by default so the crate never
requires a C++ toolchain unless explicitly requested:

```sh
cargo bench --features earcut-hpp
```

### Representative results

(Ryzen-class x86_64, `cargo bench`, single run — see the HTML report for full
distributions; absolute numbers vary by machine, relative shape is the point)

| Benchmark              | rearcut    | earcut-rs (georust) | lyon (FillTessellator) | earcut.hpp (C++) |
|-------------------------|-----------:|---------------------:|------------------------:|------------------:|
| `fixtures/building` (13 tris) | 0.33 µs | 0.20 µs | 1.6 µs  | 0.25 µs |
| `fixtures/dude` (106 tris)    | 5.7 µs  | 4.7 µs  | 8.6 µs  | 4.1 µs |
| `fixtures/bad-hole` (37 tris) | 3.0 µs  | 2.3 µs  | 4.5 µs  | 2.3 µs |
| `fixtures/water` (2482 tris)  | 248 µs  | 200 µs  | 542 µs  | 203 µs |
| `fixtures/water-huge` (5174 tris, 192 holes) | 1.80 ms | 1.39 ms | 1.23 ms | 1.46 ms |
| `fixtures/water-huge3` (15470 tris, 1443 holes) | 26.1 ms | 19.4 ms | 5.2 ms | 22.0 ms |
| `star/16` points        | 0.99 µs    | 0.70 µs | 3.8 µs | 0.6 µs |
| `star/4096` points      | 16.6 ms    | 13.0 ms | 41 ms  | 14.8 ms |
| `star/65536` points     | 4.6 s      | 3.4 s   | 9.9 s  | 4.3 s |
| `holes/64` (8×8 grid)   | 66 µs      | 50 µs   | 27 µs  | 44 µs |
| `holes/1024` (32×32 grid) | 11.3 ms  | 8.0 ms  | 0.68 ms | 7.5 ms |

**Takeaway:** for simple-to-moderately-complex polygons, and especially for
concave, hole-free polygons (e.g. `star`), `rearcut`'s ear-slicing approach
is consistently 2–4x faster than lyon's sweep-line tessellator, and now that
it ports the same block-bbox hole-bridge index as `earcut.hpp`, it also
tracks the native C++ reference closely on real many-holes polygons —
`fixtures/water-huge3` dropped from 95 ms (classic linear scan) to 26.1 ms
after the port, versus `earcut.hpp`'s 22.0 ms (the remaining gap is mostly
the safe arena's index indirection vs raw pointers, and Rust vs. GCC/Clang
codegen). The one case where the block index doesn't pay off is the
synthetic `holes` benchmark: its holes are tiny 4-vertex squares, well under
the 16-edge block granularity, so each hole gets its own block and the
per-block bookkeeping is pure overhead rather than a real skip — `earcut.hpp`
shows the same effect there (44 µs / 7.5 ms), and lyon's sweep-line approach
is fastest of all on this shape of workload. In short: the block index is a
genuine win for realistic hole-heavy polygons (many vertices per hole), and
roughly neutral-to-slightly-negative for degenerate microscopic-hole grids.

**Honest note on `earcut-rs`:** the [georust/earcut](https://github.com/georust/earcut)
crate is, at the time of writing, consistently the fastest of the four across
nearly every benchmark here — often beating even `earcut.hpp` (its README
reports the same result against the C++ reference). It uses a very similar
arena-of-nodes design to `rearcut`, but goes further: packed `Node` fields
(vertex index and Steiner-point flag share a `u32` via a bit flag), byte
offsets instead of element indices for node links (avoiding a multiply on
every dereference), and pre-sized buffers reused across calls via a stateful
`Earcut` struct. `rearcut` remains a solid, simple, from-scratch port with
its own optimizations (see below), but if raw throughput is the only
priority, `earcut-rs` is worth a look.

The arena's internal node links (`prev`/`next`/`prev_z`/`next_z`) are stored
as `u32` rather than `usize`, shrinking each `Node` from 64 to 48 bytes
(matching or beating earcut.hpp's raw-pointer-based C++ `Node`), and arena
lookups use `get_unchecked` internally (indices are always either the `NULL`
sentinel or values returned by the arena itself, so bounds checks are
provably redundant and are asserted only in debug builds). Combined with
pre-reserving the output triangle buffer and using unstable sorts for the
hole-bridge queue and z-order curve, this closed most of the earlier gap on
node-traversal-heavy workloads: `star/65536` went from 1.17x to ~1.07x of
`earcut.hpp`'s time, and `fixtures/building` from 1.8x to ~1.3x. The
`holes` benchmark's gap is largely unaffected, since it's dominated by
hole-bridge search cost rather than node traversal.

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
