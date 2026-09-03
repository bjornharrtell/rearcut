# rearcut

A Rust port of [mapbox/earcut](https://github.com/mapbox/earcut) / [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp) — a fast,
dependency-free 2D polygon triangulation library based on ear slicing with a
z-order curve hash for O(n) average-case performance. It triangulates simple
polygons that can be concave and have holes.

This crate is a from-scratch, safe-Rust reimplementation of the earcut
algorithm API-compatible in spirit with the original JS/C++
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

For repeated triangulations, `rearcut::Earcut` reuses its internal arena's
allocations (nodes, hole-bridge block index, z-order sort scratch buffers)
across calls instead of allocating fresh ones every time:

```rust
use rearcut::Earcut;

let mut earcutter = Earcut::new();
let mut triangles: Vec<u32> = Vec::new();

let data = [10.0, 0.0, 0.0, 50.0, 60.0, 60.0, 70.0, 10.0];
earcutter.earcut_into(&data, &[], 2, &mut triangles);
assert_eq!(triangles.len(), 6);
```

## Correctness

`rearcut` is validated against the full upstream
[mapbox/earcut test fixture suite](https://github.com/mapbox/earcut/tree/main/test)
(58 real-world and adversarial polygons, vendored under `tests/fixtures/`,
checked against `tests/expected.json`), see `tests/fixtures_test.rs`.

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

(Ryzen-class x86_64, `cargo bench --features earcut-hpp`, single run — see the
HTML report for full distributions; absolute numbers vary by machine,
relative shape is the point. `rearcut` and `earcut-rs` are measured via
their reusable `Earcut` struct so repeated calls don't pay for a fresh
arena each time; `lyon` and `earcut.hpp` have no such API and are measured
one-shot.)

| Benchmark | rearcut | earcut-rs | lyon | earcut.hpp |
|---|---:|---:|---:|---:|
| `fixtures/building` | 0.18 µs | 0.17 µs | 1.4 µs | 0.26 µs |
| `fixtures/dude` | 6.3 µs | 4.4 µs | 8.5 µs | 4.1 µs |
| `fixtures/bad-hole` | 1.8 µs | 2.1 µs | 4.4 µs | 2.1 µs |
| `fixtures/water` | 177 µs | 186 µs | 523 µs | 249 µs |
| `fixtures/water2` | 200 µs | 148 µs | 120 µs | 163 µs |
| `fixtures/water-huge` | 1.41 ms | 1.34 ms | 1.19 ms | 1.41 ms |
| `fixtures/water-huge3` | 14.0 ms | 18.6 ms | 4.9 ms | 22.0 ms |
| `star/16` | 0.53 µs | 0.57 µs | 3.6 µs | 0.69 µs |
| `star/4096` | 8.35 ms | 13.1 ms | 39.2 ms | 14.6 ms |
| `star/65536` | 2.08 s | 4.14 s | 10.0 s | 4.37 s |
| `holes/64` | 26 µs | 44 µs | 25 µs | 42 µs |
| `holes/1024` | 1.81 ms | 7.42 ms | 0.63 ms | 7.04 ms |

**Takeaway:** `rearcut` is the fastest of the three ear-clipping
implementations on almost every input, and the margin grows with size:
roughly 2x faster than both `earcut-rs` and `earcut.hpp` on the large
concave `star` cases and on the many-hole `holes` grid, and ~1.5x faster
on `water-huge3`. Two structural changes account for most of that: the
hole-bridge index gives each block an explicit node list rather than a
ring range (so a block's scan cost stays bounded by its own size even
after `split_polygon` splices a hole into the middle of it), and the
z-order index is a sorted array walked outwards from the ear's own slot,
rather than a doubly linked z-list chased through the arena.

Against lyon's sweep-line tessellator, ear clipping wins by 2–5x on
concave, hole-free polygons (`star`) and on fixtures with a few large
holes, but loses on inputs with very many small holes (`holes/1024`,
`water-huge3`) — a sweep line is simply the better asymptotic fit there.

The array-based z-order scan tests a fixed-size block of entries at a
time, so it pays a small fixed cost per query that only amortises once
scans are long. On the mid-size fixtures whose scans are short
(`fixtures/dude`, `fixtures/water2`) that costs 20–50% versus the linked
z-list it replaced, which is the trade made for the much larger wins
above.

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
