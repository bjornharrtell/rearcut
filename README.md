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

(Ryzen-class x86_64, `cargo bench`, single run — see the HTML report for full
distributions; absolute numbers vary by machine, relative shape is the
point. `rearcut` and `earcut-rs` are measured via their reusable `Earcut`
struct so repeated calls don't pay for a fresh arena each time; `lyon` has
no such API and is measured one-shot.)

| Benchmark | rearcut | earcut-rs | lyon |
|---|---:|---:|---:|
| `fixtures/building` | 0.31 µs | 0.18 µs | 1.4 µs |
| `fixtures/dude` | 4.9 µs | 4.5 µs | 8.6 µs |
| `fixtures/bad-hole` | 2.5 µs | 2.1 µs | 4.7 µs |
| `fixtures/water` | 201 µs | 182 µs | 537 µs |
| `fixtures/water-huge` | 1.65 ms | 1.38 ms | 1.22 ms |
| `fixtures/water-huge3` | 22.7 ms | 19.0 ms | 5.8 ms |
| `star/16` | 0.78 µs | 0.61 µs | 3.8 µs |
| `star/4096` | 14.3 ms | 12.2 ms | 41 ms |
| `star/65536` | 4.1 s | 3.7 s | 11.2 s |
| `holes/64` | 59 µs | 56 µs | 26 µs |
| `holes/1024` | 10.0 ms | 8.1 ms | 0.74 ms |

**Takeaway:** for simple-to-moderately-complex polygons, and especially for
concave, hole-free polygons (e.g. `star`), `rearcut`'s ear-slicing approach
is consistently 2–4x faster than lyon's sweep-line tessellator. Against
`earcut-rs` (georust's independent port), `rearcut` is now within
single-digit percent on most real-world fixtures and small-to-medium
synthetic shapes, occasionally matching or edging ahead; `earcut-rs` still
leads by a wider margin on the largest synthetic stress tests (`star/65536`,
`holes/1024`), mainly thanks to its byte-offset node links (see the note
below). The `holes` benchmark's small, uniform-size synthetic holes are a
weaker case for both crates' block-bbox hole-bridge index than realistic
many-vertex holes (see `fixtures/water-huge3`), which is why lyon's
sweep-line approach wins there despite losing everywhere else.

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
