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
distributions; absolute numbers vary by machine, relative shape is the point)

| Benchmark | rearcut | earcut-rs | lyon | earcut.hpp |
|---|---:|---:|---:|---:|
| `fixtures/building` | 0.33 µs | 0.20 µs | 1.6 µs | 0.25 µs |
| `fixtures/dude` | 5.7 µs | 4.7 µs | 8.6 µs | 4.1 µs |
| `fixtures/bad-hole` | 3.0 µs | 2.3 µs | 4.5 µs | 2.3 µs |
| `fixtures/water` | 248 µs | 200 µs | 542 µs | 203 µs |
| `fixtures/water-huge` | 1.80 ms | 1.39 ms | 1.23 ms | 1.46 ms |
| `fixtures/water-huge3` | 26.1 ms | 19.4 ms | 5.2 ms | 22.0 ms |
| `star/16` | 0.99 µs | 0.70 µs | 3.8 µs | 0.6 µs |
| `star/4096` | 16.6 ms | 13.0 ms | 41 ms | 14.8 ms |
| `star/65536` | 4.6 s | 3.4 s | 9.9 s | 4.3 s |
| `holes/64` | 66 µs | 50 µs | 27 µs | 44 µs |
| `holes/1024` | 11.3 ms | 8.0 ms | 0.68 ms | 7.5 ms |

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

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
