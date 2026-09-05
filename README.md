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

| Benchmark | rearcut | earcut-rs | earcut.hpp | lyon |
|---|---:|---:|---:|---:|
| `fixtures/building` | 137 ns | 151 ns | 241 ns | 1.38 µs |
| `fixtures/dude` | 4.05 µs | 4.17 µs | 3.88 µs | 8.18 µs |
| `fixtures/bad-hole` | 1.54 µs | 1.92 µs | 2.08 µs | 4.21 µs |
| `fixtures/water3b` | 702 ns | 886 ns | 1.07 µs | 2.56 µs |
| `fixtures/water4` | 45.2 µs | 50.0 µs | 46.6 µs | 66.3 µs |
| `fixtures/water2` | 150 µs | 137 µs | 157 µs | 117 µs |
| `fixtures/water` | 168 µs | 167 µs | 208 µs | 522 µs |
| `fixtures/water-huge` | 1.32 ms | 1.27 ms | 1.39 ms | 1.18 ms |
| `fixtures/water-huge3` | 18.3 ms | 18.0 ms | 21.3 ms | 4.78 ms |
| `star/16` | 523 ns | 528 ns | 576 ns | 3.42 µs |
| `star/256` | 46.6 µs | 51.6 µs | 46.6 µs | 182 µs |
| `star/4096` | 11.8 ms | 12.1 ms | 14.3 ms | 33.3 ms |
| `holes/4` | 543 ns | 724 ns | 910 ns | 2.02 µs |
| `holes/64` | 22.5 µs | 44.3 µs | 42.2 µs | 24.8 µs |
| `holes/256` | 179 µs | 617 µs | 616 µs | 117 µs |
| `holes/1024` | 2.08 ms | 8.12 ms | 7.48 ms | 605 µs |

**Takeaway:** `rearcut` is at or near parity with `earcut-rs` and
`earcut.hpp` on typical inputs, and pulls ahead as size or hole count
grows: ~4x faster than both on the many-hole `holes/1024` grid. The main structural change
behind that is the hole-bridge index: each block gets an explicit node
list rather than a ring range, so a block's scan cost stays bounded by
its own size even after `split_polygon` splices a hole into the middle
of it — this is a strict win with no trade-off, since it never scans
more per block than it did before.

Against lyon's sweep-line tessellator, ear clipping wins by 2–5x on
concave, hole-free polygons (`star`) and on fixtures with a few large
holes, but loses on inputs with very many small holes (`holes/1024`,
`water-huge3`) — a sweep line is simply the better asymptotic fit there.

Note that `.cargo/config.toml` builds with `-C target-cpu=native` (and
`build.rs` passes `-march=native` to `earcut.hpp`), so these numbers are
for host-tuned builds of every implementation. It is worth a few percent
on the `star` cases. It applies to builds run from inside this
repository only, and is not inherited by crates that depend on a
published `rearcut`.

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
