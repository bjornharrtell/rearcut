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
| `fixtures/building` | 141 ns | 151 ns | 250 ns | 1.37 µs |
| `fixtures/dude` | 4.07 µs | 4.00 µs | 3.99 µs | 8.39 µs |
| `fixtures/bad-hole` | 1.77 µs | 1.90 µs | 2.05 µs | 4.27 µs |
| `fixtures/water3b` | 858 ns | 869 ns | 1.02 µs | 2.52 µs |
| `fixtures/water4` | 45.5 µs | 47.4 µs | 47.1 µs | 64.2 µs |
| `fixtures/water2` | 155 µs | 135 µs | 158 µs | 115 µs |
| `fixtures/water` | 167 µs | 171 µs | 195 µs | 513 µs |
| `fixtures/water-huge` | 1.34 ms | 1.33 ms | 1.45 ms | 1.21 ms |
| `fixtures/water-huge3` | 19.1 ms | 18.6 ms | 21.0 ms | 4.85 ms |
| `star/16` | 555 ns | 536 ns | 572 ns | 3.45 µs |
| `star/256` | 47.9 µs | 51.6 µs | 47.1 µs | 180 µs |
| `star/4096` | 11.8 ms | 12.2 ms | 14.4 ms | 32.8 ms |
| `star/65536` | 3.55 s | 3.73 s | 4.68 s | 8.16 s |
| `holes/4` | 669 ns | 731 ns | 905 ns | 2.05 µs |
| `holes/64` | 22.9 µs | 45.3 µs | 43.0 µs | 25.5 µs |
| `holes/256` | 189 µs | 618 µs | 632 µs | 115 µs |
| `holes/1024` | 2.18 ms | 8.31 ms | 7.83 ms | 609 µs |

**Takeaway:** `rearcut` is at or near parity with `earcut-rs` and
`earcut.hpp` on typical inputs, and pulls ahead as size or hole count
grows: ~5% faster than `earcut-rs` on `star/65536`, and ~4x faster than
both on the many-hole `holes/1024` grid. The main structural change
behind that is the hole-bridge index: each block gets an explicit node
list rather than a ring range, so a block's scan cost stays bounded by
its own size even after `split_polygon` splices a hole into the middle
of it — this is a strict win with no trade-off, since it never scans
more per block than it did before.

A z-order scan implemented as a sorted array walked outward from the
ear's own slot (instead of a doubly linked z-list chased through the
arena) was also tried, and measured a further ~1.7x on `star/65536` and
~4.5x on `holes/1024`. It was **not** kept: it pays a fixed per-query
cost that only amortises on long scans, which made it 20-60% slower on
mid-size fixtures with short scans (`fixtures/dude`, `fixtures/water2`)
that are more representative of typical inputs than the synthetic
`star`/`holes` stress tests. Reverted in favor of keeping those
realistic cases fast.

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
