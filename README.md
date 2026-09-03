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
| `fixtures/building` | 145 ns | 154 ns | 247 ns | 1.40 µs |
| `fixtures/dude` | 6.18 µs | 3.99 µs | 4.02 µs | 8.72 µs |
| `fixtures/bad-hole` | 1.77 µs | 1.98 µs | 2.08 µs | 4.34 µs |
| `fixtures/water3b` | 853 ns | 866 ns | 1.01 µs | 2.60 µs |
| `fixtures/water4` | 57.9 µs | 49.1 µs | 47.0 µs | 66.1 µs |
| `fixtures/water2` | 208 µs | 135 µs | 157 µs | 118 µs |
| `fixtures/water` | 178 µs | 168 µs | 203 µs | 511 µs |
| `fixtures/water-huge` | 1.43 ms | 1.25 ms | 1.38 ms | 1.17 ms |
| `fixtures/water-huge3` | 13.8 ms | 18.0 ms | 21.5 ms | 4.60 ms |
| `star/16` | 518 ns | 540 ns | 569 ns | 3.43 µs |
| `star/256` | 47.9 µs | 51.4 µs | 45.2 µs | 176 µs |
| `star/4096` | 7.67 ms | 11.65 ms | 13.60 ms | 29.33 ms |
| `star/65536` | 1.88 s | 3.19 s | 4.03 s | 7.30 s |
| `holes/4` | 644 ns | 760 ns | 884 ns | 1.99 µs |
| `holes/64` | 26.1 µs | 43.3 µs | 40.8 µs | 24.7 µs |
| `holes/256` | 181 µs | 567 µs | 559 µs | 112 µs |
| `holes/1024` | 1.83 ms | 7.70 ms | 7.62 ms | 593 µs |

**Takeaway:** `rearcut` is the fastest of the three ear-clipping
implementations on most inputs, and its margin grows with size: ~1.7x
faster than `earcut-rs` and ~2.1x faster than `earcut.hpp` on
`star/65536`, and ~4x faster than both on the many-hole `holes/1024`
grid. Two structural changes account for most of that. The hole-bridge
index gives each block an explicit node list rather than a ring range, so
a block's scan cost stays bounded by its own size even after
`split_polygon` splices a hole into the middle of it. And the z-order
index is a sorted array walked outwards from the ear's own slot, rather
than a doubly linked z-list chased through the arena — it rejects a whole
aligned block of candidates with one branchless filter, and its entries
are sequential in memory instead of a pointer chase.

That array scan pays a small fixed cost per query which only amortises
once scans are long, so on the mid-size fixtures whose scans are short
(`fixtures/dude`, `fixtures/water2`) it is 20–50% slower than the linked
z-list it replaced. That is the trade made for the much larger wins
above.

Against lyon's sweep-line tessellator, ear clipping wins by 2–5x on
concave, hole-free polygons (`star`) and on fixtures with a few large
holes, but loses on inputs with very many small holes (`holes/1024`,
`water-huge3`) — a sweep line is simply the better asymptotic fit there.

Note that `.cargo/config.toml` builds with `-C target-cpu=native` (and
`build.rs` passes `-march=native` to `earcut.hpp`), so these numbers are
for host-tuned builds of every implementation. It is worth 3–5% on the
`star` cases. It applies to builds run from inside this repository only,
and is not inherited by crates that depend on a published `rearcut`.

## Acknowledgements

Ported from [mapbox/earcut](https://github.com/mapbox/earcut) (ISC License)
and cross-checked against [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp).
