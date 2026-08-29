//! Benchmarks comparing `rearcut` against `lyon_tessellation`'s fill tessellator, and
//! (with `--features earcut-hpp`) the reference `mapbox/earcut.hpp` C++ implementation it was
//! ported from, on the same set of real-world and synthetic polygons (a subset of the
//! upstream mapbox/earcut fixture suite, plus a few synthetic shapes of varying
//! size/complexity).
//!
//! Run with: `cargo bench` (rearcut vs lyon)
//! or:       `cargo bench --features earcut-hpp` (also includes earcut.hpp; requires a C++14
//!           compiler, invoked via `cc` in build.rs)
//! HTML report: `target/criterion/report/index.html`

use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use earcut::Earcut as GeorustEarcut;
use lyon::math::point;
use lyon::path::Path;
use lyon::tessellation::{BuffersBuilder, FillOptions, FillTessellator, FillVertex, VertexBuffers};
use serde::Deserialize;
use std::fs;
use std::path::Path as FsPath;

#[cfg(feature = "earcut-hpp")]
mod earcut_hpp_ffi {
    #[repr(C)]
    pub struct EarcutHppResult {
        pub data: *mut u32,
        pub len: usize,
    }

    unsafe extern "C" {
        pub fn earcut_hpp_triangulate(
            flat_data: *const f64,
            flat_len: usize,
            hole_indices: *const usize,
            hole_count: usize,
            dim: usize,
        ) -> EarcutHppResult;
        pub fn earcut_hpp_free(result: EarcutHppResult);
    }

    /// Safe wrapper: triangulates via the vendored `earcut.hpp` and returns the index count
    /// (the indices themselves are freed immediately, matching what the other benchmarked
    /// implementations do by only measuring/returning the triangle count).
    pub fn triangulate(vertices: &[f64], holes: &[usize], dim: usize) -> usize {
        unsafe {
            let result = earcut_hpp_triangulate(
                vertices.as_ptr(),
                vertices.len(),
                holes.as_ptr(),
                holes.len(),
                dim,
            );
            let len = result.len;
            earcut_hpp_free(result);
            len
        }
    }
}

#[derive(Deserialize)]
struct Fixture(Vec<Vec<[f64; 2]>>);

fn load_fixture(name: &str) -> Vec<Vec<[f64; 2]>> {
    let path = FsPath::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(format!("{name}.json"));
    let raw = fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {path:?}: {e}"));
    serde_json::from_str::<Fixture>(&raw)
        .unwrap_or_else(|e| panic!("parsing {path:?}: {e}"))
        .0
}

/// Generate a synthetic star-shaped polygon (a common worst-case-ish shape for ear clipping:
/// concave, non-convex, no holes) with `points` outer vertices.
fn star(points: usize, outer_r: f64, inner_r: f64) -> Vec<Vec<[f64; 2]>> {
    let n = points * 2;
    let mut ring = Vec::with_capacity(n);
    for i in 0..n {
        let r = if i % 2 == 0 { outer_r } else { inner_r };
        let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
        ring.push([r * theta.cos(), r * theta.sin()]);
    }
    vec![ring]
}

/// Generate a grid of `holes_per_side^2` small square holes punched out of a large square,
/// exercising the hole-bridging path at scale.
fn grid_with_holes(holes_per_side: usize) -> Vec<Vec<[f64; 2]>> {
    let outer_size = 1000.0;
    let margin = 20.0;
    let cell = (outer_size - 2.0 * margin) / holes_per_side as f64;
    let hole_size = cell * 0.6;

    let mut rings = vec![vec![
        [0.0, 0.0],
        [outer_size, 0.0],
        [outer_size, outer_size],
        [0.0, outer_size],
    ]];

    for row in 0..holes_per_side {
        for col in 0..holes_per_side {
            let x0 = margin + col as f64 * cell + (cell - hole_size) / 2.0;
            let y0 = margin + row as f64 * cell + (cell - hole_size) / 2.0;
            // holes must be wound opposite to the outer ring; both `rearcut::flatten` and
            // lyon's path builder handle winding normalization for us regardless
            rings.push(vec![
                [x0, y0],
                [x0, y0 + hole_size],
                [x0 + hole_size, y0 + hole_size],
                [x0 + hole_size, y0],
            ]);
        }
    }

    rings
}

/// Triangulates with a caller-provided, reused `rearcut::Earcut` and output buffer, so
/// repeated calls (as in a `b.iter` loop) measure steady-state performance without paying
/// for a fresh arena allocation every time. `vertices`/`holes` are pre-flattened once
/// outside the timing loop (mirroring how `earcut-rs` reads directly from existing ring
/// data with no separate flatten step), so only the triangulation itself is measured.
fn rearcut_triangulate(
    earcutter: &mut rearcut::Earcut,
    triangles: &mut Vec<u32>,
    vertices: &[f64],
    holes: &[usize],
    dim: usize,
) -> usize {
    earcutter.earcut_into(vertices, holes, dim, triangles);
    triangles.len()
}

/// Same as `rearcut_triangulate`, but for the reused `earcut-rs` (georust) `Earcut`.
fn georust_earcut_triangulate(
    earcutter: &mut GeorustEarcut<f64>,
    triangles: &mut Vec<u32>,
    rings: &[Vec<[f64; 2]>],
) -> usize {
    let mut hole_indices: Vec<u32> = Vec::with_capacity(rings.len().saturating_sub(1));
    let mut count = 0u32;
    for (i, ring) in rings.iter().enumerate() {
        if i > 0 {
            hole_indices.push(count);
        }
        count += ring.len() as u32;
    }
    let data = rings.iter().flatten().copied();
    earcutter.earcut(data, &hole_indices, triangles);
    triangles.len()
}

fn lyon_triangulate(rings: &[Vec<[f64; 2]>]) -> usize {
    let mut builder = Path::builder();
    for ring in rings {
        let mut points_iter = ring.iter();
        let Some(first) = points_iter.next() else {
            continue;
        };
        builder.begin(point(first[0] as f32, first[1] as f32));
        for p in points_iter {
            builder.line_to(point(p[0] as f32, p[1] as f32));
        }
        builder.end(true);
    }
    let path = builder.build();

    let mut buffers: VertexBuffers<[f32; 2], u32> = VertexBuffers::new();
    let mut tessellator = FillTessellator::new();
    let result = tessellator.tessellate_path(
        &path,
        &FillOptions::default().with_tolerance(0.05),
        &mut BuffersBuilder::new(&mut buffers, |v: FillVertex| {
            let p = v.position();
            [p.x, p.y]
        }),
    );
    result.expect("lyon tessellation failed");
    buffers.indices.len()
}

fn bench_fixtures(c: &mut Criterion) {
    let fixtures = [
        "building",
        "dude",
        "bad-hole",
        "water3b",
        "water4",
        "water2",
        "water",
        "water-huge",
        "water-huge3",
    ];

    let mut group = c.benchmark_group("fixtures");
    for name in fixtures {
        let rings = load_fixture(name);
        let vertex_count: usize = rings.iter().map(|r| r.len()).sum();
        group.throughput(Throughput::Elements(vertex_count as u64));

        let (rearcut_vertices, rearcut_holes, rearcut_dim) = rearcut::flatten(&rings);
        group.bench_with_input(
            BenchmarkId::new("rearcut", name),
            &(rearcut_vertices, rearcut_holes, rearcut_dim),
            |b, (vertices, holes, dim)| {
                let mut earcutter = rearcut::Earcut::new();
                let mut triangles = Vec::new();
                b.iter(|| {
                    black_box(rearcut_triangulate(
                        &mut earcutter,
                        &mut triangles,
                        black_box(vertices),
                        black_box(holes),
                        *dim,
                    ))
                });
            },
        );
        group.bench_with_input(BenchmarkId::new("earcut-rs", name), &rings, |b, rings| {
            let mut earcutter = GeorustEarcut::new();
            let mut triangles = Vec::new();
            b.iter(|| {
                black_box(georust_earcut_triangulate(
                    &mut earcutter,
                    &mut triangles,
                    black_box(rings),
                ))
            });
        });
        group.bench_with_input(BenchmarkId::new("lyon", name), &rings, |b, rings| {
            b.iter(|| black_box(lyon_triangulate(black_box(rings))));
        });
        #[cfg(feature = "earcut-hpp")]
        {
            let (vertices, holes, dim) = rearcut::flatten(&rings);
            group.bench_with_input(
                BenchmarkId::new("earcut.hpp", name),
                &(vertices, holes, dim),
                |b, (vertices, holes, dim)| {
                    b.iter(|| {
                        black_box(earcut_hpp_ffi::triangulate(
                            black_box(vertices),
                            black_box(holes),
                            *dim,
                        ))
                    });
                },
            );
        }
    }
    group.finish();
}

fn bench_stars(c: &mut Criterion) {
    let mut group = c.benchmark_group("star");
    for &points in &[16usize, 256, 4096, 65536] {
        let rings = star(points, 100.0, 40.0);
        group.throughput(Throughput::Elements(points as u64 * 2));

        let (rearcut_vertices, rearcut_holes, rearcut_dim) = rearcut::flatten(&rings);
        group.bench_with_input(
            BenchmarkId::new("rearcut", points),
            &(rearcut_vertices, rearcut_holes, rearcut_dim),
            |b, (vertices, holes, dim)| {
                let mut earcutter = rearcut::Earcut::new();
                let mut triangles = Vec::new();
                b.iter(|| {
                    black_box(rearcut_triangulate(
                        &mut earcutter,
                        &mut triangles,
                        black_box(vertices),
                        black_box(holes),
                        *dim,
                    ))
                });
            },
        );
        group.bench_with_input(BenchmarkId::new("earcut-rs", points), &rings, |b, rings| {
            let mut earcutter = GeorustEarcut::new();
            let mut triangles = Vec::new();
            b.iter(|| {
                black_box(georust_earcut_triangulate(
                    &mut earcutter,
                    &mut triangles,
                    black_box(rings),
                ))
            });
        });
        group.bench_with_input(BenchmarkId::new("lyon", points), &rings, |b, rings| {
            b.iter(|| black_box(lyon_triangulate(black_box(rings))));
        });
        #[cfg(feature = "earcut-hpp")]
        {
            let (vertices, holes, dim) = rearcut::flatten(&rings);
            group.bench_with_input(
                BenchmarkId::new("earcut.hpp", points),
                &(vertices, holes, dim),
                |b, (vertices, holes, dim)| {
                    b.iter(|| {
                        black_box(earcut_hpp_ffi::triangulate(
                            black_box(vertices),
                            black_box(holes),
                            *dim,
                        ))
                    });
                },
            );
        }
    }
    group.finish();
}

fn bench_holes(c: &mut Criterion) {
    let mut group = c.benchmark_group("holes");
    for &holes_per_side in &[2usize, 8, 16, 32] {
        let rings = grid_with_holes(holes_per_side);
        let hole_count = holes_per_side * holes_per_side;
        group.throughput(Throughput::Elements(hole_count as u64));

        let (rearcut_vertices, rearcut_holes, rearcut_dim) = rearcut::flatten(&rings);
        group.bench_with_input(
            BenchmarkId::new("rearcut", hole_count),
            &(rearcut_vertices, rearcut_holes, rearcut_dim),
            |b, (vertices, holes, dim)| {
                let mut earcutter = rearcut::Earcut::new();
                let mut triangles = Vec::new();
                b.iter(|| {
                    black_box(rearcut_triangulate(
                        &mut earcutter,
                        &mut triangles,
                        black_box(vertices),
                        black_box(holes),
                        *dim,
                    ))
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("earcut-rs", hole_count),
            &rings,
            |b, rings| {
                let mut earcutter = GeorustEarcut::new();
                let mut triangles = Vec::new();
                b.iter(|| {
                    black_box(georust_earcut_triangulate(
                        &mut earcutter,
                        &mut triangles,
                        black_box(rings),
                    ))
                });
            },
        );
        group.bench_with_input(BenchmarkId::new("lyon", hole_count), &rings, |b, rings| {
            b.iter(|| black_box(lyon_triangulate(black_box(rings))));
        });
        #[cfg(feature = "earcut-hpp")]
        {
            let (vertices, holes, dim) = rearcut::flatten(&rings);
            group.bench_with_input(
                BenchmarkId::new("earcut.hpp", hole_count),
                &(vertices, holes, dim),
                |b, (vertices, holes, dim)| {
                    b.iter(|| {
                        black_box(earcut_hpp_ffi::triangulate(
                            black_box(vertices),
                            black_box(holes),
                            *dim,
                        ))
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_fixtures, bench_stars, bench_holes);
criterion_main!(benches);
