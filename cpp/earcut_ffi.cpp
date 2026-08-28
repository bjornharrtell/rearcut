// Thin C ABI wrapper around mapbox/earcut.hpp, used only by the criterion benchmarks
// (benches/triangulate.rs) to compare rearcut's performance against the reference C++
// implementation it was ported from.
#include "earcut.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

extern "C" {

struct EarcutHppResult {
    uint32_t* data;
    size_t len;
};

// `flat_data` is a flat array of `dim`-tuples (x, y first); `hole_indices` lists the vertex
// index (not coordinate index) at which each hole ring starts, mirroring rearcut::earcut's
// and earcut.js's flat input convention.
EarcutHppResult earcut_hpp_triangulate(const double* flat_data, size_t flat_len,
                                        const size_t* hole_indices, size_t hole_count,
                                        size_t dim) {
    using Point = std::array<double, 2>;
    using Ring = std::vector<Point>;

    size_t num_points = dim == 0 ? 0 : flat_len / dim;

    std::vector<size_t> boundaries;
    boundaries.reserve(hole_count + 1);
    for (size_t i = 0; i < hole_count; i++) boundaries.push_back(hole_indices[i]);
    boundaries.push_back(num_points);

    std::vector<Ring> polygon;
    polygon.reserve(boundaries.size());

    size_t prev = 0;
    for (size_t b : boundaries) {
        Ring ring;
        ring.reserve(b > prev ? b - prev : 0);
        for (size_t v = prev; v < b; v++) {
            ring.push_back({flat_data[v * dim], flat_data[v * dim + 1]});
        }
        polygon.push_back(std::move(ring));
        prev = b;
    }

    std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);

    size_t n = indices.size();
    uint32_t* out = n ? static_cast<uint32_t*>(std::malloc(n * sizeof(uint32_t))) : nullptr;
    if (out) std::memcpy(out, indices.data(), n * sizeof(uint32_t));

    return EarcutHppResult{out, n};
}

void earcut_hpp_free(EarcutHppResult r) {
    std::free(r.data);
}

} // extern "C"
