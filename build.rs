fn main() {
    // Only compile the vendored earcut.hpp FFI wrapper when explicitly requested, so building
    // and testing `rearcut` itself never requires a C++ toolchain.
    if std::env::var_os("CARGO_FEATURE_EARCUT_HPP").is_some() {
        println!("cargo:rerun-if-changed=cpp/earcut_ffi.cpp");
        println!("cargo:rerun-if-changed=cpp/earcut.hpp");
        cc::Build::new()
            .cpp(true)
            .std("c++14")
            .file("cpp/earcut_ffi.cpp")
            .include("cpp")
            .opt_level(3)
            .compile("earcut_hpp_ffi");
    }
}
