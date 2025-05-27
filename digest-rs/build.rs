use std::env;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let digest_dir = PathBuf::from(env::var("HOME").unwrap()).join("vast/digest");
    
    // Build the C++ bindings
    cc::Build::new()
        .cpp(true)
        .file("src/bindings.cpp")
        .include(digest_dir.join("include"))
        .include(digest_dir.join("extern/nthash/include"))
        .flag("-std=c++17")
        .flag("-fPIC")  // Add position-independent code flag
        .compile("bindings");  // This creates libbindings.a

    // Link against nthash library
    println!("cargo:rustc-link-search=native={}/lib", digest_dir.display());
    println!("cargo:rustc-link-lib=static=nthash");
    
    // Link against our bindings library
    println!("cargo:rustc-link-lib=static=bindings");
    
    // Link against C++ standard library
    println!("cargo:rustc-link-lib=dylib=stdc++");
} 