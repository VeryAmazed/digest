use std::env;
use std::path::PathBuf;

fn main() {
    // 1. Try to get DIGEST_DIR from environment
    let digest_dir = env::var("DIGEST_DIR")
        .map(PathBuf::from)
        // 2. Fallback: try conda env or system install
        .or_else(|_| {
            // Try conda environment
            if let Ok(prefix) = env::var("CONDA_PREFIX") {
                let candidate = PathBuf::from(&prefix).join("include/digest");
                if candidate.exists() {
                    return Ok(PathBuf::from(prefix));
                }
            }
            // Try /usr/local or /usr
            for prefix in ["/usr/local", "/usr"] {
                let candidate = PathBuf::from(prefix).join("include/digest");
                if candidate.exists() {
                    return Ok(PathBuf::from(prefix));
                }
            }
            Err(env::VarError::NotPresent)
        })
        .unwrap_or_else(|_| {
            eprintln!(
                "Could not find Digest C++ library. \
                Please set the DIGEST_DIR environment variable to the install prefix \
                (containing 'include/digest' and 'lib/libnthash.a')."
            );
            std::process::exit(1);
        });

    // Build the C++ bindings
    cc::Build::new()
        .cpp(true)
        .file("src/bindings.cpp")
        .include(digest_dir.join("include"))
        .include(digest_dir.join("extern/nthash/include"))
        .flag("-std=c++17")
        .flag("-fPIC")
        .compile("bindings");

    // Link against nthash library
    println!("cargo:rustc-link-search=native={}/lib", digest_dir.display());
    println!("cargo:rustc-link-lib=static=nthash");

    // Link against our bindings library
    println!("cargo:rustc-link-lib=static=bindings");

    // Link against C++ standard library
    println!("cargo:rustc-link-lib=dylib=stdc++");
} 