# digest-rs

Rust bindings for the [digest](https://github.com/VeryAmazed/digest) C++ library, providing efficient kmer minimizer and syncmer digestion for DNA sequences.

## Requirements

- The C++ `digest` library must be available at `$DIGEST_DIR` (see `build.rs`) or installed via conda.

## Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
digest-rs = "0.1.0"
```

### Example

```rust
use digest_rs::{window_minimizer_rs, modimizer_rs, syncmer_rs};

let sequence = "ACGTACGT";
let k = 4;
let window = 2;
let mod_val = 3;

// Window minimizer
let minimizers = window_minimizer_rs(sequence, k, window)?;

// Modimizer
let modimizers = modimizer_rs(sequence, k, mod_val)?;

// Syncmer
let syncmers = syncmer_rs(sequence, k, window)?;
```

Each function returns a `Result<Vec<u32>, DigestError>` with the minimizer positions.

## Building

This crate uses a `build.rs` script to compile the C++ bindings and link against the C++ `digest` library.  
Make sure the C++ library and its dependencies are built and available at the expected paths.

## License

MIT