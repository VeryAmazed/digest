use std::ffi::{c_char, CString};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DigestError {
    #[error("Invalid k-mer size: must be greater than 3")]
    InvalidKmerSize,
    #[error("Invalid window size: must be greater than 0")]
    InvalidWindowSize,
    #[error("Invalid mod value: must be greater than 0")]
    InvalidModValue,
    #[error("C++ library error: {0}")]
    LibraryError(String),
}

extern "C" {
    fn window_minimizer(seq: *const c_char, len: usize, k: u32, large_window: u32, out: *mut u32) -> usize;
    fn modimizer(seq: *const c_char, len: usize, k: u32, mod_val: u32, out: *mut u32) -> usize;
    fn syncmer(seq: *const c_char, len: usize, k: u32, large_window: u32, out: *mut u32) -> usize;
}

pub fn window_minimizer_rs(seq: &str, k: u32, window: u32) -> Result<Vec<u32>, DigestError> {
    if k < 4 {
        return Err(DigestError::InvalidKmerSize);
    }
    if window == 0 {
        return Err(DigestError::InvalidWindowSize);
    }

    unsafe {
        let c_seq = CString::new(seq).unwrap();
        let char_count = seq.chars().count();
        let mut result = Vec::with_capacity(char_count);
        let size = window_minimizer(c_seq.as_ptr(), char_count, k, window, result.as_mut_ptr());
        result.set_len(size);
        Ok(result)
    }
}

pub fn modimizer_rs(seq: &str, k: u32, mod_val: u32) -> Result<Vec<u32>, DigestError> {
    if k < 4 {
        return Err(DigestError::InvalidKmerSize);
    }
    if mod_val == 0 {
        return Err(DigestError::InvalidModValue);
    }

    unsafe {
        let c_seq = CString::new(seq).unwrap();
        let char_count = seq.chars().count();
        let mut result = Vec::with_capacity(char_count);
        let size = modimizer(c_seq.as_ptr(), char_count, k, mod_val, result.as_mut_ptr());
        result.set_len(size);
        Ok(result)
    }
}

pub fn syncmer_rs(seq: &str, k: u32, window: u32) -> Result<Vec<u32>, DigestError> {
    if k < 4 {
        return Err(DigestError::InvalidKmerSize);
    }
    if window == 0 {
        return Err(DigestError::InvalidWindowSize);
    }

    unsafe {
        let c_seq = CString::new(seq).unwrap();
        let char_count = seq.chars().count();
        let mut result = Vec::with_capacity(char_count);
        let size = syncmer(c_seq.as_ptr(), char_count, k, window, result.as_mut_ptr());
        result.set_len(size);
        Ok(result)
    }
} 