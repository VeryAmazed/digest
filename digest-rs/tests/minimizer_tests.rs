use digest::{window_minimizer_rs, modimizer_rs, syncmer_rs};

const TEST_SEQUENCES: &[(&str, &str)] = &[
    ("basic", "ACTGCTGACTACTAGCTAGTCGATGACTGCTGACTACTAGCTAGTCGATGACACTGCTGACTACTAGCTAGTCGATGAC"),
    ("homopolymer", "AAAAAAAAAAAAAAAAAAAAAAAACCCCCTTTTTTAAAAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAA"),
    ("repeating", "ACGTACGTACGTACGTACGTACGTACGTACGTATATATATATCGCGCGCGACGTACGTACGTACGTACGTACGTACGT"),
    ("with_n", "NACTGCTGACTACTAGNATCGATCGATCGANTCGATCGATCGANACTGCTGACTACTAGNATCGATCGATCGANTCGA"),
    ("complex", "ATCGATCGATCGATCGGGGGGATCGATCGATCGATCGAAAAATCGATCGATCGATCGTTTTTATCGATCGATCGATCG"),
    ("gc_rich", "GCGCGCGCGCGCGCGCGCATGCATGCGCGCGCGCGCGCGCGCTATATGCGCGCGCGCGCGCGCGCATGCATGCGCGC"),
    ("at_rich", "ATATATATATATATATATAGCTAGCTATATATATATATATAGCTAGCTATATATATATATATATATAGCTAGCTATA"),
];

#[test]
fn test_window_minimizer_basic() {
    let seq = TEST_SEQUENCES[0].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![4, 6, 8, 12, 14, 17, 21, 24, 29, 31, 33, 37, 39, 42, 46, 51, 52, 56, 58, 60, 64, 66, 69, 73]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![8, 12, 14, 21, 31, 33, 37, 39, 46, 56, 58, 60, 64, 66, 73]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![2, 3, 8, 10, 15, 17, 21, 25, 27, 28, 33, 35, 40, 42, 46, 48, 52, 54, 55, 60, 62, 67, 69, 73]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![2, 3, 10, 15, 17, 27, 28, 35, 40, 42, 52, 54, 55, 62, 67, 69]);
}

#[test]
fn test_window_minimizer_homopolymer() {
    let seq = TEST_SEQUENCES[1].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 26, 28, 32, 34, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 57, 58, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 32, 34, 44, 45, 46, 47, 48, 49, 50, 51, 53, 57, 58, 68, 69, 70, 71, 72, 73, 74, 75]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 28, 29, 30, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 29, 30, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74]);
}

#[test]
fn test_modimizer_basic() {
    let seq = TEST_SEQUENCES[0].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![0, 1, 6, 13, 17, 24, 25, 26, 31, 38, 42, 50, 51, 52, 53, 58, 65, 69]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![3, 7, 8, 10, 11, 12, 14, 15, 16, 19, 28, 32, 33, 35, 36, 37, 39, 40, 41, 44, 55, 59, 60, 62, 63, 64, 66, 67, 68, 71]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![1, 2, 7, 9, 10, 15, 17, 26, 27, 32, 34, 35, 40, 42, 50, 53, 54, 59, 61, 62, 67, 69]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![10, 15, 20, 35, 40, 45, 48, 49, 62, 67, 72]);
}

#[test]
fn test_modimizer_homopolymer() {
    let seq = TEST_SEQUENCES[1].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![21, 28, 32, 34, 51, 58]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![20, 21, 23, 25, 26, 52, 53, 55, 57]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![21, 22, 24, 26, 27, 28, 50, 51, 52, 54]);
}

#[test]
fn test_syncmer_basic() {
    let seq = TEST_SEQUENCES[0].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![0, 2, 4, 8, 10, 12, 14, 17, 21, 24, 25, 27, 29, 33, 35, 37, 39, 42, 46, 47, 48, 52, 54, 56, 60, 62, 64, 66, 69]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![3, 5, 12, 21, 22, 24, 28, 30, 37, 46, 47, 49, 51, 55, 57, 64]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![2, 3, 4, 6, 10, 11, 13, 17, 21, 23, 27, 28, 29, 31, 35, 36, 38, 42, 46, 48, 50, 54, 55, 56, 58, 62, 63, 65, 69]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![2, 3, 6, 8, 17, 18, 27, 28, 31, 33, 42, 43, 45, 54, 55, 58, 60]);
}

#[test]
fn test_syncmer_homopolymer() {
    let seq = TEST_SEQUENCES[1].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 22, 26, 28, 30, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 53, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 22, 23, 25, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 53, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65]);
}

#[test]
fn test_error_cases() {
    // Test invalid k-mer size
    assert!(window_minimizer_rs("ACGT", 2, 5).is_err());
    assert!(modimizer_rs("ACGT", 2, 3).is_err());
    assert!(syncmer_rs("ACGT", 2, 5).is_err());
    
    // Test invalid window/mod size
    assert!(window_minimizer_rs("ACGT", 4, 0).is_err());
    assert!(modimizer_rs("ACGT", 4, 0).is_err());
    assert!(syncmer_rs("ACGT", 4, 0).is_err());
}

#[test]
fn test_window_minimizer_repeating() {
    let seq = TEST_SEQUENCES[2].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 30, 31, 33, 35, 37, 40, 44, 46, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 30, 31, 33, 35, 37, 40, 46, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 29, 30, 35, 36, 37, 38, 40, 41, 46, 47, 49, 54, 57, 58, 61, 62, 65, 66, 69, 70, 73]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 38, 40, 41, 46, 47, 49, 58, 61, 62, 65, 66, 69, 70, 73]);
}

#[test]
fn test_window_minimizer_with_n() {
    let seq = TEST_SEQUENCES[3].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![5, 7, 9, 12, 21, 23, 25, 32, 34, 36, 38, 44, 48, 50, 52, 55, 64, 66, 68]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 12, 25, 32, 34, 36, 38, 44, 48, 50, 52, 55, 68]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![3, 4, 9, 11, 20, 23, 24, 32, 33, 36, 37, 44, 46, 47, 52, 54, 63, 66, 67]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![3, 4, 11, 24, 32, 33, 36, 37, 44, 46, 47, 54]);
}

#[test]
fn test_window_minimizer_complex() {
    let seq = TEST_SEQUENCES[4].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![4, 6, 8, 10, 12, 13, 14, 19, 20, 21, 23, 25, 27, 29, 31, 33, 36, 41, 43, 45, 47, 49, 51, 53, 54, 59, 64, 66, 68, 70, 72, 74]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![8, 10, 12, 21, 23, 25, 27, 29, 31, 33, 41, 43, 45, 47, 49, 51, 53, 59, 68, 70, 72, 74]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![3, 6, 7, 10, 11, 12, 13, 17, 20, 23, 24, 27, 28, 31, 32, 34, 37, 40, 44, 47, 48, 51, 52, 54, 57, 60, 65, 68, 69, 72, 73]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![7, 10, 11, 12, 13, 17, 27, 28, 31, 32, 34, 37, 40, 48, 51, 52, 54, 57, 60, 69, 72, 73]);
}

#[test]
fn test_window_minimizer_gc_rich() {
    let seq = TEST_SEQUENCES[5].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 29, 31, 33, 35, 37, 39, 40, 42, 44, 46, 50, 52, 54, 56, 58, 60, 62, 66, 70]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 11, 13, 15, 19, 23, 33, 35, 37, 39, 40, 42, 44, 46, 56, 58, 60, 62, 66]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 20, 21, 22, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 43, 45, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 67, 68, 69]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![9, 10, 11, 12, 13, 15, 22, 32, 33, 34, 35, 36, 37, 39, 43, 45, 55, 56, 57, 58, 59, 60, 62, 69]);
}

#[test]
fn test_window_minimizer_at_rich() {
    let seq = TEST_SEQUENCES[6].1;
    
    // k=4, w=5
    let result = window_minimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71]);
    
    // k=4, w=10
    let result = window_minimizer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 11, 13, 15, 17, 19, 21, 23, 33, 35, 37, 39, 41, 43, 45, 55, 57, 59, 61, 63, 65, 67, 69, 71]);
    
    // k=5, w=5
    let result = window_minimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 69, 70]);
    
    // k=5, w=10
    let result = window_minimizer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 69, 70]);
}

#[test]
fn test_modimizer_repeating() {
    let seq = TEST_SEQUENCES[2].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 25, 26, 27, 29, 39, 49, 51, 52, 53, 55, 56, 57, 59, 60, 61, 63, 64, 65, 67, 68, 69, 71, 72, 73]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![31, 32, 33, 34, 35, 36, 37, 38]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![29, 31, 32, 33, 34, 35, 36, 37, 42, 43, 44, 45, 48]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![]);
}

#[test]
fn test_modimizer_with_n() {
    let seq = TEST_SEQUENCES[3].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![1, 2, 7, 20, 24, 33, 37, 44, 45, 50, 63, 67]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![4, 8, 9, 11, 12, 18, 22, 26, 31, 35, 39, 47, 51, 52, 54, 55, 61, 65, 69, 74]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![2, 3, 8, 10, 11, 45, 46, 51, 53, 54]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![11, 54]);
}

#[test]
fn test_modimizer_complex() {
    let seq = TEST_SEQUENCES[4].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![3, 7, 11, 13, 19, 20, 24, 28, 32, 35, 36, 44, 48, 52, 55, 56, 59, 60, 61, 65, 69, 73]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![1, 5, 9, 22, 26, 30, 34, 35, 40, 42, 46, 50, 63, 67, 71]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![12, 14, 17, 19, 38, 39, 53, 54, 56, 60]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![12, 15, 16, 18, 54, 60]);
}

#[test]
fn test_modimizer_gc_rich() {
    let seq = TEST_SEQUENCES[5].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![15, 23, 46, 62, 70]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![15, 23, 40, 42, 43, 46, 62, 70]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 42, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 65, 66, 71, 72]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![40]);
}

#[test]
fn test_modimizer_at_rich() {
    let seq = TEST_SEQUENCES[6].1;
    
    // k=4, m=3
    let result = modimizer_rs(seq, 4, 3).unwrap();
    assert_eq!(result, vec![18, 22, 40, 44, 66, 70]);
    
    // k=4, m=5
    let result = modimizer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 67, 68, 69, 71, 73]);
    
    // k=5, m=3
    let result = modimizer_rs(seq, 5, 3).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]);
    
    // k=5, m=5
    let result = modimizer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![16, 23, 38, 45, 64, 71]);
}

#[test]
fn test_syncmer_repeating() {
    let seq = TEST_SEQUENCES[2].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 26, 27, 29, 31, 33, 35, 37, 40, 42, 44, 46, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 26, 28, 31, 33, 35, 37, 40, 42, 44, 46, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 41, 42, 43, 45, 49, 50, 53, 54, 57, 58, 61, 62, 65, 66, 69]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![0, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 23, 24, 27, 28, 29, 31, 32, 37, 38, 40, 49, 50, 52, 53, 54, 56, 57, 58, 60, 61, 62, 64]);
}

#[test]
fn test_syncmer_with_n() {
    let seq = TEST_SEQUENCES[3].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![1, 3, 5, 9, 12, 17, 19, 21, 23, 24, 25, 26, 32, 34, 36, 44, 46, 48, 52, 55, 60, 62, 64, 66]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![9, 12, 17, 19, 21, 23, 25, 31, 35, 37, 39, 52, 55, 60]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![3, 4, 5, 7, 11, 19, 20, 23, 24, 32, 33, 35, 37, 46, 47, 48, 50, 54, 62, 63]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![3, 4, 11, 18, 19, 20, 22, 23, 24, 25, 32, 46, 47, 54]);
}

#[test]
fn test_syncmer_complex() {
    let seq = TEST_SEQUENCES[4].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![0, 2, 4, 6, 8, 10, 12, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27, 29, 31, 33, 36, 37, 39, 41, 43, 45, 47, 49, 51, 53, 54, 55, 59, 60, 62, 64, 66, 68, 70]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 24, 25, 27, 29, 31, 32, 33, 34, 36, 38, 40, 41, 42, 43, 44, 45, 47, 49, 50, 59, 61, 62, 63, 64, 65]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![2, 3, 6, 7, 8, 9, 13, 17, 19, 20, 23, 24, 27, 28, 30, 34, 37, 40, 43, 44, 47, 48, 50, 53, 57, 60, 61, 64, 65, 68, 69]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![1, 2, 3, 4, 8, 17, 18, 19, 20, 22, 23, 24, 25, 34, 37, 40, 42, 43, 44, 45, 48, 57, 60, 63, 64]);
}

#[test]
fn test_syncmer_gc_rich() {
    let seq = TEST_SEQUENCES[5].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![1, 3, 5, 7, 9, 11, 15, 19, 23, 25, 27, 29, 31, 33, 35, 36, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 62, 66]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 10, 19, 23, 24, 25, 26, 27, 28, 29, 30, 31, 40, 42, 44, 46, 47, 48, 49, 50, 51, 52, 53, 57]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 39, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 62, 63, 64, 65]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 13, 15, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 60, 62]);
}

#[test]
fn test_syncmer_at_rich() {
    let seq = TEST_SEQUENCES[6].1;
    
    // k=4, w=5
    let result = syncmer_rs(seq, 4, 5).unwrap();
    assert_eq!(result, vec![1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]);
    
    // k=4, w=10
    let result = syncmer_rs(seq, 4, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 39, 41, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 60, 62]);
    
    // k=5, w=5
    let result = syncmer_rs(seq, 5, 5).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 65, 66]);
    
    // k=5, w=10
    let result = syncmer_rs(seq, 5, 10).unwrap();
    assert_eq!(result, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 35, 39, 40, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61]);
} 