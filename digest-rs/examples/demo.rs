use digest::{window_minimizer_rs, modimizer_rs, syncmer_rs};

fn main() {
    // Longer test sequences
    let sequences = vec![
        "ACTGCTGACTACTAGCTAGTCGATGACTGCTGACTACTAGCTAGTCGATGAC",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA",
    ];

    // Test parameters
    let k_values = vec![4, 5];
    let window_sizes = vec![5, 10];
    let mod_values = vec![3, 5];

    println!("=== Window Minimizer Tests ===");
    for seq in &sequences {
        println!("\nSequence: {}", seq);
        for k in &k_values {
            for window in &window_sizes {
                println!("k={}, window={}", k, window);
                match window_minimizer_rs(seq, *k, *window) {
                    Ok(positions) => println!("  Positions: {:?}", positions),
                    Err(e) => println!("  Error: {}", e),
                }
            }
        }
    }

    println!("\n=== Modimizer Tests ===");
    for seq in &sequences {
        println!("\nSequence: {}", seq);
        for k in &k_values {
            for mod_val in &mod_values {
                println!("k={}, mod={}", k, mod_val);
                match modimizer_rs(seq, *k, *mod_val) {
                    Ok(positions) => println!("  Positions: {:?}", positions),
                    Err(e) => println!("  Error: {}", e),
                }
            }
        }
    }

    println!("\n=== Syncmer Tests ===");
    for seq in &sequences {
        println!("\nSequence: {}", seq);
        for k in &k_values {
            for window in &window_sizes {
                println!("k={}, window={}", k, window);
                match syncmer_rs(seq, *k, *window) {
                    Ok(positions) => println!("  Positions: {:?}", positions),
                    Err(e) => println!("  Error: {}", e),
                }
            }
        }
    }
} 