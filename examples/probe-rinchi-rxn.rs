//! Probe: parse an RXNfile through `parse_rxnfile` and print the five
//! output lines exactly as the reference `rinchi_cmdline` prints them
//! (RInChI, RAuxInfo, Long-/Short-/Web-RInChIKey). Development-time
//! diff harness against the built IUPAC oracle:
//!
//!   ./rinchi_cmdline reaction.rxn
//!   cargo run --example probe-rinchi-rxn -- reaction.rxn
//!
//! Line-by-line diff decides encode-path parity. `/force-equilibrium`
//! (set_directionality) mirrors the oracle flag of the same name.
use std::process::Command;

fn main() {
    let force = std::env::args().any(|arg| arg == "/force-equilibrium");
    let path = std::env::args()
        .nth(1)
        .filter(|arg| arg != "/force-equilibrium")
        .map(std::path::PathBuf::from)
        .expect("an RXNfile path argument");
    let ours = {
        let input = std::fs::read_to_string(&path).expect("read the reaction file");
        // The reference CLI's input detection: $RXN reads as an RXNfile,
        // a RInChI header as string input, everything else as RDfile.
        let first_line = input.lines().next().unwrap_or("").trim_end_matches('\r');
        let mut reaction = if first_line == "$RXN" {
            chemomer::parse_rxnfile(&input)
        } else {
            chemomer::parse_rdfile(&input)
        }
        .expect("chemomer parses");
        if force {
            reaction.set_directionality(chemomer::rinchi::Directionality::Equilibrium);
        }
        vec![
            reaction.rinchi_string().expect("rinchi"),
            reaction.rinchi_auxinfo().expect("rauxinfo"),
            reaction.rinchi_long_key().expect("long key"),
            reaction.rinchi_short_key().expect("short key"),
            reaction.rinchi_web_key().expect("web key"),
        ]
    };
    for line in &ours {
        println!("{line}");
    }

    let Ok(oracle) = std::env::var("RINCHI_CMDLINE") else {
        println!("cross-check: RINCHI_CMDLINE not set; oracle not compared");
        return;
    };
    let mut command = Command::new(&oracle);
    command.arg(&path);
    if force {
        command.arg("/force-equilibrium");
    }
    let output = command.output().expect("run the reference CLI");
    assert!(
        output.status.success(),
        "reference failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let theirs: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|line| {
            line.starts_with("RInChI=")
                || line.starts_with("RAuxInfo=")
                || line.starts_with("Long-RInChIKey=")
                || line.starts_with("Short-RInChIKey=")
                || line.starts_with("Web-RInChIKey=")
        })
        .map(str::to_owned)
        .collect();
    let mut differences = 0;
    for (index, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
        if ours == theirs {
            continue;
        }
        differences += 1;
        println!("--- line {index} differs");
        println!("ours:     {ours}");
        println!("oracle:   {theirs}");
    }
    if differences == 0 {
        println!("cross-check: all five lines match the oracle byte-for-byte");
    } else {
        println!("cross-check: {differences} line difference(s) against the oracle");
    }
}
