//! Probe: decode a RInChI string and print the recovered roles —
//! `reactants[i] <inchi>` / `products[i] <inchi>` / `agents[i] <inchi>`
//! lines, then the directionality. With `RINCHI_CMDLINE` pointing at a
//! built reference `rinchi_cmdline` binary (the IUPAC oracle, read at
//! development time only, never vendored), the probe also cross-checks
//! our decode against the reference's own: its `/rxn` reconstruction
//! of the same string is reparsed and every reconstructed component
//! must land in the same role carrying the same standard InChI. The
//! reference's reconstruction drops the stereo layers (/b /t /m /s),
//! so components whose identifier carries stereo report a known
//! difference — the reference's own decode∘encode loses them too;
//! our identifier-carrying decode keeps the string exact.
use std::collections::BTreeMap;
use std::process::Command;

fn role_inchis(reaction: &chemomer::RinchiReaction) -> [BTreeMap<String, usize>; 3] {
    [
        counts_of(reaction.reactants()),
        counts_of(reaction.products()),
        counts_of(reaction.agents()),
    ]
}

fn counts_of(components: &[chemomer::RinchiComponent]) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for component in components {
        let inchi = component.inchi_string().expect("component inchi");
        *counts.entry(inchi).or_default() += 1;
    }
    counts
}

fn main() {
    let rinchi = std::env::args().nth(1).expect("a RInChI string argument");
    let reaction = chemomer::rinchi::reaction_from_rinchi(&rinchi).expect("decodes");
    let ours = role_inchis(&reaction);
    for (role, counts) in ["reactants", "products", "agents"].iter().zip(&ours) {
        for (inchi, count) in counts {
            println!("{role}[x{count}] {inchi}");
        }
    }
    println!("directionality: {:?}", reaction.directionality());

    let Ok(oracle) = std::env::var("RINCHI_CMDLINE") else {
        println!("cross-check: RINCHI_CMDLINE not set; reference decode not compared");
        return;
    };
    // The reference decode: its /rxn reconstruction of the string,
    // reparsed and re-identified through the engine.
    let dir = std::env::temp_dir().join("chemomer-rinchi-decode-probe");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let input = dir.join("input.rinchi");
    std::fs::write(&input, &rinchi).expect("write input");
    let output = Command::new(&oracle)
        .arg(&input)
        .arg("/rxn")
        .output()
        .expect("run the reference CLI");
    assert!(
        output.status.success(),
        "reference decode failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rxn = String::from_utf8_lossy(&output.stdout).to_string();
    let reference = chemomer::parse_rxnfile(&rxn).expect("reference RXN reparses");
    let theirs = role_inchis(&reference);

    let mut differences = 0;
    for ((role, ours), theirs) in ["reactants", "products", "agents"]
        .iter()
        .zip(&ours)
        .zip(&theirs)
    {
        if ours == theirs {
            continue;
        }
        for inchi in ours.keys().chain(theirs.keys()).collect::<Vec<_>>() {
            let ours_count = ours.get(inchi);
            let theirs_count = theirs.get(inchi);
            if ours_count != theirs_count {
                differences += 1;
                println!("diff {role}: ours {ours_count:?} / reference {theirs_count:?} — {inchi}");
            }
        }
    }
    if differences == 0 {
        println!("cross-check: all three roles match the reference decode exactly");
    } else {
        println!("cross-check: {differences} component difference(s) against the reference");
    }
}
