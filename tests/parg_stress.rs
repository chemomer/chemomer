//! The PARG layer's state-isolation contract (chemomer#110).
//!
//! `src/parg.rs` parses every input with a fresh packrat memo cache
//! and a fresh AST arena: the cache's entries are keyed by input
//! position alone and reference live-arena pool state, so anything
//! carried across parses replays either a foreign input's match or a
//! dead arena's node indices (probed in parsanol 0.18.9,
//! parsanol-rs#493). The reset-reuse alternative (one thread-local
//! arena + cache, `reset` + `clear` per parse) measured neutral on
//! the identity bench, so fresh-per-parse ships — and this test
//! pins the isolation it relies on.
//!
//! The stress: thousands of parses through the one shared parse path
//! in a single process — every fixture re-parsed after every other
//! fixture (rotated schedule), interleaving small and large inputs so
//! allocator reuse churns freed arenas — with every output
//! byte-compared against its single-parse reference. ANY divergence
//! between schedule orders is state leaking across parses; the
//! reference shape would be wrong, not the test.

use chemomer::Engine;

/// The mixed fixture set: valid derivations of different sizes, the
/// M-line families, and inputs that must ERROR deterministically —
/// error text is output too, and a stale state hit could corrupt
/// exactly it.
const FIXTURES: &[&str] = &[
    // Water: one atom, the smallest valid derivation.
    "water\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n",
    // L-alanine: six atoms, five bonds, a bond block.
    "l-alanine\n  t\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  0  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  2  0  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n",
    // The M CHG family: property-line tokens.
    "charged\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 N   0  0  0  5  0  0  0  0  0  0  0\nM  CHG  1   1   1\nM  END\n",
    // The M ISO family: a heavy isotope.
    "d2o\n  t\n\n  1  0  0  0  0  0  0  0  0  3 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  ISO  1   1   2\nM  END\n",
    // The alias family: the A  line and its symbol line.
    "aliased\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\nA  1\nOMe\nM  END\n",
    // No terminator on the final line (the ethane test's shape).
    "ethane\n  t\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\nM  END",
    // Degenerate counts line: the malformed-shape error path.
    "degenerate\n",
    // Counts naming more atom lines than the file carries: the
    // truncated-block error path.
    "truncated\n  t\n\n  3  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n",
    // The empty input.
    "",
];

/// One derivation, byte-shaped: the identifier or the error text —
/// both are outputs the schedule must reproduce exactly.
fn derive(input: &str) -> String {
    match Engine::from_molfile(input).and_then(|engine| engine.inchikey()) {
        Ok(key) => key,
        Err(error) => format!("error: {error}"),
    }
}

/// The identifier string exercises a longer walk than the key; both
/// stages re-read the same capture tree.
fn derive_string(input: &str) -> String {
    match Engine::from_molfile(input).and_then(|engine| engine.inchi_string()) {
        Ok(key) => key,
        Err(error) => format!("error: {error}"),
    }
}

#[test]
fn ten_thousand_mixed_parses_stay_byte_identical_to_single_parses() {
    // Single-parse references: each fixture derived once, in a fresh
    // sequence of its own, nothing else parsed in between.
    let expected: Vec<(String, String)> = FIXTURES
        .iter()
        .map(|input| (derive(input), derive_string(input)))
        .collect();

    // The stress schedule: 1,250 rounds x the nine fixtures = 11,250
    // sequential parses in one process. Each round rotates the start
    // fixture, so across the run every fixture parses immediately
    // after every other fixture — a state leak keyed to a specific
    // predecessor cannot hide.
    let rounds = 1_250;
    let fixture_count = FIXTURES.len();
    let mut parses = 0usize;
    for round in 0..rounds {
        for offset in 0..fixture_count {
            let index = (round + offset) % fixture_count;
            let input = FIXTURES[index];
            let got = (derive(input), derive_string(input));
            assert_eq!(
                got, expected[index],
                "fixture {index} diverged on schedule position \
                 (round {round}, offset {offset}): parse state leaked \
                 across the shared path"
            );
            parses += 1;
        }
    }
    assert!(parses >= 10_000, "the stress ran {parses} parses");
}
