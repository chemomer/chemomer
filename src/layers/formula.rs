//! The `/formula` layer (M1): Hill order with the InChI twist
//! (carbon first when present, then hydrogen, then alphabetical;
//! counts omitted when 1), with the anion/cation hydrogen
//! standardization applied via the shared rules.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::LayerEmitter;

pub struct FormulaLayer;

impl LayerEmitter for FormulaLayer {
    fn name(&self) -> &'static str {
        "formula"
    }

    fn order(&self) -> u8 {
        crate::layers::order::FORMULA
    }

    fn segment_kind(&self) -> Option<char> {
        None
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_typable_charges(structure)?;
        let mut composition = structure.composition();
        let adjustment = standardization::hydrogen_adjustment(structure);
        if adjustment != 0 {
            add_hydrogen(&mut composition, adjustment);
        }
        Ok(Some(hill_string(&composition)))
    }
}

fn add_hydrogen(composition: &mut Vec<(String, u32)>, delta: i32) {
    for entry in composition.iter_mut() {
        if entry.0 == "H" {
            entry.1 = (entry.1 as i64 + delta as i64).max(0) as u32;
            return;
        }
    }
    if delta > 0 {
        composition.push(("H".to_string(), delta as u32));
    }
}

// Hill order with the InChI twist: carbon first (when present), then
// hydrogen, then the rest alphabetically; counts omitted when 1.
// Without carbon: pure alphabetical.
fn hill_string(composition: &[(String, u32)]) -> String {
    let has_carbon = composition.iter().any(|(e, _)| e == "C");
    let rank = |element: &str| match (has_carbon, element) {
        (true, "C") => 0,
        (true, "H") => 1,
        (true, _) => 2,
        (false, _) => 0,
    };
    let mut entries: Vec<&(String, u32)> = composition.iter().collect();
    entries.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)).then(a.0.cmp(&b.0)));
    entries
        .iter()
        .filter(|(_, count)| *count > 0)
        .map(|(element, count)| {
            if *count > 1 {
                format!("{element}{count}")
            } else {
                element.clone()
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;
    use crate::layers::LayerEmitter;

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
    }

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            h_pinned: false,
            valence_code: 0,
            isotope: 0,
            charge: 0,
        }
    }

    #[test]
    fn hill_orders_carbon_hydrogen_then_alphabetical() {
        assert_eq!(
            hill_string(&[("N".into(), 1), ("C".into(), 2), ("H".into(), 6)]),
            "C2H6N"
        );
        // Without carbon: pure alphabetical (water is H2O, not OH2).
        assert_eq!(hill_string(&[("O".into(), 1), ("H".into(), 2)]), "H2O");
    }

    #[test]
    fn counts_are_omitted_when_one() {
        assert_eq!(hill_string(&[("C".into(), 1), ("H".into(), 4)]), "CH4");
    }

    #[test]
    fn ammonium_loses_its_proton() {
        // [NH4+]: code-forced 4 H on the N -> H3N with /p+1 (the
        // fixture identity/ammonium).
        let mut n = atom("N");
        n.h_count = 4;
        n.charge = 1;
        let s = structure(vec![n], vec![]);
        assert_eq!(FormulaLayer.emit(&s).unwrap(), Some("H3N".to_string()));
    }

    #[test]
    fn acetate_gains_its_mobile_hydrogen() {
        // [CH3COO-]: the anionic O carries no H; the +1 H arrives via
        // standardization -> C2H4O2 (fixture identity/acetate).
        let mut o = atom("O");
        o.charge = -1;
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), o],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        assert_eq!(FormulaLayer.emit(&s).unwrap(), Some("C2H4O2".to_string()));
    }

    #[test]
    fn a_carbanion_keeps_its_composition() {
        // [CH3-]: the charge is on carbon — no proton accounting.
        let mut c = atom("C");
        c.h_count = 3;
        c.charge = -1;
        let s = structure(vec![c], vec![]);
        assert_eq!(FormulaLayer.emit(&s).unwrap(), Some("CH3".to_string()));
    }

    #[test]
    fn beyond_the_plus_minus_1_gate_is_rejected() {
        // The gate bounds PER-ATOM charges (GetAtomChargeType's
        // guard); multi-charge TOTALS (all atoms within ±1) assemble
        // through the pinned path.
        let mut o = atom("O");
        o.charge = -2;
        let s = structure(vec![o, atom("C")], vec![(0, 1, 1)]);
        let error = FormulaLayer.emit(&s).unwrap_err();
        assert!(error.message().contains("beyond the ±1 corpus gate"));
    }
}
