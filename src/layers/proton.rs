//! The `/p` layer: the net proton shift the (DE)PROTONATION
//! normalization applied — protons removed from the submitted
//! structure read positive (added back), protons added read negative
//! — or, when no normalization fired, the ±1 heuristic marker.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::LayerEmitter;

#[derive(Default)]
pub struct ProtonLayer {
    shift: i32,
}

impl ProtonLayer {
    pub fn with_shift(shift: i32) -> Self {
        ProtonLayer { shift }
    }
}

impl LayerEmitter for ProtonLayer {
    fn name(&self) -> &'static str {
        "proton"
    }

    fn order(&self) -> u8 {
        crate::layers::order::P
    }

    fn segment_kind(&self) -> Option<char> {
        Some('p')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_typable_charges(structure)?;
        if self.shift != 0 {
            // The segment is the full InChI text: `p+1`, `p-1`, ...
            let sign = if self.shift > 0 { '+' } else { '-' };
            return Ok(Some(format!("p{sign}{}", self.shift.abs())));
        }
        Ok(standardization::proton_layer(structure))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;
    use crate::layers::LayerEmitter;

    fn charged_atom(element: &str, charge: i32, h_count: u32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            charge,
            h_pinned: false,
            isotope: 0,
        }
    }

    #[test]
    fn proton_layers_are_signed_segments() {
        let s = Structure {
            atoms: vec![charged_atom("N", 1, 4)],
            bonds: vec![],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(
            ProtonLayer::default().emit(&s).unwrap(),
            Some("p+1".to_string())
        );

        // A lone O-H anion is not carbon-carried: the parent path —
        // /p-1 with the parent water (fixture identity/hydroxide:
        // H2O /h1H2 /p-1). Contrast the carbon-carried methoxide,
        // which is fixed (/q-1).
        let s = Structure {
            atoms: vec![charged_atom("O", -1, 1)],
            bonds: vec![],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(
            ProtonLayer::default().emit(&s).unwrap(),
            Some("p-1".to_string())
        );
    }

    #[test]
    fn bond_charged_cations_emit_no_proton_layer() {
        // Quaternary N: four heavy bonds, no H — the + is the identity,
        // so no /p layer (fixture identity/tetramethylammonium).
        let s = Structure {
            atoms: vec![
                charged_atom("N", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(ProtonLayer::default().emit(&s).unwrap(), None);
    }

    #[test]
    fn all_heavy_onium_centers_keep_their_charge_in_q() {
        // Sulfonium (fixture identity/breadth3-814) and phosphonium
        // (identity/breadth3-500): the charge satisfies the neutral
        // valence exactly, no proton moved, so the oracle records
        // /q+1 with no /p layer.
        let s = Structure {
            atoms: vec![
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("S", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(ProtonLayer::default().emit(&s).unwrap(), None);

        let s = Structure {
            atoms: vec![
                charged_atom("P", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(ProtonLayer::default().emit(&s).unwrap(), None);
    }

    #[test]
    fn protonated_heteroatom_onium_centers_emit_p_plus_one() {
        // A group-15/16 center with an actual riding proton still
        // has it removable: R2SH+ and R3NH+ carry /p+1 (oracle
        // battery C[SH+](C), C[NH+](C)C).
        let s = Structure {
            atoms: vec![
                charged_atom("C", 0, 0),
                charged_atom("S", 1, 1),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(
            ProtonLayer::default().emit(&s).unwrap(),
            Some("p+1".to_string())
        );

        let s = Structure {
            atoms: vec![
                charged_atom("C", 0, 0),
                charged_atom("N", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1), (1, 3, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        assert_eq!(
            ProtonLayer::default().emit(&s).unwrap(),
            Some("p+1".to_string())
        );
    }

    #[test]
    fn onium_centers_carry_no_phantom_hydrogens() {
        // The formula and /h layers read the same table: a sulfonium
        // sulfur counts no implicit H (C4H11S, the neutral parent).
        let s = Structure {
            atoms: vec![
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
                charged_atom("S", 1, 0),
                charged_atom("C", 0, 0),
                charged_atom("C", 0, 0),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        };
        let hydrogens = s.hydrogens();
        assert_eq!(hydrogens[2], 0);
    }
}
