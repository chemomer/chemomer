//! The `/q` layer: the fixed charges the normalization did not move —
//! an all-heavy onium center (R4N+, R3S+, R3O+, R4P+), a lone
//! alkoxide, a carbanion — as one signed total.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::LayerEmitter;

pub struct ChargeLayer;

impl LayerEmitter for ChargeLayer {
    fn name(&self) -> &'static str {
        "charge"
    }

    fn order(&self) -> u8 {
        crate::layers::order::Q
    }

    fn segment_kind(&self) -> Option<char> {
        Some('q')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_typable_charges(structure)?;
        Ok(standardization::charge_layer(structure))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

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
    fn all_heavy_onium_centers_emit_q() {
        // Tetramethylammonium: the + is the identity (fixture
        // identity/tetramethylammonium; oracle battery C[N+](C)(C)C).
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
            taut_ext: Default::default(),
            chiral_flag: false,
        };
        assert_eq!(ChargeLayer.emit(&s).unwrap(), Some("q+1".to_string()));
    }

    #[test]
    fn lone_alkoxide_and_carbanion_emit_negative_q() {
        // Methoxide: [O-]C — no mobile-H partner, so the anion is
        // fixed (oracle battery: q-1). Carbanion [CH3-]: same
        // (fixture identity/breadth-74).
        let s = Structure {
            atoms: vec![charged_atom("O", -1, 0), charged_atom("C", 0, 0)],
            bonds: vec![(0, 1, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            taut_ext: Default::default(),
            chiral_flag: false,
        };
        assert_eq!(ChargeLayer.emit(&s).unwrap(), Some("q-1".to_string()));

        let s = Structure {
            atoms: vec![charged_atom("C", -1, 0)],
            bonds: vec![],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            taut_ext: Default::default(),
            chiral_flag: false,
        };
        assert_eq!(ChargeLayer.emit(&s).unwrap(), Some("q-1".to_string()));
    }

    #[test]
    fn internal_salts_annihilate_before_the_layer_runs() {
        // Nitromethane C[N+](=O)[O-]: the pair cancels — no /q (and
        // the /h layer records fixed H only).
        let mut s = Structure {
            atoms: vec![
                charged_atom("C", 0, 0),
                charged_atom("N", 1, 0),
                charged_atom("O", 0, 0),
                charged_atom("O", -1, 0),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            taut_ext: Default::default(),
            chiral_flag: false,
        };
        standardization::normalize_internal_salts(&mut s);
        assert!(s.atoms.iter().all(|atom| atom.charge == 0));
        assert_eq!(ChargeLayer.emit(&s).unwrap(), None);
        assert_eq!(standardization::proton_layer(&s), None);
    }

    #[test]
    fn distant_ions_never_annihilate() {
        // Two charges far apart keep their classifications: the N+
        // is a proton-related pnicogen deficit (/p+1), the lone O-
        // is fixed (/q-1).
        let mut s = Structure {
            atoms: (0..8)
                .map(|i| {
                    if i == 0 {
                        charged_atom("N", 1, 0)
                    } else if i == 7 {
                        charged_atom("O", -1, 0)
                    } else {
                        charged_atom("C", 0, 0)
                    }
                })
                .collect(),
            bonds: (0..7).map(|i| (i, i + 1, 1)).collect(),
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            taut_ext: Default::default(),
            chiral_flag: false,
        };
        standardization::normalize_internal_salts(&mut s);
        assert_eq!(s.atoms[0].charge, 1);
        assert_eq!(s.atoms[7].charge, -1);
        assert_eq!(ChargeLayer.emit(&s).unwrap(), Some("q-1".to_string()));
        assert_eq!(standardization::proton_layer(&s), Some("p+1".to_string()));
    }
}
