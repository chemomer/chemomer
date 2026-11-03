//! The `/c` layer (M2): the canonical connection segment.

use crate::canonical;
use crate::error::Result;
use crate::graph::Structure;
use crate::layers::{order, LayerEmitter};

pub struct ConnectionLayer;

impl LayerEmitter for ConnectionLayer {
    fn name(&self) -> &'static str {
        "connection"
    }

    fn order(&self) -> u8 {
        order::CONNECTION
    }

    fn segment_kind(&self) -> Option<char> {
        Some('c')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        if structure.atoms.len() < 2 {
            return Ok(None);
        }
        let (adjacency, invariants) = canonical::graph_of(structure);
        let counts = crate::layers::standardization::adjusted_hydrogens(structure);
        let mask = canonical::tautomer::endpoint_mask(structure, &counts);
        let groups = canonical::main_canon_groups(structure, &counts);
        let table = canonical::canonical_table(
            &adjacency,
            &invariants,
            &canonical::tautomer::without_mobile_h(&counts, &mask),
            &groups,
        );
        Ok(canonical::connection_segment(&table).map(|segment| format!("c{segment}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layers::LayerEmitter;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn single_atom_molecules_have_no_connection_layer() {
        let structure = crate::parse_molfile(WATER).expect("parses");
        let segments = ConnectionLayer.emit(&structure).expect("emits");
        assert_eq!(segments, None);
    }
}
