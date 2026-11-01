//! The tetrahedral stereo layers (`/t`, `/m`, `/s`) — M4 slice 1.
//! One derived segment set, three registered emitters: the InChI
//! layer order runs `/t`, then `/m`, then `/s`, each its own segment
//! kind. Input contract and parity model live in
//! [`crate::canonical::stereo`]; the generation options
//! ([`crate::options::Options`]) consult here at emit — `SNon`
//! withholds the whole set, `SRel`/`SRac` reshape the `/m`/`/s` pair.

use crate::canonical;
use crate::error::Result;
use crate::graph::Structure;
use crate::layers::{order, LayerEmitter};
use crate::options::Options;

/// Computes the stereo segments once per structure, under the
/// engine's options. The default path is the pre-options one,
/// bit-identical; `SNon` emits nothing (the perception's numbering
/// verdict keeps running elsewhere — the aux numbering still speaks
/// the stereo-informed canonicalization, matching the oracle's SNon
/// recordings).
fn segments_of(structure: &Structure, options: &Options) -> Vec<(char, String)> {
    if options.no_stereo {
        return Vec::new();
    }
    let counts = crate::layers::standardization::adjusted_hydrogens(structure);
    let (adjacency, invariants) = canonical::graph_of(structure);
    let mask = canonical::tautomer::endpoint_mask(structure, &counts);
    let no_taut = canonical::tautomer::without_mobile_h(&counts, &mask);
    let symmetry = canonical::base_symmetry(&adjacency, &invariants);
    let groups = canonical::group_membership(structure, &counts);
    let ranks = canonical::canonical_ranks(&adjacency, &invariants, &no_taut, &groups);
    canonical::stereo::segments(structure, &counts, &ranks, &symmetry)
}

macro_rules! stereo_layer {
    ($name:ident, $kind:expr, $position:expr, $doc:literal) => {
        #[doc = $doc]
        pub struct $name {
            options: Options,
        }

        impl $name {
            #[must_use]
            pub fn with_options(options: Options) -> Self {
                Self { options }
            }
        }

        impl LayerEmitter for $name {
            fn name(&self) -> &'static str {
                stringify!($name)
            }

            fn order(&self) -> u8 {
                $position
            }

            fn segment_kind(&self) -> Option<char> {
                Some($kind)
            }

            fn emit(&self, structure: &Structure) -> Result<Option<String>> {
                Ok(segments_of(structure, &self.options)
                    .into_iter()
                    .find(|(kind, _)| *kind == $kind)
                    .map(|(_, segment)| segment))
            }
        }
    };
}

stereo_layer!(TLayer, 't', order::T, "The `/t` stereocenter segment.");
stereo_layer!(MLayer, 'm', order::M, "The `/m` inversion-flag segment.");
stereo_layer!(SLayer, 's', order::S, "The `/s` stereo-type segment.");

/// The `/b` double-bond segment — its own derived set (the bond
/// stereo signs, not the tetrahedral centers); model in
/// [`crate::canonical::bond_stereo`]. `SNon` withholds it like the
/// tetrahedral set; the relative/racemic modes leave `/b` standing
/// (witnessed on the oracle: `-SRel` on an E/Z fixture keeps `/b`).
pub struct BLayer {
    options: Options,
}

impl BLayer {
    #[must_use]
    pub fn with_options(options: Options) -> Self {
        Self { options }
    }
}

impl LayerEmitter for BLayer {
    fn name(&self) -> &'static str {
        "BLayer"
    }

    fn order(&self) -> u8 {
        order::B
    }

    fn segment_kind(&self) -> Option<char> {
        Some('b')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        if self.options.no_stereo {
            return Ok(None);
        }
        let counts = crate::layers::standardization::adjusted_hydrogens(structure);
        let (adjacency, invariants) = canonical::graph_of(structure);
        let mask = canonical::tautomer::endpoint_mask(structure, &counts);
        let groups = canonical::group_membership(structure, &counts);
        let ranks = canonical::canonical_ranks(
            &adjacency,
            &invariants,
            &canonical::tautomer::without_mobile_h(&counts, &mask),
            &groups,
        );
        Ok(canonical::bond_stereo::segments(structure, &counts, &ranks)
            .into_iter()
            .next()
            .map(|(_, segment)| segment))
    }
}
