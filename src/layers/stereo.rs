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
/// recordings); `SRel`/`SRac` reshape the set per the reference's
/// GetSp3RelRacAbs: `/m` never prints, `/s` speaks `2`/`3`, and the
/// rows are the stereo-informed verdict's winning pass (the
/// lexicographically smaller of the absolute/inverted passes,
/// ichimak2.c:1224–1250). When the inversion verdict ties — the meso
/// or trivial-inversion case — no `/s` prints at all (SP3_ONLY), and
/// with no perceivable center the whole set stays empty (`SP3_NONE`:
/// the oracle records no `/s2` for a bond-only-stereo fixture under
/// `-SRel`).
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
    let emission = canonical::stereo::Emission {
        keep_undefined: options.include_unknown_stereo,
        label_unknown: options.distinct_unknown_undefined,
    };
    let segments = canonical::stereo::segments(structure, &counts, &ranks, &symmetry, emission);
    let stereo_type = options.requested_stereo_type();
    if stereo_type.is_none() {
        return segments;
    }
    // The relative/racemic reshape. The verdict mirrors the aux's
    // input assembly (tautomeric structures strip the mobile H and
    // drop the group rows) so the identifier and the aux speak one
    // verdict; no perceivable center means no set at all.
    let tautomeric = canonical::tautomer::mobile_groups(structure, &counts, &ranks)
        .iter()
        .any(|group| group.count > group.charge);
    let Some(verdict) = canonical::stereo::informed_verdict(structure, &counts, &ranks, tautomeric)
    else {
        return Vec::new();
    };
    if verdict.sign == 0 {
        // SP3_ONLY: the rows tie — the /t stands (the default path's
        // rendering is the winning one), no /m and no /s follow.
        return segments
            .into_iter()
            .filter(|(kind, _)| *kind == 't')
            .collect();
    }
    let mut out = Vec::new();
    if !verdict.main_pairs.is_empty() {
        // The winning pass's rows: the centers renumber per the
        // winning ranks, the parities speak 1 = odd / 2 = even, the
        // first center normalized `-` (the winning rows are the
        // lexicographic minimum, so its parity is odd).
        let first = verdict.main_pairs.first().map(|&(_, parity)| parity);
        let text = verdict
            .main_pairs
            .iter()
            .enumerate()
            .map(|(index, &(rank, parity))| {
                let mark = match first {
                    Some(first_parity) if index > 0 && parity != first_parity => '+',
                    _ => '-',
                };
                format!("{rank}{mark}")
            })
            .collect::<Vec<_>>()
            .join(",");
        out.push(('t', format!("t{text}")));
    }
    if let Some(kind) = stereo_type {
        out.push(('s', format!("s{kind}")));
    }
    out
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
        let emission = canonical::stereo::Emission {
            keep_undefined: self.options.include_unknown_stereo,
            label_unknown: self.options.distinct_unknown_undefined,
        };
        Ok(
            canonical::bond_stereo::segments(structure, &counts, &ranks, emission)
                .into_iter()
                .next()
                .map(|(_, segment)| segment),
        )
    }
}
