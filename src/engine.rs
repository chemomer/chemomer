//! The engine façade: parses a molfile and runs the registered
//! layers. Closed for modification, open for extension — new layers
//! join the registry in `layers::default_layers` without this file
//! changing.

use crate::error::{err, Result};
use crate::graph::Structure;
use crate::layers::{default_layers, standardization, LayerEmitter, MANDATORY_KINDS};
use crate::molfile;
use crate::options::Options;

pub struct Engine {
    structure: Structure,
    /// The submitted structure exactly as parsed — explicit hydrogens
    /// and input charges intact. The AuxInfo reversibility block
    /// echoes this, not the normalized form.
    original: Structure,
    /// The net proton shift the (DE)PROTONATION standardization and
    /// the hard group-H removal applied — the `/p` layer's value on
    /// the pinned path (zero when the layers read the heuristics).
    proton_shift: i32,
    /// The generation options this engine runs under — the stereo
    /// layers consult them at emit, the hydrogen perception obeys
    /// `do_not_add_h`, and the assembler reads the version prefix.
    options: Options,
    /// The Fixed-H analysis under `options.fixed_h`: the pre-taut
    /// snapshot's canonical numbering, H deltas, formula and charge
    /// (the /f layer and the aux /F: block read it).
    fixed_h: Option<crate::canonical::fixedh::FixedH>,
    layers: Vec<Box<dyn LayerEmitter>>,
}

impl Engine {
    /// Parses the V2000 molfile subset, normalizes the bonded ion
    /// pairs (terminal N⁺–O⁻), folds explicit hydrogen atoms onto
    /// their carriers, and equips the default layer registry — the
    /// standard identifier path.
    pub fn from_molfile(input: &str) -> Result<Self> {
        Self::with_options(input, Options::default())
    }

    /// [`Engine::from_molfile`] under generation `options` — the
    /// reference's option flags (see [`Options`]). The default options
    /// reproduce the standard identifier bit-for-bit.
    pub fn with_options(input: &str, options: Options) -> Result<Self> {
        let original = molfile::parse(input)?;
        let mut structure = original.clone();
        // The H-count codes as read (the reader maps MDL code n to
        // n-1 hydrogen): DoNotAddH subtracts them again after the
        // fold so only DRAWN hydrogens survive.
        let h_codes: Vec<u32> = structure.atoms.iter().map(|atom| atom.h_count).collect();
        standardization::normalize_internal_salts(&mut structure);
        structure.fold_explicit_hydrogens();
        if options.do_not_add_h {
            // DoNotAddH: only drawn hydrogens exist — every atom is
            // pinned at its drawn count (the folded explicit H atoms
            // minus the H-count codes, which the reference's
            // calculate_valences gate also ignores), so the
            // valence-gap inference never runs and the (DE)PROTONATION
            // standardization reads the drawn-only state from here on.
            for (atom, code) in structure.atoms.iter_mut().zip(&h_codes) {
                atom.h_count = atom.h_count.saturating_sub(*code);
                atom.h_pinned = true;
            }
        }
        // The Fixed-H SNAPSHOT (ichimake.c:3904): the preprocessed
        // structure copied BEFORE the taut pass — the (DE)PROTONATION
        // and the hard group-H removal never touch it.
        let snapshot = options.fixed_h.then(|| structure.clone());
        // The (DE)PROTONATION normalization: charged structures read
        // as their conjugate form before any layer runs; the pinned
        // H/charge state is final.
        let mut proton_shift = 0;
        if let Some(std) = crate::canonical::protons::standardize(&structure) {
            proton_shift = std.removed_protons as i32 - std.added_protons as i32;
            structure.protons_hard_moved = std.hard_moved;
            for (atom, h) in structure.atoms.iter_mut().zip(&std.h) {
                atom.h_count = *h;
                atom.h_pinned = true;
            }
            for (atom, charge) in structure.atoms.iter_mut().zip(&std.charge) {
                atom.charge = *charge;
            }
            // The hard removal's π alternation: the bond orders the
            // route moved ride with the proton (the reference keeps
            // them — BNS_EF_SAVE_ALL), and the tautomer machinery
            // reads the new orders from here on. The ride also leaves
            // its mobility MARKS (the reference's ALT12NS /
            // ALTERN|TRIPLE bond types): the salt-role sets read
            // those, not the plain orders.
            let mut ride_marks = vec![0u8; structure.bonds.len()];
            for (a, b, order, mark) in &std.bonds {
                for (index, bond) in structure.bonds.iter_mut().enumerate() {
                    let (from, to) = (bond.0, bond.1);
                    if (from == *a && to == *b) || (from == *b && to == *a) {
                        bond.2 = *order;
                        ride_marks[index] = *mark;
                    }
                }
            }
            structure.ride_marks = ride_marks;
        }
        let group_removals = remove_group_protons(&mut structure);
        proton_shift += group_removals;
        // The group-H removal is the 5.2b hard pass's group route; a
        // successful removal sets the same normalization flag.
        if group_removals > 0 {
            structure.protons_hard_moved = true;
        }
        // SUCF resolves per structure (the reference's on-the-fly mode
        // mutation, runichi2.c:538–552): the molfile chiral flag set
        // forces ABSOLUTE — the chiral-flag mode stays and the stereo
        // type reads /s1 with the /m standing; the flag clear forces
        // RELATIVE — the reference sets REQ_MODE_RELATIVE_STEREO
        // (clearing the racemic bit), so the relative reshape plus the
        // InChI=1/ prefix apply.
        let mut options = options;
        if options.chiral_flag_stereo && !structure.chiral_flag {
            options = options.with_relative_stereo();
        }
        // The Fixed-H analysis (ichican2.c 6246-6400): deltas,
        // CanonGraph07 numbering seeded by (nSymmRankBase, NumHfixed),
        // the snapshot's formula and charge.
        let fixed_h = match &snapshot {
            Some(snapshot) => {
                let centers = fixed_h_centers(&structure);
                let main_inverted = fixed_h_main_inverted(&structure, &options);
                Some(crate::canonical::fixedh::analyze(
                    snapshot,
                    &structure,
                    &centers,
                    main_inverted,
                )?)
            }
            None => None,
        };
        Ok(Engine {
            structure,
            original,
            proton_shift,
            options: options.clone(),
            fixed_h: fixed_h.clone(),
            layers: default_layers(&options, proton_shift, fixed_h),
        })
    }

    /// The generation options this engine runs under.
    pub fn options(&self) -> &Options {
        &self.options
    }

    /// The normalized structure — diagnostics for the triage tools.
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// The net proton shift the standardization applied (the `/p`
    /// layer's pinned-path value).
    pub fn proton_shift(&self) -> i32 {
        self.proton_shift
    }

    /// The derived segments, tagged by layer kind — everything the
    /// registered layers emit. Milestones add layers; the gates match
    /// tagged segments against the oracle's recording by kind.
    pub fn segments(&self) -> Result<Vec<(Option<char>, String)>> {
        let mut out = Vec::new();
        for layer in &self.layers {
            if let Some(segment) = layer.emit(&self.structure)? {
                out.push((layer.segment_kind(), segment));
            }
        }
        Ok(out)
    }

    /// Assembles the full InChI string, segments sorted by declared
    /// layer order. Errors until the registered layers cover the
    /// mandatory kinds (formula + connection): a partial layer set
    /// must never present itself as an identifier.
    pub fn inchi_string(&self) -> Result<String> {
        let registered: Vec<Option<char>> = self.layers.iter().map(|l| l.segment_kind()).collect();
        for mandatory in &MANDATORY_KINDS {
            if !registered.contains(mandatory) {
                return err(format!(
                    "inchi_string: mandatory layer kind {:?} not registered yet — \
                     the engine's layer set is milestone-partial",
                    mandatory
                        .map(|c| c.to_string())
                        .unwrap_or_else(|| "formula".to_string())
                ));
            }
        }

        let mut sorted: Vec<&Box<dyn LayerEmitter>> = self.layers.iter().collect();
        sorted.sort_by_key(|layer| layer.order());
        // The version prefix: the standard flag is the assembler's
        // only contract with the options — the non-std modes
        // (SRel/SRac/SUCF/SUU/SLUUD) read `InChI=1/` (the reference's
        // bStdFormat), everything else stays `InChI=1S/`.
        let mut out = String::from(if self.options.is_standard() {
            "InChI=1S/"
        } else {
            "InChI=1/"
        });
        let mut appended = false;
        for layer in sorted {
            if let Some(segment) = layer.emit(&self.structure)? {
                if appended {
                    out.push('/');
                }
                out.push_str(&segment);
                appended = true;
            }
        }
        Ok(out)
    }

    /// The InChIKey over the assembled identifier — the truncated
    /// SHA-256 key, computed natively (see [`crate::key`]).
    pub fn inchikey(&self) -> Result<String> {
        crate::key::inchikey(&self.inchi_string()?)
    }

    /// The key output the reference CLI prints for the assembled
    /// identifier under the options' key flags (ichimain.c:1334–1412):
    /// `InChIKey=<key>`, each requested extension on its own
    /// `XHash1=`/`XHash2=` line (newline-joined); under `MergeHash` the
    /// single `InChIHash=<key>-<xtra1>-<xtra2>` line with a bare dash
    /// for an absent extension. Without the flags this is the bare
    /// `InChIKey=` line. The identifier itself never changes — see
    /// [`Options`].
    pub fn inchikey_output(&self) -> Result<String> {
        let (key, xtra1, xtra2) = crate::key::key_extensions(&self.inchi_string()?)?;
        let xtra1 = self.options.xhash1.then_some(xtra1);
        let xtra2 = self.options.xhash2.then_some(xtra2);
        Ok(if self.options.merge_hash {
            crate::key::merged_line(&key, xtra1.as_deref(), xtra2.as_deref())
        } else {
            let mut out = format!("InChIKey={key}");
            if let Some(x1) = &xtra1 {
                out.push_str("\nXHash1=");
                out.push_str(x1);
            }
            if let Some(x2) = &xtra2 {
                out.push_str("\nXHash2=");
                out.push_str(x2);
            }
            out
        })
    }

    /// The AuxInfo line beside the identifier: canonical numbering,
    /// equivalence classes, charges, and the echo of the submitted
    /// structure (see [`crate::aux`]). The options ride along — SNon
    /// clears the inverted-stereo block; the Fixed-H analysis appends
    /// the `/F:` repetition block when the identifier's /f carried
    /// non-isotopic content.
    pub fn aux_string(&self) -> Result<String> {
        let fixed_block = self
            .fixed_h
            .as_ref()
            .map(|fixed| {
                crate::layers::fixedh::FixedHLayer::new(fixed.clone(), self.options.clone())
                    .non_iso_content(&self.structure)
                    .unwrap_or(false)
            })
            .unwrap_or(false);
        crate::aux::aux_info(
            &self.original,
            &self.structure,
            self.proton_shift,
            &self.options,
            self.fixed_h.as_ref(),
            fixed_block,
        )
    }
}

/// The wedge-defined stereo centers over the working structure — the
/// Fixed-H stereo refinement's input (the main pass's center inputs).
fn fixed_h_centers(working: &Structure) -> Vec<(usize, Vec<usize>, u32)> {
    let counts = standardization::adjusted_hydrogens(working);
    let ranks = crate::layers::hydrogen::canonical_ranks_for(working);
    crate::canonical::stereo::center_inputs(working, &counts, &ranks)
}

/// Whether the main aux verdict adopted its inverted pass (the aux
/// `/it:im` condition) — the /F: `it:m` marker mirrors it.
fn fixed_h_main_inverted(working: &Structure, options: &Options) -> bool {
    if options.no_stereo {
        return false;
    }
    let counts = standardization::adjusted_hydrogens(working);
    let ranks = crate::layers::hydrogen::canonical_ranks_for(working);
    let tautomeric = crate::canonical::tautomer::mobile_groups(working, &counts, &ranks)
        .iter()
        .any(|group| group.count > group.charge);
    crate::canonical::stereo::informed_verdict(working, &counts, &ranks, tautomeric)
        .map(|verdict| verdict.sign != 0 && crate::canonical::stereo::inverted_marker(&verdict))
        .unwrap_or(false)
}

/// The hard group-H removal (the reference's HardRemoveAcidicProtons,
/// TECHMAN 5.2b): a cationic structure whose registered mobile-H
/// group still holds a removable proton sheds it onto the group's
/// acidic oxygen as a fixed −1 — the cation stays put, the pinned
/// total reads neutral, and the shift surfaces as `/p+1` with no `/q`
/// and no mobile group. Each removal is gated by the reference's
/// REM_PROTON route test over the two virtual group stars (the
/// H-donor star and the acidic =O/N-O star): only a proton whose
/// donor-to-acceptor skeleton route the balanced search can cover
/// leaves the structure.
fn remove_group_protons(structure: &mut Structure) -> i32 {
    let mut removals = 0;
    loop {
        // The quota: one removal per unit of positive charge. The +2
        // family lands when the charged super-group merge follows.
        let total = structure.total_charge() + removals;
        if total <= 0 {
            break;
        }
        // The fixed cation the shift leaves behind: an N/P(+) with no
        // hydrogen of its own (the quaternary centers included).
        let plus = structure.atoms.iter().position(|atom| {
            atom.charge == 1
                && atom.h_count == 0
                && matches!(
                    crate::canonical::protons::el_number(&atom.element),
                    Some(7 | 15)
                )
        });
        let Some(plus) = plus else {
            break;
        };
        // The derived hydrogen counts (implicit H never shows in
        // `h_count`), the canonical ranks and the registered mobile
        // groups over them.
        let hydrogens = structure.hydrogens();
        let ranks = crate::layers::hydrogen::canonical_ranks_for(structure);
        let groups = crate::canonical::tautomer::mobile_groups(structure, &hydrogens, &ranks);
        let mut membership: Vec<Option<usize>> = vec![None; structure.atoms.len()];
        for (index, group) in groups.iter().enumerate() {
            for &rank in &group.ranks {
                if let Some(atom) = ranks.iter().position(|&r| r == rank) {
                    membership[atom] = Some(index);
                }
            }
        }
        // The holder: the group member still carrying its hydrogen.
        let Some(holder) = (0..structure.atoms.len())
            .find(|&atom| membership[atom].is_some() && hydrogens[atom] > 0)
        else {
            break;
        };
        // The acceptor: the group's double-bonded, hydrogen-free
        // member — where the moved proton turns into a fixed charge.
        let sums = structure.bond_order_sums();
        let degrees = structure.degrees();
        let members: Vec<usize> = groups[membership[holder].unwrap()]
            .ranks
            .iter()
            .filter_map(|&rank| ranks.iter().position(|&r| r == rank))
            .collect();
        let acceptor = members
            .iter()
            .copied()
            .find(|&atom| hydrogens[atom] == 0 && sums[atom] > degrees[atom] as f64)
            .unwrap_or(holder);
        let _ = plus;
        // The route: donor star vs acidic star over the whole
        // structure — no route, no removal.
        let route = rem_proton_route_for(structure, holder, acceptor);
        if std::env::var_os("INCHI_TRACE_PROTON").is_some() {
            eprintln!("GROUP-REM plus={plus} holder={holder} acceptor={acceptor} route={route}");
        }
        if !route {
            break;
        }
        // Apply: freeze the derived counts (the ±1 heuristics must not
        // re-read the half-moved state), then move the proton.
        if !structure.atoms.iter().all(|atom| atom.h_pinned) {
            for (atom, h) in structure.atoms.iter_mut().zip(&hydrogens) {
                atom.h_count = *h;
                atom.h_pinned = true;
            }
        }
        structure.atoms[holder].h_count = hydrogens[holder] - 1;
        structure.atoms[acceptor].charge -= 1;
        removals += 1;
    }
    removals
}

/// The REM_PROTON route test over a fresh network for this pass: one
/// t-group star over the H-donor heteroatoms, one over the acidic
/// =O/N-O oxygens, both as the reference builds them for the hard
/// removal (AR_HARD classes, neutral-valence gate, member-edge
/// capacity raises).
fn rem_proton_route_for(structure: &Structure, holder: usize, acceptor: usize) -> bool {
    let context = crate::canonical::protons::ProtonContext::new(structure);
    let (donors, acids) = crate::canonical::protons::hard_removal_stars(&context);
    crate::canonical::protons::rem_proton_route(&context, &donors, &acids, holder, acceptor).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Structure;
    use crate::layers::order;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    struct FixedSegment {
        order: u8,
        kind: Option<char>,
        text: &'static str,
    }

    impl LayerEmitter for FixedSegment {
        fn name(&self) -> &'static str {
            "fixed"
        }

        fn order(&self) -> u8 {
            self.order
        }

        fn segment_kind(&self) -> Option<char> {
            self.kind
        }

        fn emit(&self, _structure: &Structure) -> Result<Option<String>> {
            Ok(Some(self.text.to_string()))
        }
    }

    #[test]
    fn a_milestone_partial_layer_set_refuses_to_assemble_an_inchi() {
        // No connection layer registered: the set is milestone-partial
        // and must not present a prefix as an InChI.
        let engine = Engine {
            structure: crate::parse_molfile(WATER).expect("parses"),
            original: Structure::default(),
            proton_shift: 0,
            options: Options::default(),
            fixed_h: None,
            layers: vec![Box::new(FixedSegment {
                order: order::FORMULA,
                kind: None,
                text: "H2O",
            })],
        };
        let error = engine.inchi_string().unwrap_err();
        assert!(error.message().contains("mandatory layer kind"));
    }

    #[test]
    fn water_assembles_from_the_registered_layers() {
        let engine = Engine::from_molfile(WATER).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn the_key_output_mirrors_the_oracle_cli_lines() {
        // ichimain.c: -Key prints InChIKey=; each XHash flag adds its
        // line; -MergeHash collapses everything onto one InChIHash=
        // line with a bare dash for an absent extension. The oracle's
        // recorded water output (clean 1.07.5 binary).
        let x1 = "e014f429aafac730dc565eecf4ec4d2c4263e05c04ae2ef6";
        let x2 = "80fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
        let engine = Engine::from_molfile(WATER).expect("parses");
        assert_eq!(
            engine.inchikey_output().unwrap(),
            "InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        );
        let engine = Engine::with_options(WATER, Options::default().with_xhash1().with_xhash2())
            .expect("parses");
        assert_eq!(
            engine.inchikey_output().unwrap(),
            format!("InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N\nXHash1={x1}\nXHash2={x2}")
        );
        let engine =
            Engine::with_options(WATER, Options::default().with_xhash1().with_merge_hash())
                .expect("parses");
        assert_eq!(
            engine.inchikey_output().unwrap(),
            format!("InChIHash=XLYOFNOQVPJJNP-UHFFFAOYSA-N-{x1}-")
        );
    }

    const H2: &str = "h2\n  chemomer\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    1.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\nM  END\n";

    const LONE_H: &str = "h\n  chemomer\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const HYDRIDE: &str = "h-\n  chemomer\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  CHG  1   1  -1\nM  END\n";

    #[test]
    fn molecular_hydrogen_coalesces_to_the_oracle_identifier() {
        // The oracle reads H-H as one atom carrying one implicit H
        // (InChI=1S/H2/h1H); the old fold removed BOTH hydrogens and
        // underflowed the empty partition.
        let engine = Engine::from_molfile(H2).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2/h1H");
        assert_eq!(
            engine.aux_string().unwrap(),
            "AuxInfo=1/0/N:1/rA:2nHH/rB:s1;/rC:1,1,0;2,1,0;"
        );
    }

    #[test]
    fn a_lone_uncharged_hydrogen_reads_as_h2() {
        // The valence gap supplies its implicit partner.
        let engine = Engine::from_molfile(LONE_H).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2/h1H");
    }

    #[test]
    fn hydride_stays_a_bare_charged_atom() {
        let engine = Engine::from_molfile(HYDRIDE).expect("parses");
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H/q-1");
    }

    #[test]
    fn assembly_follows_declared_layer_order_not_registry_order() {
        // Registry scrambled on purpose; the assembler must sort by
        // `order()` and join only the segments that emit.
        let engine = Engine {
            structure: Structure::default(),
            original: Structure::default(),
            proton_shift: 0,
            options: Options::default(),
            fixed_h: None,
            layers: vec![
                Box::new(FixedSegment {
                    order: order::P,
                    kind: Some('p'),
                    text: "p+1",
                }),
                Box::new(FixedSegment {
                    order: order::CONNECTION,
                    kind: Some('c'),
                    text: "c1-2",
                }),
                Box::new(FixedSegment {
                    order: order::FORMULA,
                    kind: None,
                    text: "H2O",
                }),
            ],
        };
        assert_eq!(engine.inchi_string().unwrap(), "InChI=1S/H2O/c1-2/p+1");
    }
}
