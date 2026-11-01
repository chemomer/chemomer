//! Layer emitters — the engine's extension axis (OCP). Each InChI
//! layer is its own module implementing [`LayerEmitter`]; the engine
//! runs its registered layers without knowing any of their internals.
//! Adding a layer (connection, mobile-H, stereo) is a new module plus
//! one registration in [`default_layers`].

pub mod charge;
pub mod connection;
pub mod fixedh;
pub mod formula;
pub mod hydrogen;
pub mod isotope;
pub mod proton;
pub mod standardization;
pub mod stereo;

use crate::error::Result;
use crate::graph::Structure;

/// InChI segment positions (standard InChI layer order). Registry
/// order is irrelevant; the assembler sorts by these.
pub mod order {
    pub const FORMULA: u8 = 0;
    pub const CONNECTION: u8 = 1;
    pub const H: u8 = 2;
    pub const Q: u8 = 3;
    pub const P: u8 = 4;
    pub const B: u8 = 5;
    pub const T: u8 = 6;
    pub const M: u8 = 7;
    pub const S: u8 = 8;
    pub const I: u8 = 9;
    pub const O: u8 = 10;
}

/// The layer kinds a registered set must cover before the engine may
/// assemble a full InChI string: a bare formula plus a connection
/// layer is the minimum that identifies anything.
pub const MANDATORY_KINDS: [Option<char>; 2] = [None, Some('c')];

pub trait LayerEmitter {
    /// The layer's name.
    fn name(&self) -> &'static str;

    /// The layer's position in the InChI segment order (see
    /// [`order`]). Registry order is irrelevant; the assembler sorts
    /// by this.
    fn order(&self) -> u8;

    /// The layer's segment identifier in the InChI string: `None`
    /// for the formula layer (the bare segment right after
    /// `InChI=1S/`), `Some('p')` for `/p`, `Some('c')` for `/c`, ...
    /// Gates match emitted segments to the oracle's recording BY
    /// THIS KIND, not by position — a milestone-partial engine emits
    /// a subset of the layers, and InChI orders them itself.
    fn segment_kind(&self) -> Option<char>;

    /// Emits the layer's segment, or `Ok(None)` when the structure
    /// carries nothing this layer records.
    fn emit(&self, structure: &Structure) -> Result<Option<String>>;
}

/// The engine's default layer set. `proton_shift` is the net proton
/// change the (DE)PROTONATION normalization applied (see
/// [`crate::canonical::protons`]); `options` ride into the layers
/// that consult them at emit (the stereo layers, per
/// [`crate::options::Options`]); the Fixed-H analysis (present under
/// `options.fixed_h`) feeds the full-InChI `/f` layer.
pub fn default_layers(
    options: &crate::options::Options,
    proton_shift: i32,
    fixed_h: Option<crate::canonical::fixedh::FixedH>,
) -> Vec<Box<dyn LayerEmitter>> {
    let mut layers: Vec<Box<dyn LayerEmitter>> = vec![
        Box::new(formula::FormulaLayer),
        Box::new(connection::ConnectionLayer),
        Box::new(hydrogen::HydrogenLayer),
        Box::new(charge::ChargeLayer),
        Box::new(proton::ProtonLayer::with_shift(proton_shift)),
        Box::new(stereo::BLayer::with_options(options.clone())),
        Box::new(stereo::TLayer::with_options(options.clone())),
        Box::new(stereo::MLayer::with_options(options.clone())),
        Box::new(stereo::SLayer::with_options(options.clone())),
        Box::new(isotope::IsotopeLayer),
    ];
    if let Some(fixed) = fixed_h {
        layers.push(Box::new(fixedh::FixedHLayer::new(fixed, options.clone())));
    }
    layers
}
