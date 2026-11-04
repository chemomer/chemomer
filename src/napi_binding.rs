//! The Node binding (feature `napi`): the engine for Node consumers
//! via napi-rs.

use napi_derive::napi;

/// Derives the identifier from a V2000 molfile.
#[napi]
pub fn from_molfile(molfile: String) -> napi::Result<String> {
    crate::Engine::from_molfile(&molfile)
        .and_then(|engine| engine.inchi_string())
        .map_err(|error| napi::Error::from_reason(error.message().to_string()))
}

/// Derives the AuxInfo line from a V2000 molfile.
#[napi]
pub fn aux_from_molfile(molfile: String) -> napi::Result<String> {
    crate::Engine::from_molfile(&molfile)
        .and_then(|engine| engine.aux_string())
        .map_err(|error| napi::Error::from_reason(error.message().to_string()))
}

/// Derives the InChIKey from a V2000 molfile.
#[napi]
pub fn inchikey(molfile: String) -> napi::Result<String> {
    crate::Engine::from_molfile(&molfile)
        .and_then(|engine| engine.inchikey())
        .map_err(|error| napi::Error::from_reason(error.message().to_string()))
}

#[cfg(test)]
mod tests {
    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn the_napi_surface_derives_the_water_identity() {
        // The napi wrappers keep the plain functions callable: the
        // native test drives the same entry points the JS glue does.
        assert_eq!(
            super::from_molfile(WATER.to_string()).expect("identifier"),
            "InChI=1S/H2O/h1H2"
        );
        assert_eq!(
            super::aux_from_molfile(WATER.to_string()).expect("aux"),
            "AuxInfo=1/0/N:1/rA:1nO/rB:/rC:12,12,0;"
        );
        assert_eq!(
            super::inchikey(WATER.to_string()).expect("key"),
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        );
    }

    #[test]
    fn napi_errors_carry_the_engine_message() {
        let error = super::from_molfile("not a molfile".to_string()).unwrap_err();
        assert!(!error.reason.is_empty());
    }
}
