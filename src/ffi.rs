//! The C ABI — M6. A dependency-free foreign-function surface so any
//! runtime (Ruby FFI, Python ctypes, Node ffi-napi, Swift, …) can
//! drive the engine without a bespoke binding. Ownership: every
//! returned string is heap-allocated and must be released with
//! [`chemomer_string_free`]; on error the function returns null and
//! [`chemomer_last_error`] carries the message for the calling
//! thread. The crate-wide `#![deny(unsafe_code)]` is lifted only
//! inside this module — the boundary where unsafety is the job.

use std::cell::RefCell;
use std::ffi::{c_char, CStr, CString};

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
}

fn set_error(message: String) {
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = CString::new(message).ok();
    });
}

/// Releases a string returned by this ABI.
///
/// # Safety
/// The pointer must have been returned by one of this ABI's string
/// functions and not yet freed.
#[no_mangle]
pub unsafe extern "C" fn chemomer_string_free(string: *mut c_char) {
    if string.is_null() {
        return;
    }
    drop(CString::from_raw(string));
}

/// The message of the last error on the calling thread, or null.
/// The storage stays owned by the thread until the next error.
#[no_mangle]
pub extern "C" fn chemomer_last_error() -> *const c_char {
    LAST_ERROR.with(|slot| {
        slot.borrow()
            .as_ref()
            .map(|message| message.as_ptr())
            .unwrap_or(std::ptr::null())
    })
}

/// Derives the identifier from a V2000 molfile. Null on error (see
/// [`chemomer_last_error`]).
///
/// # Safety
/// `molfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_from_molfile(molfile: *const c_char) -> *mut c_char {
    let identifier = derive(molfile, |engine| engine.inchi_string());
    match identifier {
        Ok(text) => match CString::new(text) {
            Ok(text) => text.into_raw(),
            Err(error) => {
                set_error(format!("identifier contains a NUL: {error}"));
                std::ptr::null_mut()
            }
        },
        Err(error) => {
            set_error(error.message().to_string());
            std::ptr::null_mut()
        }
    }
}

/// Derives the AuxInfo line from a V2000 molfile. Null on error.
///
/// # Safety
/// `molfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_aux_from_molfile(molfile: *const c_char) -> *mut c_char {
    let aux = derive(molfile, |engine| engine.aux_string());
    match aux {
        Ok(text) => match CString::new(text) {
            Ok(text) => text.into_raw(),
            Err(error) => {
                set_error(format!("aux contains a NUL: {error}"));
                std::ptr::null_mut()
            }
        },
        Err(error) => {
            set_error(error.message().to_string());
            std::ptr::null_mut()
        }
    }
}

/// Derives the InChIKey from a V2000 molfile. Null on error.
///
/// # Safety
/// `molfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_inchikey(molfile: *const c_char) -> *mut c_char {
    let key = derive(molfile, |engine| engine.inchikey());
    match key {
        Ok(key) => match CString::new(key) {
            Ok(key) => key.into_raw(),
            Err(error) => {
                set_error(format!("key contains a NUL: {error}"));
                std::ptr::null_mut()
            }
        },
        Err(error) => {
            set_error(error.message().to_string());
            std::ptr::null_mut()
        }
    }
}

/// The ABI version: 1.
#[no_mangle]
pub extern "C" fn chemomer_ffi_version() -> u32 {
    1
}

// -- the reaction-identity surface (the RInChI layer) -------------------
// The RXNfile and RDfile input edges derive the five oracle lines
// (identifier, RAuxInfo, the three key shapes); the decode direction
// reads a RInChI string back into its components as
// [`crate::RinchiReaction::component_table`] text. Same contract as
// the identity surface: null on error, the message on
// [`chemomer_last_error`], returned strings released with
// [`chemomer_string_free`].

/// Derives the RInChI identifier from a V2000 RXNfile. Null on error.
///
/// # Safety
/// `rxnfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_from_rxn(rxnfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rxnfile, crate::parse_rxnfile, |reaction| {
        reaction.rinchi_string()
    });
    exported(value, "rinchi")
}

/// Derives the RAuxInfo line from a V2000 RXNfile. Null on error.
///
/// # Safety
/// `rxnfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_aux_from_rxn(rxnfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rxnfile, crate::parse_rxnfile, |reaction| {
        reaction.rinchi_auxinfo()
    });
    exported(value, "rauxinfo")
}

/// Derives the short (hashed) RInChIKey from a V2000 RXNfile. Null on
/// error.
///
/// # Safety
/// `rxnfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_short_key_from_rxn(rxnfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rxnfile, crate::parse_rxnfile, |reaction| {
        reaction.rinchi_short_key()
    });
    exported(value, "rinchikey")
}

/// Derives the long RInChIKey from a V2000 RXNfile. Null on error.
///
/// # Safety
/// `rxnfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_long_key_from_rxn(rxnfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rxnfile, crate::parse_rxnfile, |reaction| {
        reaction.rinchi_long_key()
    });
    exported(value, "long key")
}

/// Derives the web RInChIKey from a V2000 RXNfile. Null on error.
///
/// # Safety
/// `rxnfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_web_key_from_rxn(rxnfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rxnfile, crate::parse_rxnfile, |reaction| {
        reaction.rinchi_web_key()
    });
    exported(value, "web key")
}

/// Derives the RInChI identifier from an MDL RDfile — the second
/// input edge, whose `$DATUM $MFMT` agent scan differs from the RXN
/// reader; the five assembly paths are the ones the RXN exports
/// carry, so the identifier line gates the reader's transport.
/// Null on error.
///
/// # Safety
/// `rdfile` must be a NUL-terminated UTF-8 C string.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_from_rd(rdfile: *const c_char) -> *mut c_char {
    let value = derive_reaction(rdfile, crate::parse_rdfile, |reaction| {
        reaction.rinchi_string()
    });
    exported(value, "rinchi")
}

/// Decodes a RInChI string (with its RAuxInfo companion; pass the
/// empty string when unknown) into the component table — the decode
/// direction's wire form documented on
/// [`crate::RinchiReaction::component_table`]. Null on error.
///
/// # Safety
/// `rinchi` and `rauxinfo` must be NUL-terminated UTF-8 C strings.
#[no_mangle]
pub unsafe extern "C" fn chemomer_rinchi_decode(
    rinchi: *const c_char,
    rauxinfo: *const c_char,
) -> *mut c_char {
    let value = decode_args(rinchi, rauxinfo)
        .and_then(|(rinchi, rauxinfo)| {
            crate::rinchi::reaction_from_rinchi_with_auxinfo(&rinchi, &rauxinfo)
        })
        .and_then(|reaction| reaction.component_table());
    exported(value, "component table")
}

/// Reads the two decode arguments and drives one assembler over the
/// decoded reaction — the two-string counterpart of [`derive`].
unsafe fn decode_args(
    rinchi: *const c_char,
    rauxinfo: *const c_char,
) -> crate::Result<(String, String)> {
    let read = |name: &str, pointer: *const c_char| -> crate::Result<String> {
        if pointer.is_null() {
            return crate::error::err(format!("{name} pointer is null"));
        }
        CStr::from_ptr(pointer)
            .to_str()
            .map(str::to_owned)
            .map_err(|_| crate::error::Error::new(format!("{name} is not valid UTF-8")))
    };
    let rinchi = read("rinchi", rinchi)?;
    let rauxinfo = read("rauxinfo", rauxinfo)?;
    Ok((rinchi, rauxinfo))
}

/// Reads the C string and drives the reaction reader plus one
/// assembler over it — the RXN/RD counterpart of [`derive`].
unsafe fn derive_reaction<T>(
    input: *const c_char,
    read: impl FnOnce(&str) -> crate::Result<crate::RinchiReaction>,
    f: impl FnOnce(&crate::RinchiReaction) -> crate::Result<T>,
) -> crate::Result<T> {
    if input.is_null() {
        return crate::error::err("reaction input pointer is null");
    }
    let text = CStr::from_ptr(input)
        .to_str()
        .map_err(|_| crate::error::Error::new("reaction input is not valid UTF-8"))?;
    let reaction = read(text)?;
    f(&reaction)
}

/// Materializes an assembled value as the returned C string under the
/// shared contract: null on error, the message on the thread-local
/// slot.
fn exported(value: crate::Result<String>, what: &str) -> *mut c_char {
    match value {
        Ok(text) => match CString::new(text) {
            Ok(text) => text.into_raw(),
            Err(error) => {
                set_error(format!("{what} contains a NUL: {error}"));
                std::ptr::null_mut()
            }
        },
        Err(error) => {
            set_error(error.message().to_string());
            std::ptr::null_mut()
        }
    }
}

unsafe fn derive<T>(
    molfile: *const c_char,
    f: impl FnOnce(&crate::Engine) -> crate::Result<T>,
) -> crate::Result<T> {
    if molfile.is_null() {
        return Err(crate::error::err::<()>("molfile pointer is null").unwrap_err());
    }
    let text = CStr::from_ptr(molfile)
        .to_str()
        .map_err(|_| crate::error::err::<()>("molfile is not valid UTF-8").unwrap_err())?;
    let engine = crate::Engine::from_molfile(text)?;
    f(&engine)
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    // The corpus's water-is-water RXNfile (inchi-tests v0.3.41,
    // rinchi/water-is-water): one water drawn each side.
    const WATER_RXN: &str = "$RXN\nwater-is-water\n  chemomer recording\n\n  1  1\n$MOL\n\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$MOL\n\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const WATER_RAUXINFO_GROUP: &str = "0/N:1/rA:1nO/rB:/rC:12,12,0;";

    #[test]
    fn the_abi_derives_identifier_and_key() {
        let molfile = CString::new(WATER).unwrap();
        unsafe {
            let identifier = chemomer_from_molfile(molfile.as_ptr());
            assert!(!identifier.is_null());
            let text = CStr::from_ptr(identifier).to_string_lossy().into_owned();
            assert_eq!(text, "InChI=1S/H2O/h1H2");
            chemomer_string_free(identifier);

            let key = chemomer_inchikey(molfile.as_ptr());
            let text = CStr::from_ptr(key).to_string_lossy().into_owned();
            assert_eq!(text, "XLYOFNOQVPJJNP-UHFFFAOYSA-N");
            chemomer_string_free(key);
        }
    }

    #[test]
    fn errors_surface_through_the_thread_local() {
        unsafe {
            let identifier = chemomer_from_molfile(std::ptr::null());
            assert!(identifier.is_null());
            let message = CStr::from_ptr(chemomer_last_error())
                .to_string_lossy()
                .into_owned();
            assert!(message.contains("null"), "{message}");
        }
    }

    #[test]
    fn invalid_molfiles_report_the_parse_error() {
        let molfile = CString::new("not a molfile").unwrap();
        unsafe {
            let identifier = chemomer_from_molfile(molfile.as_ptr());
            assert!(identifier.is_null());
            let message = CStr::from_ptr(chemomer_last_error())
                .to_string_lossy()
                .into_owned();
            assert!(message.contains("molfile"), "{message}");
        }
    }

    #[test]
    fn the_abi_derives_the_rinchi_lines_from_an_rxnfile() {
        let rxnfile = CString::new(WATER_RXN).unwrap();
        unsafe {
            let read = |pointer: *mut c_char| {
                assert!(!pointer.is_null());
                let text = CStr::from_ptr(pointer).to_string_lossy().into_owned();
                chemomer_string_free(pointer);
                text
            };
            assert_eq!(
                read(chemomer_rinchi_from_rxn(rxnfile.as_ptr())),
                "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+"
            );
            assert_eq!(
                read(chemomer_rinchi_aux_from_rxn(rxnfile.as_ptr())),
                format!("RAuxInfo=1.00.1/{WATER_RAUXINFO_GROUP}<>{WATER_RAUXINFO_GROUP}")
            );
            assert_eq!(
                read(chemomer_rinchi_short_key_from_rxn(rxnfile.as_ptr())),
                "Short-RInChIKey=SA-FUHFF-XLYOFNOQVP-XLYOFNOQVP-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ"
            );
            assert_eq!(
                read(chemomer_rinchi_long_key_from_rxn(rxnfile.as_ptr())),
                "Long-RInChIKey=SA-FUHFF-XLYOFNOQVPJJNP-UHFFFAOYSA-N--XLYOFNOQVPJJNP-UHFFFAOYSA-N"
            );
            assert_eq!(
                read(chemomer_rinchi_web_key_from_rxn(rxnfile.as_ptr())),
                "Web-RInChIKey=XLYOFNOQVPJJNPIXP-NUHFFFADPSCTJSA"
            );
        }
    }

    #[test]
    fn the_abi_decodes_a_rinchi_into_the_component_table() {
        let rinchi = CString::new("RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+").unwrap();
        let rauxinfo = CString::new(format!(
            "RAuxInfo=1.00.1/{WATER_RAUXINFO_GROUP}<>{WATER_RAUXINFO_GROUP}"
        ))
        .unwrap();
        let water = "InChI=1S/H2O/h1H2\tAuxInfo=1/0/N:1/rA:1nO/rB:/rC:12,12,0;\t\
                     XLYOFNOQVPJJNP-UHFFFAOYSA-N\n";
        unsafe {
            let table = chemomer_rinchi_decode(rinchi.as_ptr(), rauxinfo.as_ptr());
            assert!(!table.is_null());
            let text = CStr::from_ptr(table).to_string_lossy().into_owned();
            chemomer_string_free(table);
            assert_eq!(text, format!("+\n{water}<>\n{water}<>\n"));

            // The bare decode (empty companion) leaves the aux fields
            // empty; a null pointer reports through the thread-local
            // slot like every string export.
            let bare = chemomer_rinchi_decode(rinchi.as_ptr(), CString::new("").unwrap().as_ptr());
            let text = CStr::from_ptr(bare).to_string_lossy().into_owned();
            chemomer_string_free(bare);
            assert!(text.contains("InChI=1S/H2O/h1H2\t\tXLYOFNOQVPJJNP-UHFFFAOYSA-N"));

            let failed = chemomer_rinchi_decode(rinchi.as_ptr(), std::ptr::null());
            assert!(failed.is_null());
            let message = CStr::from_ptr(chemomer_last_error())
                .to_string_lossy()
                .into_owned();
            assert!(message.contains("rauxinfo"), "{message}");
        }
    }

    #[test]
    fn the_abi_reports_reaction_input_errors() {
        let rxnfile = CString::new("not an rxnfile").unwrap();
        unsafe {
            let identifier = chemomer_rinchi_from_rxn(rxnfile.as_ptr());
            assert!(identifier.is_null());
            let message = CStr::from_ptr(chemomer_last_error())
                .to_string_lossy()
                .into_owned();
            assert!(message.contains("$RXN"), "{message}");

            let identifier = chemomer_rinchi_from_rxn(std::ptr::null());
            assert!(identifier.is_null());
            let message = CStr::from_ptr(chemomer_last_error())
                .to_string_lossy()
                .into_owned();
            assert!(message.contains("null"), "{message}");
        }
    }
}
