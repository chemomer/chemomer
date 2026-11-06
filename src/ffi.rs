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
}
