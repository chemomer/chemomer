//! The napi link setup (feature `napi`): the Node-API symbols come
//! from the host process at load time, so the macOS cdylib link keeps
//! them undefined (the `napi-build` crate's `setup()` mac arm,
//! inlined here to keep the core dependency-free at build time too).
//! Linux shared objects and the Windows dyn-symbol path (built into
//! `napi-sys`) need no link arguments.

fn main() {
    if std::env::var_os("CARGO_FEATURE_NAPI").is_none() {
        return;
    }
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-cdylib-link-arg=-Wl,-undefined,dynamic_lookup");
    }
}
