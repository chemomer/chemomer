# chemomer

The production InChI engine: a native Rust implementation of InChI
derivation — canonical ranking, every standard layer, and the
InChIKey — with no third-party binary at runtime.

It is gated by the shared oracle corpus
[inchi-tests](https://github.com/chemomer/inchi-tests), pinned in CI
at `v0.3.1` (5559 identity fixtures). The corpus is recorded from
IUPAC `inchi-1` (MIT) at
development time; the IUPAC binary is a test oracle only and is never
shipped or executed by this engine.

## Status

M0–M7 complete and merged. The full 1559-fixture corpus gates green:
formula, `/p` proton accounting, canonical ranking with `/c`, the
mobile-H groups, the tetrahedral layers (`/t` `/m` `/s`), double-bond
stereo (`/b`), and the InChIKey — asserted natively against every
recorded key (1559/1559). The 64 remaining disagreements are triaged
with dated, decoded reasons; every triage entry is a debt to retire,
not a settlement.

## Runtime surfaces

- **Library** — `Engine::from_molfile(molfile)` → `segments()`,
  `inchi_string()`, `inchikey()`. Zero runtime dependencies.
- **C ABI** — the crate builds a `cdylib`: `chemomer_from_molfile`,
  `chemomer_inchikey`, `chemomer_string_free`,
  `chemomer_last_error` — drivable from Ruby FFI, Python ctypes,
  Node ffi-napi, Swift.
- **Feature-gated bindings** — `wasm` (wasm-bindgen exports); the
  native `chemomer` Node module ships from
  [chemomer/chemomer-ts](https://github.com/chemomer/chemomer-ts), the
  Python module from
  [chemomer/chemomer-py](https://github.com/chemomer/chemomer-py).
  Each gates its own dependency tree; the default build compiles none
  of them.

The Ruby gem consumes the engine through the C ABI:
`AsciiChem::Inchi::NativeEngine` in
[asciichem-ruby](https://github.com/asciichem/asciichem-ruby).

## Commands

```sh
INCHI_TESTS=/path/to/inchi-tests cargo test   # unit tests + oracle gate + key gate
cargo bench                                   # criterion: identity pipeline
cargo +nightly fuzz run parse_molfile         # fuzz/ (manual, nightly)
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for the module map and
extension axes.
