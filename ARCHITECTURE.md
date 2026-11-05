# chemomer — architecture

The Rust InChI engine. Design principles: **encapsulation** (every
concern behind a module boundary; callers see types, not column
arithmetic or rule tables), **open/closed** (new InChI layers,
invariant schemes, and input formats join by adding modules and
registering implementations — never by editing the engine core), and
**tests extracted from old InChI**: the oracle corpus in
[inchi-tests](https://github.com/chemomer/inchi-tests) is recorded
from the IUPAC reference implementation and gates this engine
byte-for-byte. The oracle is a development-time test oracle only —
nothing from the InChI Trust ships here.

## The derivation pipeline

```
molfile text ──▶ molfile::parse ──▶ graph::Structure ──▶ layers (registry) ──▶ segments
     (input edge)     (encapsulated)     (graph semantics)   (OCP axis)       (InChI order)
```

Each layer owns consecutive *decisions* of the derivation; the engine
only runs the registry and tags emitted segments by layer kind.
InChI's own segment order (formula, /c, /h, /q, /p, /i, /t…) is a
property of the format, not of the registry order — gates match
segments **by layer kind**, so a milestone-partial engine claims
exactly the layers it can decide.

## Module map

| Module | Responsibility | Extension axis |
|---|---|---|
| `error` | `Error` / `Result`; gate-visible messages. | — |
| `molfile` | V2000 subset reader → `Structure`. Column arithmetic and unsupported-construct reporting live here and nowhere else. | New input formats (SDF, SMILES) = new reader modules behind the same `Structure`. |
| `graph` | `Structure` / `AtomRecord`, bond orders, valences, implicit-H and composition semantics. Graph queries (bond-order sums, per-atom H) are graph concerns, not layer concerns. | New element/valence data lands in one table. |
| `layers` | One module per InChI layer; the `LayerEmitter` trait is the OCP seam: `name()`, `segment_kind()` (None = formula, Some('p')…), `emit(&Structure) -> Option<String>`. Registered: formula, connection, hydrogen (fixed + mobile groups), proton, stereo (`/t` `/m` `/s`), bond stereo (`/b`). | A new layer = a new module + one registration in `default_layers()`. |
| `layers::standardization` | The ±1 proton rules decoded from the oracle (anion/cation proton accounting) — shared by the layers the rules couple (formula adjusts H, proton emits /p). One home, no duplication. | M3 replaces its arithmetic with full mobile-H rules. |
| `canonical` | Native port of the reference's canonicalization paths (MIT source read at development time, never compiled or called): invariants → partition → equitable refinement → two-pass individualization search → connection table → `/c` DFS emission; `tautomer` ports the ichitaut endpoint semantics (valence arithmetic, donor/acceptor roles, the π-run walk grammar); `stereo` ports `set_stereo_atom_parity` (3D lift, H closure, triple product); `bond_stereo` the double-bond signs; `emit::parse_connection_segment` replays a segment back for the property test. | The seams stay corpus-arbitrated; new reference machinery joins as a sibling module. |
| `key` | The InChIKey: dependency-free SHA-256, the reference's split (major formula+c/h/q, proton segment, minor block doubled), the base-26 triplets/dublets (E** and TAA..TTV skipped). | The layout constants are the corpus-arbitrated seam. |
| `ffi` | The C ABI: owned C strings, a thread-local last error, `cdylib` build. The one module allowed `unsafe`. | New ABI functions join here; the engine surface they wrap never changes. |
| `wasm_binding` / `napi_binding` | Runtime bindings behind the `wasm` / `napi` features; the default build compiles none of them. The Python binding ships from chemomer/chemomer-py. | A new runtime = a new feature + module. |
| `engine` | The façade: `Engine::from_molfile`, `segments()`, `inchi_string()`, `inchikey()` — runs the registry, tags segments by kind. Closed for modification. | — |

## Encapsulation boundaries

- **Callers never see molfile column arithmetic** — `molfile::parse`
  is the only code that indexes a molfile line.
- **Layers never see each other.** The formula layer cannot consult
  the proton layer; coupled rules live in `standardization` and are
  consumed as plain values.
- **The engine never branches on layer identity.** It runs the
  registry and tags results; `segment_kind` carries the identity.
- Rule tables (valences, proton standardization) are private to the
  module that owns the semantics.

## Tests extracted from old InChI

The suite (inchi-tests) is the extracted behaviour of the old
InChI — every fixture is `(molfile, InChI, InChIKey)` recorded from
the IUPAC reference binary:

1. **Gate (`tests/gate.rs`)** — the whole-suite contract: for every
   fixture, every emitted segment must match the oracle's recording
   of that layer kind, byte-for-byte, and the native hasher must
   reproduce the recorded InChIKey. The gate compares every fixture
   and classifies the result — pass, expected-failure skip (triage,
   dated and decoded), RETIREMENT (a triaged fixture that passes
   blocks CI until its entry is withdrawn), isotopic, beyond-±1.
   No silent skips, no measurement overrides.
2. **Unit tests** — per module: molfile offsets, valence/H inference,
   Hill ordering, proton-rule cases (ammonium vs quaternary N vs
   carbanion), each decoded rule pinned by the fixture that
   discovered it.
3. **Extraction at scale** — `inchi-tests/scripts/import-pubchem-breadth.rb`
   pulls PubChem structures through the gem's molfile writer, verifies
   the oracle, and grows the corpus; disagreements become triage
   entries or bug fixes (this loop found four real identity bugs in
   its first run).

## Milestones → modules (all landed)

| Milestone | Modules | State |
|---|---|---|
| M0 corpus | inchi-tests (5559 fixtures at v0.3.1) | done |
| M1 formula + proton | `layers::formula`, `layers::proton`, `layers::standardization` | done — kind-matched byte agreement |
| M2 /c connection | `canonical` (invariants, refinement, ranking), `layers::connection` | done — disagreements triaged with decoded reasons |
| M3 mobile-H groups | `canonical::tautomer`, `layers::hydrogen` | done — 92 triage entries retired by the template slice; 64 dated entries remain |
| M4 stereo /t /b /m /s | `canonical::stereo`, `canonical::bond_stereo`, `layers::stereo` | done — the stereo (20/20) and E/Z (10/10) tranches gate |
| M5 InChIKey | `key` | done — every recorded key reproduces (native SHA-256) |
| M6 bindings | `ffi` (C ABI, cdylib) + the `wasm` / `napi` features (Python via chemomer-py) | done — default build compiles zero runtime dependencies |
| M7 engine claims | the Ruby `NativeEngine`, the site claims matrix | done — gem and matrix merged |
| Hardening | `emit::parse_connection_segment` + proptest, `benches/identity.rs`, `fuzz/` | done — fuzz runs manually on nightly |
