# Grammar artifact checksums

The compiled artifacts are committed HERE — this directory is the
single source of truth for both the PARG text contracts and the
compiled artifacts the engine embeds. `src/parg.rs` loads each
artifact through `PargArtifact::from_json`, which verifies the
artifact's own envelope checksum at every load, so a drifted artifact
refuses to parse rather than parsing differently. Re-verify the file
hashes below whenever a contract changes.

## The current pin

| File | sha256 |
|---|---|
| `molfile.parg` (the MDL molfile scan: V2000 + V3000 framing + SDfile records, v1.2.0) | `ea3b621f639cf7cde6829c0c7691edebc8193ea8076393edc3baec6f59cea48b` |
| `molfile.json` (the compiled artifact) | `accc25a5cd835c94d2d859bd9dc605f9ef7019c8fb63a0d03f987ed2a9a8fa64` |
| `inchi.parg` (the standard-InChI layer scan) | `c1d0a68c6b888b086e75fe559174f5e22b7d35f88a49a2166f53c28b95c888c8` |
| `inchi.json` (the compiled artifact) | `e5121fbea692c8a0e7be2cab468026c798caff0ce17566eb29fc780d836ea2ed` |
| `rinchi.parg` (the RInChI string scan) | `e38f0d9dc65f5d8a43265a16436cffaf051e954eca8e7f4867a43124a21cbc7d` |
| `rinchi.json` (the compiled artifact) | `bca7d7290d61e5deac6ecac6a4273cf1c97d0f4008e0ee6773c2b67a99487d76` |
| `smiles.parg` (the SMILES scan, v1.1.0) | `abc377bdeba856c93a7460092e090a18333b7355a03f37821d9315ccbf82ef26` |
| `smiles.json` (the compiled artifact) | `1c270d3f1ada1900b15c93678cd22e0f2a516f726484b5732c75a1308242a152` |

**molfile v1.1.0 (rs#110 round 2)** — the capture surface shrank to
what the engine reads: the header lines, the line-shape wrappers
(`as line` / `as body_line` / the per-shape names), the unread counts
granules (k6/k9/k15/k33/krest), the unread atom granules
(k20/k30/k34/k40/k51), the u0/u1 unit wrappers, and the per-character
`ch` captures are all gone. The consumed surface (molfile/counts/body,
the read granules, units/kt) is byte-identical; the grammar still
consumes every column verbatim (the inline suite is unchanged and
passes). The sample profile that drove this: malloc/free dominated the
water bench (~55%), the arena Hash stores per named capture and the
flatten walk's String keys were the bulk of it.

**smiles v1.1.0 (chemomer#131)** — the component dot: `entry smiles:
formula` wraps the chains — `chain 0*(("." as dot chain) as
component)` — and the reader resets its current atom at each
component boundary, so the next chain's first atom bonds to nothing
(the molfile multi-fragment semantics). The identity layer's
connected-component analysis needed ZERO new code: the 34
identity-multicomponent fixtures (mixtures, salts: `[Na+].[Cl-]`)
derive the recorded 1.07.5-oracle InChIs byte-for-byte through the
SMILES leg (tests/smiles_gate.rs rides the tranche). A trailing dot
and a leading dot still reject at parse. Recompiled with the pinned
1.3.94 compiler.

**smiles v1.0.0 (chemomer#121)** — the day-1 SMILES reader: chains,
branches, single-digit ring closures, organic + bracket atoms,
charges, isotopes (`as iso`), explicit H counts, aromatic markers.
The reader walks the capture tree into the same `graph::Structure`
the molfile path produces (isotopic H folded onto carriers at parse,
the pipeline shared verbatim through `Engine::from_read`). The gate:
**1,537 identity-fixture SMILES derive the recorded 1.07.5-oracle
InChIs byte-for-byte** (tests/smiles_gate.rs); 6 recorded skips — the
4 tetrahedral cases + 2 directional-bond cases (the contract's phase
2) and the 2 transcription disagreements (sulfuric-acid/sulfate-2:
the fixture smiles carry two S=O where the recorded molfile carries
one — the recording pipeline owns the fixture fix).

**molfile v1.2.0 (chemomer#122)** — the V3000 CTAB line framing and
the SDfile record framing join the contract. `v30_line` frames the
`M  V30 ` prefix structurally (a block line lacking it falls to
m_line / other_line, where the engine's prefix check reports it); the
`sdf` entry frames the `$$$$` terminators (`$$$$` + optional
whitespace + line end — the `trim_end() == "$$$$"` read — and the
`!`-guard against any character is the end-of-input anchor) with
record lines riding ordered k0 captures. The semantics stay
engine-side: the `-` continuation joins and tokenization of the
joined logical line, the COUNTS/atom/bond block walk, and the
empty-record filter. Both flavors' bodies are consumed through the
one `mol` entry's body shapes; the version dispatch reads the counts
capture's k33 granule. The order-dependence lints are unchanged from
v1.1.0's baseline (the body_line classification is ordered choice by
design).

**smiles v1.0.0, the repeated-sign charge spelling (chemomer#126)** —
the charge rule gains the `++` / `--` alternative (each repeat adds
one): real-world writers spell ferric `[Fe++]`, the walk's
parse_charge already folded the repeats, only the grammar rejected
them. The inline suite carries the accepts; the asciichem corpus's
`structure/smiles/ferric-plusplus-alt` joins the recorded identity
intersection on the next recording run.

**smiles v1.0.0, the phase-2 stereo ride (chemomer#128)** — the
contract text gains no productions (chirality and the directional
bonds were always captured); the header's phase note reflects the
carried state and both hashes refresh. Engine-side,
`synthesize_stereo` turns the captured markers into drawn geometry —
the parity atoms' ligands laid out trigonally with a wedge realizing
the @/@@ hand per the SMILES neighbor order, the directional bonds'
substituents lifted to the side their spelling draws — so the wedge
perception and the cis/trans side product read a SMILES structure
like a drawn molfile. The gate: the alanine enantiomers derive the
recorded oracle bytes and run IN the gate now (2 skips became runs);
the two butene rows stay skipped — their fixture molfiles were
recorded without coordinates, so the oracle dropped the E/Z layer
(both record the same stereo-less string); un-skipping them needs a
corpus re-record from coordinate-bearing molfiles (the /b signs are
oracle-verified over drawings: trans → +, cis → −).

## Compiling

The artifacts are produced by the `Parsanol::PARG` compiler
(the `parsanol` gem). The engine consumes them through
`parsanol = "0.18.9"`; the artifact checksum canonicalization the Rust
verifier reads matches the 1.3.94 compiler, so re-compile with that
version pinned:

```sh
gem install parsanol -v 1.3.94
parsanol _1.3.94_ parg test    grammar/molfile.parg   # the contract's inline suite
parsanol _1.3.94_ parg compile grammar/molfile.parg -o grammar/molfile.json
parsanol _1.3.94_ parg test    grammar/inchi.parg
parsanol _1.3.94_ parg compile grammar/inchi.parg   -o grammar/inchi.json
parsanol _1.3.94_ parg test    grammar/rinchi.parg
parsanol _1.3.94_ parg compile grammar/rinchi.parg  -o grammar/rinchi.json
```

Re-compiling and committing both files is part of any grammar change:
update the contract, run the contract's test block, re-compile,
re-hash, update this file, then ride the change with the engine's
gates (`cargo test` — the identity corpus re-reads every fixture
through the new artifact).

## Scope

| Contract | Scans | Stays engine-side |
|---|---|---|
| `molfile.parg` | the V2000 counts line, atom block, bond block, property tail (the alias, M CHG/M ISO and M Sgroup families); the V3000 `M  V30 ` line framing; the SDfile `$$$$` record framing | the V3000 continuation joins, the COUNTS/atom/bond block walk and the CTAB tail blocks, the empty-record filter, the charge-code table, the polymer semantics |
| `inchi.parg` | the InChI=1S/ prefix, the Hill-order formula, the /c /h /q /p /i /b /t /m /s token streams, the exchanged-H record | the branch-stack bond walk, the atom-range/order checks, the element table, the reconstruct stage |
| `rinchi.parg` | the RInChI/RAuxInfo headers, the `<>` role groups, the `!` components, the /u //d trailing tags | the identifier validation (the key machinery), the reaction framing (RXN/RD files), every emitter |

Extracting a chemomer-grammar contract repository (the
asciichem-grammar pattern) can follow; the artifact envelope's
`source` + `source_hash` fields already carry the contract text.
