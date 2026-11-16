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
| `molfile.parg` (the MDL V2000 molfile scan, v1.1.0) | `77755b0dc1fb38c16aeb32638b4d1184b588419125f3aea6c55088108a595132` |
| `molfile.json` (the compiled artifact) | `9af4ca187dccbe617bc5d022180ad93581a3cb20f008114ba4a57bf61d697cf9` |
| `inchi.parg` (the standard-InChI layer scan) | `c1d0a68c6b888b086e75fe559174f5e22b7d35f88a49a2166f53c28b95c888c8` |
| `inchi.json` (the compiled artifact) | `e5121fbea692c8a0e7be2cab468026c798caff0ce17566eb29fc780d836ea2ed` |
| `rinchi.parg` (the RInChI string scan) | `e38f0d9dc65f5d8a43265a16436cffaf051e954eca8e7f4867a43124a21cbc7d` |
| `rinchi.json` (the compiled artifact) | `bca7d7290d61e5deac6ecac6a4273cf1c97d0f4008e0ee6773c2b67a99487d76` |
| `smiles.parg` (the day-1 SMILES scan, v1.0.0) | `522e868047f4460faadf319e4f3e53b9e803a8ce0d1b9c33960b9c5abcafabf0` |
| `smiles.json` (the compiled artifact) | `ad5fb4709277895560033d101e55747aa3de1aa9972332594f08f0e1027431bf` |

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
| `molfile.parg` | the V2000 counts line, atom block, bond block, property tail (the alias, M CHG/M ISO and M Sgroup families) | the V3000 CTAB flavor, the SDfile record framing, the charge-code table, the polymer semantics |
| `inchi.parg` | the InChI=1S/ prefix, the Hill-order formula, the /c /h /q /p /i /b /t /m /s token streams, the exchanged-H record | the branch-stack bond walk, the atom-range/order checks, the element table, the reconstruct stage |
| `rinchi.parg` | the RInChI/RAuxInfo headers, the `<>` role groups, the `!` components, the /u //d trailing tags | the identifier validation (the key machinery), the reaction framing (RXN/RD files), every emitter |

Extracting a chemomer-grammar contract repository (the
asciichem-grammar pattern) can follow; the artifact envelope's
`source` + `source_hash` fields already carry the contract text.
