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
| `molfile.parg` (the MDL V2000 molfile scan) | `13cc4b98561144259772da5555470b95adf1c44306377b6d00f651e0f8e1ee3c` |
| `molfile.json` (the compiled artifact) | `40325ac8ad9e852b4d6eab7294f8cc80fc991acf35be87de55d03b246a513407` |
| `inchi.parg` (the standard-InChI layer scan) | `c1d0a68c6b888b086e75fe559174f5e22b7d35f88a49a2166f53c28b95c888c8` |
| `inchi.json` (the compiled artifact) | `e5121fbea692c8a0e7be2cab468026c798caff0ce17566eb29fc780d836ea2ed` |
| `rinchi.parg` (the RInChI string scan) | `e38f0d9dc65f5d8a43265a16436cffaf051e954eca8e7f4867a43124a21cbc7d` |
| `rinchi.json` (the compiled artifact) | `bca7d7290d61e5deac6ecac6a4273cf1c97d0f4008e0ee6773c2b67a99487d76` |

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
