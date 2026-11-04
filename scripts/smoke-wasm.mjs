// The wasm surface smoke test: a V2000 molfile in, the water identity
// out — run against the wasm-pack output (`pkg/`, target nodejs) after
// `wasm-pack build --target nodejs --features wasm`.
import { from_molfile, aux_from_molfile, inchikey } from "../pkg/chemomer.js";

const WATER = [
  "water",
  "  AsciiChem",
  "",
  "  1  0  0  0  0  0  0  0  0  0999 V2000",
  "   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0",
  "M  END",
  "",
].join("\n");

const inchi = from_molfile(WATER);
const aux = aux_from_molfile(WATER);
const key = inchikey(WATER);

const failures = [];
if (inchi !== "InChI=1S/H2O/h1H2") {
  failures.push(`identifier: got ${JSON.stringify(inchi)}`);
}
if (aux !== "AuxInfo=1/0/N:1/rA:1nO/rB:/rC:12,12,0;") {
  failures.push(`aux: got ${JSON.stringify(aux)}`);
}
if (key !== "XLYOFNOQVPJJNP-UHFFFAOYSA-N") {
  failures.push(`key: got ${JSON.stringify(key)}`);
}
if (failures.length > 0) {
  console.error("chemomer wasm smoke test FAILED:");
  for (const failure of failures) console.error(`  ${failure}`);
  process.exit(1);
}
console.log(`identifier  ${inchi}`);
console.log(`aux         ${aux}`);
console.log(`key         ${key}`);
console.log("chemomer wasm smoke test passed");
