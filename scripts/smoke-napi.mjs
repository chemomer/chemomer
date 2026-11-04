// The napi surface smoke test: a V2000 molfile in, the water identity
// out — `npm test` drives this against the freshly built binding.
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const require = createRequire(import.meta.url);
const { fromMolfile, auxFromMolfile, inchikey } = require(
  fileURLToPath(new URL("../index.js", import.meta.url)),
);

const WATER = [
  "water",
  "  AsciiChem",
  "",
  "  1  0  0  0  0  0  0  0  0  0999 V2000",
  "   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0",
  "M  END",
  "",
].join("\n");

const inchi = fromMolfile(WATER);
const aux = auxFromMolfile(WATER);
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
  console.error("chemomer napi smoke test FAILED:");
  for (const failure of failures) console.error(`  ${failure}`);
  process.exit(1);
}
console.log(`identifier  ${inchi}`);
console.log(`aux         ${aux}`);
console.log(`key         ${key}`);
console.log("chemomer napi smoke test passed");
