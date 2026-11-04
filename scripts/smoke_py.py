"""The pyo3 surface smoke test: a V2000 molfile in, the water identity
out — run against the built extension module (maturin develop, or a
wheel installed from dist/)."""

import sys

import chemomer

WATER = "\n".join(
    [
        "water",
        "  AsciiChem",
        "",
        "  1  0  0  0  0  0  0  0  0  0999 V2000",
        "   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0",
        "M  END",
        "",
    ]
)

inchi, key = chemomer.from_molfile(WATER)
aux = chemomer.aux_from_molfile(WATER)

failures = []
if inchi != "InChI=1S/H2O/h1H2":
    failures.append(f"identifier: got {inchi!r}")
if aux != "AuxInfo=1/0/N:1/rA:1nO/rB:/rC:12,12,0;":
    failures.append(f"aux: got {aux!r}")
if key != "XLYOFNOQVPJJNP-UHFFFAOYSA-N":
    failures.append(f"key: got {key!r}")
if failures:
    print("chemomer pyo3 smoke test FAILED:", file=sys.stderr)
    for failure in failures:
        print(f"  {failure}", file=sys.stderr)
    sys.exit(1)

print(f"identifier  {inchi}")
print(f"aux         {aux}")
print(f"key         {key}")
print("chemomer pyo3 smoke test passed")
