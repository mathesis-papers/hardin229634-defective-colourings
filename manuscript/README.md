# Extensions of Shrock's and Yurishchev's formulas to Hardin's defective-colouring recurrences

Haobo Ma (ChronoAI Pte Ltd, Singapore) and Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

## Contents and mathematical scope

- `note.tex` and `companions.tex`: the note and its complete written arguments.
- `verify.py`: the full last-row/budget identity for the five-column count, its generating function, coprimality modulo 101, and the two narrow-strip identities.
- `verify_siblings.py`: colour-and-reversal orbit transitions, integer characteristic polynomials, eleven full transfer-vector identities, generating functions, least eventual orders, and sharp starting indices.
- `checker-output.json` and `sibling-output.json`: exact-arithmetic values and polynomial certificates. Polynomial coefficients are in ascending order.
- `metadata.json`: title, creators, scope, license, and the five-column Lean declaration.

The literal width-three, width-four, and width-five recurrence statements are respectively `HardinA229632Strip.result`, `HardinA229633Strip.result`, and `HardinA229634Strip.result`, at commit `9849c0d81bfe0a953a10ad37400f5b1b61f028b5`. The two-row recurrence is `HardinA229638Strip.result` at the same commit. Their links and precise ranges are in the PDF. The generating functions, minimal orders, boundary residuals, general defect-budget theorem, and remaining companion results have the written proofs and exact finite identities specified in the note.

Hardin's empirical formulas are in OEIS A229634, Formula, revision 9, and A229637, Formula, revision 8, both dated 27 April 2021. The two-column and two-row formulas are derivations from earlier exact formulas of Shrock (2000) and Yurishchev (2001), respectively. The note gives their precise specializations and free-boundary scope. Additional empirical formulas are credited individually to Colin Barker and R. J. Mathar, and the one-row identities to N. J. A. Sloane, Barker, and Amiram Eldar. The zero-defect five-column companion A232917 explains the sextic base denominator. The bibliography distinguishes these concrete antecedents from the general strip-transfer method of Salas and Sokal.

## Rebuild the PDF

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard LaTeX installation with `amsmath`, `amssymb`, `amsthm`, `booktabs`, `array`, `geometry`, `lmodern`, and `hyperref`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded in `note.tex`; no BibTeX step is required.

## Check the finite identities

Python 3.9 or later, standard library only:

```sh
python3 verify.py --output checker-output.json
python3 verify_siblings.py --output sibling-output.json
```

The first command checks 729 coordinates and the second checks eleven orbit-vector identities. Every residual coordinate is zero. For every nonzero table sequence the reduced denominator is a cube of its zero-defect determinant, with coprimality certified modulo 101. The boundary residuals in Corollary 4 are nonzero.

Optional comparisons with the public OEIS tables:

```sh
curl -fL https://oeis.org/A229634/b229634.txt -o b229634.txt
curl -fL https://oeis.org/A229637/b229637.txt -o b229637.txt
python3 verify.py --table b229634.txt --output checker-output.json
python3 verify_siblings.py --table b229637.txt --output sibling-output.json
```

The table SHA256 values are `c866bceb77e2549b1eafd4d5e9f242d939c7ae8d102b6e43b91fec9b31a387a5` and `1979a56002878b5c7441bff4a8ed8c24217fd83f5d37db8067aeda60a45189be`, respectively. The comparisons cover all 210 five-column values and 227 of the 287 companion-table values; the latter are the entries within the supplied strip ranges. These optional commands reproduce the additional table-comparison fields in the included JSON outputs. The core identities and polynomials do not require downloading tables.

The folded-cube graph interpretations in A052482 and the other widths of A229606 are outside the note's coverage.

## License

Copyright 2026 Haobo Ma and Wenlin Zhang. This note and its accompanying original code are licensed under Creative Commons Attribution 4.0 International (CC BY 4.0, `cc-by-4.0`). The cited sources retain their respective attribution and terms.
