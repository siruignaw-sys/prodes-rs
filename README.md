# prodes-rs

A from-scratch Rust reimplementation of [`prodes`](https://github.com/tneijenhuis/prodes), a Python tool that calculates features from 3D protein structures with a focus on the protein surface — solvent-accessible surface area, isoelectric point, and surface electrostatic potential.

`prodes-rs` computes the same quantities as the original Python tool, verified where possible against its actual output, but is not a line-by-line port. It's built on [`pdbtbx`](https://crates.io/crates/pdbtbx) for PDB parsing rather than reimplementing a parser from scratch, and uses a functional style over borrowed data rather than `prodes`'s mutable, object-oriented design. The goal has been numerical fidelity to `prodes`'s outputs, not structural fidelity to its internals — the two codebases solve the same problem in ways suited to their respective languages.

## Relationship to `prodes` and citation

`prodes` and this reimplementation exist to support the same underlying research: predicting protein retention in ion-exchange chromatography from structural features, developed by Neijenhuis et al. If you use `prodes-rs`, please cite the original tool and the paper it was built to support:

> Neijenhuis, T., Le Bussy, O., Geldhof, G., Klijn, M. E., & Ottens, M. (2024). Predicting protein retention in ion-exchange chromatography using an open source QSPR workflow. *Biotechnology Journal*, 19, e2300708. https://doi.org/10.1002/biot.202300708

The two-feature model (isoelectric point + sum of negative surface electrostatics) this tool is built to support comes from:

> Neijenhuis, T., Le Bussy, O., Geldhof, G., Klijn, M. E., & Ottens, M. (2025). Using generalized quantitative structure–property relationship (QSPR) models to predict host cell protein retention in ion-exchange chromatography. *Journal of Chemical Technology and Biotechnology*, 101, 1420–1428. https://doi.org/10.1002/jctb.70026

`prodes` itself is MIT licensed and requests citation of the first paper above; this project follows the same convention.

## Status

Implemented and tested:

- **Isoelectric point** (`core.rs`) — Henderson-Hasselbalch charge model over a standard ionizable-residue pKa table, including N/C-termini, alternate-location conformer deduplication, and correct handling of trailing HETATM records (e.g. crystallographic waters) that would otherwise be mistaken for a chain's true C-terminus. Multi-atom ionizable groups (e.g. arginine's guanidinium, aspartate's carboxylate) correctly split their group charge across their constituent atoms rather than over-counting it. Validated against the full 13-protein training set from Neijenhuis et al. (2025), Table 1: 11 of 13 proteins land within ~0.15 pH units of the paper's published theoretical pI, several within 0.01–0.05 (bovine serum albumin, PDB 4F5S: 5.494 computed vs. 5.5 published). The remaining discrepancies have identified, understood causes — see "Known differences from `prodes`" below.
- **Solvent-accessible surface area** (`calculations/sasa.rs`) — Shrake-Rupley algorithm with a spatial-hash grid for neighbor search. Sample point density per atom scales with each atom's own probe-inflated surface area, matching `prodes`'s approach rather than using a fixed point count. Cross-checked directly against a live run of the real Python `prodes` on the same structures: within ~1.5% of the reference point count on 1GDW, and within ~0.5% on 4F5S chain A (29,287 points vs. 29,112; 28,500 negative-potential points vs. 28,632).
- **Surface electrostatic potential** (`calculations/electrostatics.rs`, `calculations/distance_functions.rs`) — Coulomb potential summed at a coarsened grid of surface points (`calculations/geometry.rs::property_points_on_surface`), filtered to negative-potential points and summed, matching `prodes`'s `SurfEpNegSumAverage` feature. The underlying Coulomb formulas are verified bit-for-bit against the Python reference on matched inputs, and the aggregate feature value has been directly compared against a live `prodes` run on a real structure (see "Known differences from `prodes`" below for the one identified source of divergence).

Not yet implemented:

- The "Formal" (binarized ±1/n) charge mode — only the continuous Henderson-Hasselbalch ("Average") charge mode is implemented, which is what the two-feature retention model actually uses.
- `prodes`'s "Shell" electrostatics feature family (potential projected onto a plane with a two-dielectric protein/solvent boundary) — a separate, more elaborate feature set in the original tool that the two-feature model doesn't require.
- Hydrophobicity/lipophilicity features, surface shape descriptors, and custom pKa file import (PROPKA/H++/pypka) — all present in `prodes`, none currently ported.

This has been verified through unit tests on synthetic structures, hand-computed numeric checks against the underlying physics, and direct comparisons against live runs of the real Python `prodes` on real structures — not through an exhaustive parity test suite across `prodes`'s full feature set.

### Known differences from `prodes`

On residues with two full alternate-location conformers spanning their entire side chain — a real occurrence in disordered regions of crystal structures, not a hypothetical edge case — `prodes-rs` selects one conformer and assigns it the residue's correct charge. The reference Python `prodes`, on the same input, silently assigns zero charge to both conformers of such a residue instead. This has been observed twice, with opposite charge signs, both consistent with the same underlying mechanism:

- **4F5S (bovine serum albumin), residues 81 and 185**, two disordered arginines (positively charged, pKa 13.8). Verified exhaustively: an atom-by-atom diff against a live `prodes` run on the same chain found this was the *only* divergence across all 4,672 atoms. `prodes-rs`'s computed pI comes out *higher* than a live `prodes` re-run on this structure, consistent with `prodes` missing positive charge.
- **6PO0 (catalase), chain A residue 376**, a disordered cysteine (negatively charged, pKa 8.33; confirmed via 0.50/0.50 split occupancy between its two conformers). Not verified to the same exhaustive standard as 4F5S — this is based on a targeted search for disordered residues rather than a full atom-by-atom diff — but directionally consistent: `prodes-rs`'s computed pI (6.42) comes out *lower* than the paper's own published theoretical pI (6.80), consistent with the published value missing a negative charge contribution.

Net effect: the direction of the resulting pI gap depends on the charge sign of whichever residue is affected — positive residues missing from `prodes`'s count push its pI too high, negative residues push it too low. `SurfEpNegSumAverage` is affected the same way, since it depends on the same per-atom charge assignment. In every case observed so far, the divergence traces back to this one limitation in the reference tool's alternate-conformer handling, not an error in `prodes-rs`.

## Usage

```
cargo build --release
./target/release/prodes-rs <path-to-pdb> <path-to-output-csv> [--ph <comma-separated pH values>]
```

`--ph` defaults to `7`. Example, computing features across the four pH values used in the retention model's training data:

```
prodes-rs my_protein.pdb features.csv --ph 7,8,9,10
```

Each run appends one row per pH value to the output CSV (writing the header first if the file doesn't already exist), so repeated invocations across multiple structures accumulate into a single growing feature table:

```
id,ph,pi,surf_ep_neg_sum_average
my_protein,7,6.02,-12.4
my_protein,8,6.02,-58.1
...
```

The protein's chain A is used by default (falling back to the first chain present if none is named "A"). The isoelectric point is pH-independent and repeated on every row for convenience; the surface electrostatics feature is recomputed per pH value.

No example PDB is bundled with this repository — `prodes-rs` is a general-purpose tool, not tied to any particular research dataset. Structures can be sourced from the [RCSB PDB](https://www.rcsb.org/) or the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).

### Known input quirks

Two issues have come up when running real PDB/AlphaFold files through this tool, neither a bug in `prodes-rs` itself:

- **`SEQADV` records documenting a residue deletion relative to a reference sequence can fail to parse.** A `DELETION`-type `SEQADV` line leaves several fields legitimately blank (there's no PDB residue to report a number for), and `pdbtbx` expects a valid integer there — the parse fails with `InvalidatingError`, unconditionally, regardless of strictness settings. Confirmed on PDB entry 4PEP (porcine pepsin), which has one. Fix: strip the offending line(s) from the file before parsing — `SEQADV` is pure sequence-database metadata, never referenced by any computation here, so removing it has no effect on results:
  ```powershell
  (Get-Content file.pdb) | Where-Object { $_ -notmatch "^SEQADV" } | Set-Content file.pdb
  ```
- **AlphaFold models a UniProt entry's full sequence, including signal peptides and propeptides that aren't part of the mature protein.** If the source UniProt entry is a precursor (check its feature table for `SIGNAL`/`PROPEP` annotations), the downloaded structure will include residues that were never part of what was actually purified or measured — inflating the residue count and, since signal peptides tend to carry a positive charge (real biological signal: the "positive-inside rule"), skewing the isoelectric point. Confirmed on trypsin inhibitor A (UniProt P01070, AlphaFold model AF-P01070-F1): the full 216-residue precursor gives a computed pI of 4.76, while trimming to just the mature chain (UniProt-annotated residues 25–205) gives 4.41 — matching Neijenhuis et al.'s published 4.4 almost exactly. Fix: check the source UniProt entry's feature table before using an AlphaFold structure, and trim to the mature chain's residue range if the entry is a precursor.

## Project layout

```
src/
├── main.rs                          CLI entry point
├── io.rs                            PDB loading
├── core.rs                          pKa table, charge assignment, isoelectric point
└── calculations/
    ├── geometry.rs                  Point3D, distance, spatial grid cells, surface-point coarsening
    ├── standard_equations.rs        Henderson-Hasselbalch charge equations
    ├── sasa.rs                      Shrake-Rupley solvent-accessible surface area
    ├── distance_functions.rs        Coulomb potential
    └── electrostatics.rs            Surface electrostatic potential, negative-sum feature
```

## Testing

```
cargo test
```

Tests are self-contained where possible — synthetic structures are constructed inline via `pdbtbx`'s raw-string parsing rather than depending on bundled fixture files, in keeping with this repository not carrying research-specific data.

## License

MIT — see `LICENSE`.
