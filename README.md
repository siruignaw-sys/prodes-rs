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

- **Isoelectric point** (`core.rs`) — Henderson-Hasselbalch charge model over a standard ionizable-residue pKa table, including N/C-termini, alternate-location conformer deduplication, and correct handling of trailing HETATM records (e.g. crystallographic waters) that would otherwise be mistaken for a chain's true C-terminus. Multi-atom ionizable groups (e.g. arginine's guanidinium, aspartate's carboxylate) correctly split their group charge across their constituent atoms rather than over-counting it.
- **Solvent-accessible surface area** (`calculations/sasa.rs`) — Shrake-Rupley algorithm with a spatial-hash grid for neighbor search. Sample point density per atom scales with each atom's own probe-inflated surface area, matching `prodes`'s approach rather than using a fixed point count; cross-checked against the real `prodes` output on a real PDB structure (1GDW), landing within ~1.5% of the reference point count.
- **Surface electrostatic potential** (`calculations/electrostatics.rs`, `calculations/distance_functions.rs`) — Coulomb potential summed at a coarsened grid of surface points (`calculations/geometry.rs::property_points_on_surface`), filtered to negative-potential points and summed, matching `prodes`'s `SurfEpNegSumAverage` feature. The underlying Coulomb formulas are verified bit-for-bit against the Python reference on matched inputs.

Not yet implemented:

- The "Formal" (binarized ±1/n) charge mode — only the continuous Henderson-Hasselbalch ("Average") charge mode is implemented, which is what the two-feature retention model actually uses.
- `prodes`'s "Shell" electrostatics feature family (potential projected onto a plane with a two-dielectric protein/solvent boundary) — a separate, more elaborate feature set in the original tool that the two-feature model doesn't require.
- Hydrophobicity/lipophilicity features, surface shape descriptors, and custom pKa file import (PROPKA/H++/pypka) — all present in `prodes`, none currently ported.

This has been verified through unit tests on synthetic structures, hand-computed numeric checks against the underlying physics, and targeted comparisons against the real Python `prodes` output where practical — not through an exhaustive parity test suite across `prodes`'s full feature set.

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
