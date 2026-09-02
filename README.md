# prodes-rs

A from-scratch Rust reimplementation of [`prodes`](https://github.com/tneijenhuis/prodes), a tool that calculates features from 3D protein structures with a focus on the protein surface — solvent-accessible surface area (SASA), ionizable-atom detection, and (in progress) surface electrostatics and isoelectric point.

This project mirrors `prodes`'s Python module layout module-for-module, reimplemented in Rust.

## Relationship to `prodes`

`prodes-rs` is not affiliated with the original authors; it's an independent reimplementation, built and verified against the original tool's output as a ground-truth reference. If you use this tool, please also cite the original:

> Neijenhuis, T., Le Bussy, O., Geldhof, G., Klijn, M. E., & Ottens, M. (2024). Predicting protein retention in ion-exchange chromatography using an open source QSPR workflow. *Biotechnology Journal*, 19, e2300708. https://doi.org/10.1002/biot.202300708

See also [`prodes-fork`](https://github.com/datacatalysis/prodes) (Mark Teese / 22DataCatalysis GmbH), an independent ~170x-faster fork of the original Python tool.

## Status

**Implemented and tested:**
- PDB parsing (via [`pdbtbx`](https://crates.io/crates/pdbtbx))
- Shrake-Rupley solvent-accessible surface area (SASA), with a spatial-grid neighbor search for performance
- Ionizable-atom detection: the seven ionizable amino acid side chains (Arg, Asp, Cys, Glu, His, Lys, Tyr) plus the N- and C-termini, each with pKa and charge sign
- Henderson-Hasselbalch fractional charge calculation at an arbitrary pH
- Correct handling of alternate-location (altLoc) conformers in crystal structures, so charge and surface calculations aren't double-counted on residues with multiple resolved conformations

**Not yet implemented:**
- Isoelectric point (pI) calculation
- Coulomb-potential-based surface electrostatics ("sum of negative surface electrostatics")

## Usage

```sh
cargo run --release -- path/to/structure.pdb
```

Analysis is restricted to chain A of the given structure. This is a deliberate choice: many PDB files contain multiple chains that are crystallographic copies of a single biological monomer rather than a true multi-chain assembly, and the PDB `MASTER` record doesn't reliably distinguish the two cases.

Run the test suite with:

```sh
cargo test
```

## Project layout

```
src/
├── main.rs                              orchestration: load PDB, select/sort chain A, run analyses
├── io.rs                                PDB loading
├── core.rs                              pKa table, ionizable-atom detection, charge assignment
└── calculations/
    ├── geometry.rs                      3D point/vector math, Fibonacci-sphere sampling
    ├── sasa.rs                          Shrake-Rupley SASA calculation
    ├── standard_equations.rs            Henderson-Hasselbalch charge equations
    └── distance_functions.rs            Coulomb potential (in progress)
```

## License

MIT — see [LICENSE](LICENSE).
