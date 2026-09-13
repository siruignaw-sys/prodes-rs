use pdbtbx::*;
use crate::calculations::standard_equations::{pos_charge, neg_charge, pos_charge_formal, neg_charge_formal};
use clap::ValueEnum;

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ChargeMode {
    Average,
    Formal
}

pub fn pka_data(residue_name: &str) -> Option<(f64, bool, Vec<&'static str>)> {
    match residue_name {
        "ARG" => Some((13.8, true, vec!["NE", "NH1", "NH2"])),
        "ASP" => Some((3.86, false, vec!["OD1", "OD2"])),
        "CYS" => Some((8.33, false, vec!["SG"])),
        "GLU" => Some((4.25, false, vec!["OE1", "OE2"])),
        "HIS" => Some((6.0, true, vec!["ND1", "NE2"])),
        "LYS" => Some((10.5, true, vec!["NZ"])),
        "TYR" => Some((10.0, false, vec!["OH"])),
        _ => None
    }
}

pub fn termini_check(atom_name: &str, first: bool, last: bool) -> Option<(f64, bool)> {
    if atom_name == "N" && first {
        Some((9.69, true))
    }
    else if atom_name == "C" && last {
        Some((2.34, false))
    }
    else {
        None
    }
}

pub fn ionizable_atoms(chain: &Chain) -> Vec<(&Atom, f64, bool, usize)> {
    let mut result = Vec::new();
    let polymer_residues: Vec<&Residue> = chain.residues().filter(|r| is_polymer(r)).collect();
    let residue_count = polymer_residues.len(); 

    for (i, residue) in polymer_residues.into_iter().enumerate() {
        let is_first = i == 0;
        let is_last = i == residue_count - 1;
    
        let Some(conformer) = residue.conformers().find(|c| c.alternative_location().is_none())
            .or_else(|| residue.conformers().next())
            else {
                continue;
            };
        if let Some((pka, positive, atom_names)) = residue.name().and_then(pka_data) {
            let matched: Vec<&Atom> = conformer.atoms()
                .filter(|a| atom_names.contains(&a.name()))
                .collect();
            let group_size = matched.len();
            for atom in matched {
                result.push((atom, pka, positive, group_size));
            }
        }

        for atom in conformer.atoms() {
            if let Some((pka, positive)) = termini_check(atom.name(), is_first, is_last) {
                result.push((atom, pka, positive, 1));
            }
        }
    }
    result
}

pub fn is_polymer(residue: &Residue) -> bool {
    residue.conformers().next()
        .and_then(|c| c.atoms().next())
        .map(|a| !a.hetero())
        .unwrap_or(false)
}

pub fn charges_at_ph(chain: &Chain, ph: f64, mode: ChargeMode) -> Vec<(&Atom, f64)> {
    ionizable_atoms(chain)
        .into_iter()
        .map(|(atom, pka, positive, group_size)| {
            match mode {
                ChargeMode::Formal => {
                    let charge = if positive { pos_charge_formal(pka, ph) } else { neg_charge_formal (pka, ph) };
                    (atom, charge / group_size as f64)
                },
                ChargeMode::Average => {
                    let charge = if positive { pos_charge(pka, ph) } else { neg_charge (pka, ph) };
                    (atom, charge / group_size as f64)

                }
            }
        })
    .collect()
}

pub fn total_charge_at_ph(chain: &Chain, ph: f64, mode: ChargeMode) -> f64 {
    charges_at_ph(chain, ph, mode).iter().map(|(_, q)| q).sum()
}

pub fn isoelectric_point(chain: &Chain, mode: ChargeMode) -> f64 {
    let mut low = 0.0;
    let mut high = 14.0;

    while high - low > 0.0000001 {
        let mid = (low + high) / 2.0;
        let charge = total_charge_at_ph(chain, mid, mode);
        if charge > 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    (low + high) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    const TOY_CHAIN_PDB: &str = "ATOM      1  N   ASP A   1       0.000   0.000   1.000  1.00  0.00           N
    \nATOM      2  CA  ASP A   1       0.000   0.000   2.000  1.00  0.00           C
    \nATOM      3  C   ASP A   1       0.000   0.000   3.000  1.00  0.00           C
    \nATOM      4  O   ASP A   1       0.000   0.000   4.000  1.00  0.00           O
    \nATOM      5  CB  ASP A   1       0.000   0.000   5.000  1.00  0.00           C
    \nATOM      6  OD1 ASP A   1       0.000   0.000   6.000  1.00  0.00           O
    \nATOM      7  OD2 ASP A   1       0.000   0.000   7.000  1.00  0.00           O
    \nATOM      8  N  AARG A   2       0.000   0.000   8.000  0.50  0.00           N
    \nATOM      9  CA AARG A   2       0.000   0.000   9.000  0.50  0.00           C
    \nATOM     10  C  AARG A   2       0.000   0.000  10.000  0.50  0.00           C
    \nATOM     11  O  AARG A   2       0.000   0.000  11.000  0.50  0.00           O
    \nATOM     12  CB AARG A   2       0.000   0.000  12.000  0.50  0.00           C
    \nATOM     13  CG AARG A   2       0.000   0.000  13.000  0.50  0.00           C
    \nATOM     14  CD AARG A   2       0.000   0.000  14.000  0.50  0.00           C
    \nATOM     15  NE AARG A   2       0.000   0.000  15.000  0.50  0.00           N
    \nATOM     16  CZ AARG A   2       0.000   0.000  16.000  0.50  0.00           C
    \nATOM     17  NH1AARG A   2       0.000   0.000  17.000  0.50  0.00           N
    \nATOM     18  NH2AARG A   2       0.000   0.000  18.000  0.50  0.00           N
    \nATOM     19  N  BARG A   2       0.000   0.000  19.000  0.50  0.00           N
    \nATOM     20  CA BARG A   2       0.000   0.000  20.000  0.50  0.00           C
    \nATOM     21  C  BARG A   2       0.000   0.000  21.000  0.50  0.00           C
    \nATOM     22  O  BARG A   2       0.000   0.000  22.000  0.50  0.00           O
    \nATOM     23  CB BARG A   2       0.000   0.000  23.000  0.50  0.00           C
    \nATOM     24  CG BARG A   2       0.000   0.000  24.000  0.50  0.00           C
    \nATOM     25  CD BARG A   2       0.000   0.000  25.000  0.50  0.00           C
    \nATOM     26  NE BARG A   2       0.000   0.000  26.000  0.50  0.00           N
    \nATOM     27  CZ BARG A   2       0.000   0.000  27.000  0.50  0.00           C
    \nATOM     28  NH1BARG A   2       0.000   0.000  28.000  0.50  0.00           N
    \nATOM     29  NH2BARG A   2       0.000   0.000  29.000  0.50  0.00           N
    \nATOM     30  N   LYS A   3       0.000   0.000  30.000  1.00  0.00           N
    \nATOM     31  CA  LYS A   3       0.000   0.000  31.000  1.00  0.00           C
    \nATOM     32  C   LYS A   3       0.000   0.000  32.000  1.00  0.00           C
    \nATOM     33  O   LYS A   3       0.000   0.000  33.000  1.00  0.00           O
    \nATOM     34  CB  LYS A   3       0.000   0.000  34.000  1.00  0.00           C
    \nATOM     35  CG  LYS A   3       0.000   0.000  35.000  1.00  0.00           C
    \nATOM     36  CD  LYS A   3       0.000   0.000  36.000  1.00  0.00           C
    \nATOM     37  CE  LYS A   3       0.000   0.000  37.000  1.00  0.00           C
    \nATOM     38  NZ  LYS A   3       0.000   0.000  38.000  1.00  0.00           N
    \nTER\nEND\n";
    
    const TOY_CHAIN_WITH_WATER_PDB: &str = "ATOM      1  N   ASP A   1       0.000   0.000   1.000  1.00  0.00           N
    \nATOM      2  CA  ASP A   1       0.000   0.000   2.000  1.00  0.00           C
    \nATOM      3  C   ASP A   1       0.000   0.000   3.000  1.00  0.00           C
    \nATOM      4  O   ASP A   1       0.000   0.000   4.000  1.00  0.00           O
    \nATOM      5  CB  ASP A   1       0.000   0.000   5.000  1.00  0.00           C
    \nATOM      6  OD1 ASP A   1       0.000   0.000   6.000  1.00  0.00           O
    \nATOM      7  OD2 ASP A   1       0.000   0.000   7.000  1.00  0.00           O
    \nATOM      8  N  AARG A   2       0.000   0.000   8.000  0.50  0.00           N
    \nATOM      9  CA AARG A   2       0.000   0.000   9.000  0.50  0.00           C
    \nATOM     10  C  AARG A   2       0.000   0.000  10.000  0.50  0.00           C
    \nATOM     11  O  AARG A   2       0.000   0.000  11.000  0.50  0.00           O
    \nATOM     12  CB AARG A   2       0.000   0.000  12.000  0.50  0.00           C
    \nATOM     13  CG AARG A   2       0.000   0.000  13.000  0.50  0.00           C
    \nATOM     14  CD AARG A   2       0.000   0.000  14.000  0.50  0.00           C
    \nATOM     15  NE AARG A   2       0.000   0.000  15.000  0.50  0.00           N
    \nATOM     16  CZ AARG A   2       0.000   0.000  16.000  0.50  0.00           C
    \nATOM     17  NH1AARG A   2       0.000   0.000  17.000  0.50  0.00           N
    \nATOM     18  NH2AARG A   2       0.000   0.000  18.000  0.50  0.00           N
    \nATOM     19  N  BARG A   2       0.000   0.000  19.000  0.50  0.00           N
    \nATOM     20  CA BARG A   2       0.000   0.000  20.000  0.50  0.00           C
    \nATOM     21  C  BARG A   2       0.000   0.000  21.000  0.50  0.00           C
    \nATOM     22  O  BARG A   2       0.000   0.000  22.000  0.50  0.00           O
    \nATOM     23  CB BARG A   2       0.000   0.000  23.000  0.50  0.00           C
    \nATOM     24  CG BARG A   2       0.000   0.000  24.000  0.50  0.00           C
    \nATOM     25  CD BARG A   2       0.000   0.000  25.000  0.50  0.00           C
    \nATOM     26  NE BARG A   2       0.000   0.000  26.000  0.50  0.00           N
    \nATOM     27  CZ BARG A   2       0.000   0.000  27.000  0.50  0.00           C
    \nATOM     28  NH1BARG A   2       0.000   0.000  28.000  0.50  0.00           N
    \nATOM     29  NH2BARG A   2       0.000   0.000  29.000  0.50  0.00           N
    \nATOM     30  N   LYS A   3       0.000   0.000  30.000  1.00  0.00           N
    \nATOM     31  CA  LYS A   3       0.000   0.000  31.000  1.00  0.00           C
    \nATOM     32  C   LYS A   3       0.000   0.000  32.000  1.00  0.00           C
    \nATOM     33  O   LYS A   3       0.000   0.000  33.000  1.00  0.00           O
    \nATOM     34  CB  LYS A   3       0.000   0.000  34.000  1.00  0.00           C
    \nATOM     35  CG  LYS A   3       0.000   0.000  35.000  1.00  0.00           C
    \nATOM     36  CD  LYS A   3       0.000   0.000  36.000  1.00  0.00           C
    \nATOM     37  CE  LYS A   3       0.000   0.000  37.000  1.00  0.00           C
    \nATOM     38  NZ  LYS A   3       0.000   0.000  38.000  1.00  0.00           N
    \nHETATM   39  O   HOH A 101       0.000   0.000  39.000  1.00  0.00           O
    \nHETATM   40  O   HOH A 102       0.000   0.000  40.000  1.00  0.00           O
    \nTER\nEND\n";

    #[test]
    fn trailing_hetatm_waters_do_not_hide_the_true_c_terminus() {
        let mut options = ReadOptions::new();
        options.set_level(StrictnessLevel::Loose);
        options.set_format(Format::Pdb);
        let cursor = Cursor::new(TOY_CHAIN_WITH_WATER_PDB.as_bytes());
        let (pdb, _warnings) = options.read_raw(BufReader::new(cursor)).unwrap();
        let mut chain = pdb.chains().find(|c| c.id() == "A").unwrap().clone();
        chain.sort();

        // sanity check the fixture: two HOH residues really are in this chain, after LYS
        let water_count = chain.residues().filter(|r| r.name() == Some("HOH")).count();
        assert_eq!(water_count, 2, "test fixture should have 2 waters");

        let charged = ionizable_atoms(&chain);
        let c_term: Vec<_> = charged.iter().filter(|(_, pka, _, _)| *pka == 2.34).collect();
        assert_eq!(c_term.len(), 1, "LYS's C should still be flagged as C-terminus despite trailing waters");

        // should match the water-free baseline exactly: waters contribute nothing
        assert_eq!(charged.len(), 8);
    }

    fn parse_toy_chain() -> Chain {
        let mut options = ReadOptions::new();
        options.set_level(StrictnessLevel::Loose);
        options.set_format(Format::Pdb);
        let cursor = Cursor::new(TOY_CHAIN_PDB.as_bytes());
        let (pdb, _warnings) = options.read_raw(BufReader::new(cursor)).unwrap();
        let mut chain = pdb.chains().find(|c| c.id() == "A").unwrap().clone();
        chain.sort();
        chain
    }

    #[test]
    fn arg_with_two_full_conformers_is_not_double_counted() {
        let chain = parse_toy_chain();
        let arg = chain.residues().find(|r| r.serial_number() == 2).unwrap();
        assert_eq!(arg.conformer_count(), 2, "test fixture should have 2 ARG conformers");

        let charged = ionizable_atoms(&chain);
        let arg_atoms: Vec<_> = charged.iter().filter(|(_, pka, _, _)| *pka == 13.8).collect();
        assert_eq!(arg_atoms.len(), 3, "ARG side chain should contribute NE/NH1/NH2 once, not twice");
    }

    #[test]
    fn termini_only_hit_first_and_last_residue() {
        let chain = parse_toy_chain();
        let charged = ionizable_atoms(&chain);

        let n_term: Vec<_> = charged.iter().filter(|(_, pka, _, _)| *pka == 9.69).collect();
        let c_term: Vec<_> = charged.iter().filter(|(_, pka, _, _)| *pka == 2.34).collect();
        assert_eq!(n_term.len(), 1, "exactly one N-terminus hit expected");
        assert_eq!(c_term.len(), 1, "exactly one C-terminus hit expected");
    }

    #[test]
    fn total_ionizable_atom_count_is_correct() {
        let chain = parse_toy_chain();
        let charged = ionizable_atoms(&chain);
        // ASP side chain (2) + N-term (1) + ARG side chain deduped (3) + LYS side chain (1) + C-term (1)
        assert_eq!(charged.len(), 8);
    }
}
