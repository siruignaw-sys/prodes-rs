use pdbtbx::*;

pub fn get_pdb(data: &str) -> PDB {
    let mut options = ReadOptions::new();
    options.set_level(StrictnessLevel::Loose);
    let (pdb, _error) = options.read(data).unwrap();

    return pdb; 
}
