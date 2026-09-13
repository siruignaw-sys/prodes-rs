use pdbtbx::*;

pub fn get_pdb(data: &str) -> Result<PDB, Vec<PDBError>> {
    let mut options = ReadOptions::new();
    options.set_level(StrictnessLevel::Loose);
    let (pdb, _warnings) = options.read(data)?;
    Ok(pdb)
}
