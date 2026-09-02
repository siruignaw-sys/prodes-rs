use pdbtbx::*;
use std::collections::HashMap;
mod calculations;
mod core;
mod io;
use crate::calculations::*;
use crate::io::*;
use crate::core::*;
use crate::calculations::sasa::{exposure, radius};

fn main() {
    let mut pdb = get_pdb("data/4F5S.pdb").expect("Failed to load PDB");
    let chain_index = pdb.chains().position(|c| c.id() == "A").unwrap_or(0);
    let chain = pdb.chains_mut().nth(chain_index).expect("PDB has no chains");

    chain.sort(); 

    ionizable_atoms(&chain);
    
    let atoms: Vec<&Atom> = chain.atoms().collect();
    
    let mut max_radius = 0.0;
    for atom in atoms.iter() {
        if let Some(rad) = radius(atom.element()){
            if rad + 1.4 > max_radius {
                max_radius = rad + 1.4;
            }
        }
    }
    let cell_size = max_radius * 2.0;
    let (exposed_points, fraction_exposed) = exposure(atoms, cell_size);
    println!("{}", exposed_points.len());

    
}










