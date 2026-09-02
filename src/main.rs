use pdbtbx::*;
use std::collections::HashMap;
mod calculations;
mod core;
mod io;
use crate::calculations::*;
use crate::io::*;
use crate::calculations::sasa::{exposure, radius};

fn main() {
    let pdb = get_pdb("data/4F5S.pdb");
    let atoms: Vec<&Atom> = if let Some(chain) = pdb.chains().find(|c| c.id() == "A") {
        chain.atoms().collect()
    }
    else {
        pdb.atoms().collect()
    };
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










