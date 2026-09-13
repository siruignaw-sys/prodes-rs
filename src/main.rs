use pdbtbx::*;
mod calculations;
mod core;
mod io;
use crate::calculations::geometry::*;
use crate::io::*;
use crate::core::*;
use crate::calculations::sasa::*;
use crate::calculations::electrostatics::*;
use clap::Parser;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;


#[derive(Parser)]
struct Args {
    pdb_path: String,

    out_file: String,

    #[arg(long, default_value = "7")]
    ph: String,
}

fn main() {
    let args = Args::parse();
    let ph_values: Vec<f64> = args.ph.split(',').map(|s| s.trim().parse().unwrap()).collect();
    let path = args.pdb_path;
    
    let out_file = args.out_file;
    
    let file_exists = Path::new(&out_file).exists();

    let header = "id,ph,pi,surf_ep_neg_sum_average\n";

    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&out_file)
        .expect("failed to open output file");
    
    if !file_exists {
        file.write_all(header.as_bytes()).expect("failed to write header");
    }

    let mut pdb = get_pdb(&path).expect("Failed to load pdb");
    let chain_index = pdb.chains().position(|c| c.id() == "A").unwrap_or(0);
    let chain = pdb.chains_mut().nth(chain_index).expect("PDB has no chains");

    chain.sort(); 

    let atoms: Vec<&Atom> = chain.atoms().collect();
    let iso_point = isoelectric_point(&chain);
    let mut max_radius = 0.0;
    for atom in atoms.iter() {
        if let Some(rad) = radius(atom.element()){
            if rad + 1.4 > max_radius {
                max_radius = rad + 1.4;
            }
        }
    }
    let cell_size = max_radius * 2.0;
    
    let (exposed_points, _fraction_exposed) = exposure(atoms, cell_size);
    let property_points = property_points_on_surface(&exposed_points, 1.0);
  
    let id = Path::new(&path).file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
    
    for ph in ph_values {
        let charged_atoms = charges_at_ph(&chain, ph);
        
        let potentials: Vec<f64> = property_points.iter()
            .map(|p| electrostatic_potential(p, &charged_atoms))
            .collect();

        let surf_ep_neg_sum_val = surf_ep_neg_sum(&property_points, &charged_atoms);
        let row = format!("{},{},{},{}\n", id, ph, iso_point, surf_ep_neg_sum_val);
        file.write_all(row.as_bytes()).expect("failed to write row");
    }
}










