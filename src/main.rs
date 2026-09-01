use pdbtbx::*;
use std::f64::consts::PI;


#[derive(Debug)]
struct Point3D {
    x: f64,
    y: f64,
    z: f64
}


fn main() {
    println!("{}", fibonacci_sphere(Point3D{x: 0.0, y: 0.0, z: 0.0}, 20.0).len());
}

fn get_stuff() {
    let mut options = ReadOptions::new();
    options.set_level(StrictnessLevel::Loose);

    let (pdb, _error) = options.read("data/4F5S.pdb").unwrap();
    
    if let Some(chain) = pdb.chains().find(|c| c.id() == "A") {
 
        let total_atoms = chain.atoms().count();
        println!("Total atoms: {}", total_atoms);

        let total_residues = chain.residues().count();
        println!("Total residues: {}", total_residues);


    }
    else {
        println!("Chain A not found");
    }
}

fn radius(atom: Option<&Element>) -> Option<f64> {
     if let Some(a) = atom {
         return a.atomic_radius().van_der_waals;
     }
     return None;
}

fn fibonacci_sphere(center: Point3D, radius: f64) -> Vec<Point3D>{
    let mut points = Vec::new(); 
    let golden_angle = PI * (3.0 - 5.0_f64.sprt());
    for i in 0..99 {
        let i_f = i as f64;
        
        let y = 1.0 - (i_f / 98.0) * 2.0;
        let ring_radius = (1.0 - y * y).max(0.0).sqrt();

        let theta = i_f * golden_angle;

        points.push(
            Point3D {
                x: radius * ring_radius * theta.cos() + center.x, 
                y: y * radius + center.y, 
                z: radius * ring_radius*theta.sin() + center.z
            }
        );
    }
    points
}
