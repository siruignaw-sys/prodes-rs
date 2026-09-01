use pdbtbx::*;
use std::f64::consts::PI;


#[derive(Debug)]
struct Point3D {
    x: f64,
    y: f64,
    z: f64
}


fn main() {
    let pdb = get_pdb("data/4F5S.pdb");
    let atoms: Vec<&Atom> = pdb.atoms().collect();
    let mut exposed_points: Vec<Point3D> = Vec::new();
    let mut max_radius = 0;
    for atom in atoms {
        let rad = radius(atom.element().expect("NO ELEMENT")) + 1.4;
        if rad > max_radius {
            max_radius = rad;
        }
    }
    let cell_size = max_radius * 2;

    for (i, atom) in atoms.iter().enumerate() {
        let rad = match radius(atom.element()) {
            Some(r) => r,
            None => continue,
        };

        let (ax, ay, az) = atom.pos();
        let center = Point3D{ x: ax, y: ay, z: az };
        let own_radius = rad + 1.4;
        let points = fibonacci_sphere(center, own_radius);

        let mut exposed = 0;
        for point in &points {
            let mut blocked = false;
            for (j, neighbor) in atoms.iter().enumerate() {
                if i == j {
                    continue;
                }
                if let Some(n_rad) = radius(neighbor.element()) {
                    let (nx, ny, nz) = neighbor.pos();
                    let n_center = Point3D{ x: nx, y: ny, z: nz };
                    let neighbor_radius = n_rad + 1.4;

                    if distance(point, &n_center) < neighbor_radius {
                        blocked = true;
                        break;
                    }
                }
            }

            if !blocked {
                exposed += 1;
                exposed_points.push(point); 
            }
        }
        let fraction_exposed = exposed as f64 / points.len() as f64;
        println!("atom {} ({:?}): {:.2} exposed", i, atom.element(), fraction_exposed);
    }
}

fn distance(p1: &Point3D, p2: &Point3D) -> f64 {
    return ((p1.x-p2.x).powi(2) + (p1.y-p2.y).powi(2) + (p1.z-p2.z).powi(2)).sqrt()
}

fn get_pdb(data: &str) -> PDB {
    let mut options = ReadOptions::new();
    options.set_level(StrictnessLevel::Loose);

    let (pdb, _error) = options.read(data).unwrap();
    return pdb; 
}

fn radius(atom: Option<&Element>) -> Option<f64> {
     if let Some(a) = atom {
         return a.atomic_radius().van_der_waals;
     }
     return None;
}

fn fibonacci_sphere(center: Point3D, radius: f64) -> Vec<Point3D>{
    let mut points = Vec::new(); 
    let golden_angle = PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..=99 {
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


