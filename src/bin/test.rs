use pdbtbx::*;
use std::f64::consts::PI;
use std::collections::HashMap;


#[derive(Debug, Clone, Copy)]
struct Point3D {
    x: f64,
    y: f64,
    z: f64
}


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

    let mut map: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, atom) in atoms.iter().enumerate() {
        if radius(atom.element()).is_none() {
            continue;
        }
        let (ax, ay, az) = atom.pos();
        let center = Point3D{ x: ax, y: ay, z: az };
        let (xc, yc, zc) = point_to_cell(center, cell_size);
        map.entry((xc, yc, zc)).or_default().push(i);
    }

    let mut exposed_points = Vec::new();

    for (i, atom) in atoms.iter().enumerate() {
        let rad = match radius(atom.element()) {
            Some(r) => r,
            None => continue
        };
        let own_radius = rad + 1.4;
        let (ax, ay, az) = atom.pos();
        let a_center = Point3D{x: ax, y: ay, z: az};
        let cell = point_to_cell(a_center, cell_size);
        
        let mut candidates: Vec<usize> = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = map.get(&(cell.0 + dx, cell.1 + dy, cell.2 + dz)) {
                        candidates.extend(v.iter().copied());
                    }               
                }
            }
        }

        let points = fibonacci_sphere(a_center, own_radius);
        let total_points = points.len();
        let mut exposed = 0;
        
        for point in points {
            let mut blocked = false;
            for &idx in &candidates {
                if idx == i {
                    continue;
                }
                let neighbor = atoms[idx];
                if let Some(n_rad) = radius(neighbor.element()) {
                    let (nx, ny, nz) = neighbor.pos();
                    let n_center = Point3D{x: nx, y: ny, z: nz};
                    let neighbor_radius = n_rad + 1.4;
                    if distance(&point, &n_center) < neighbor_radius {
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
        let fraction_exposed = exposed as f64 / total_points as f64;
        println!("atom {} ({:?}): {:.2} exposed", i, atom.element(), fraction_exposed);
    }

    println!("total surface points collected: {}", exposed_points.len());
}

fn point_to_cell(point: Point3D, cell_size: f64) -> (i64, i64, i64) {
    let x_c = (point.x / cell_size).floor() as i64;
    let y_c = (point.y / cell_size).floor() as i64;
    let z_c = (point.z / cell_size).floor() as i64;
    return (x_c, y_c, z_c);
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


