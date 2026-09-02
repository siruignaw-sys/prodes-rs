use pdbtbx::*;
use std::collections::HashMap;
use crate::calculations::geometry::{Point3D, point_to_cell, fibonacci_sphere, distance};
use std::f64::consts::PI;

pub fn exposure(atoms: Vec<&Atom>, cell_size: f64) -> (Vec<Point3D>, Vec<(f64, &Atom)>) {
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
    let mut fraction_exposed: Vec<(f64, &Atom)> = Vec::new();

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


        let n_points = (own_radius.powi(2) * 4.0 * PI * 2.0) as usize;

        let points = fibonacci_sphere(a_center, own_radius, n_points);
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
        fraction_exposed.push((exposed as f64 / total_points as f64, atom));
    }
    return (exposed_points, fraction_exposed);
}

pub fn radius(atom: Option<&Element>) -> Option<f64> {
     if let Some(a) = atom {
         return a.atomic_radius().van_der_waals;
     }
     return None;
}
