use crate::calculations::geometry::*;
use std::f64::consts::PI;
use std::collections::HashMap;
use std::collections::HashSet;
use crate::calculations::distance_functions::*;
use pdbtbx::*;
use crate::core::*;

#[derive(Debug)]
pub struct Stats {
    pub mean: f64,
    pub trimean: f64,
    pub median: f64,
    pub sum: f64,
    pub std: f64,
}

pub struct ShellFeatures {
    pub shell_ep_max: f64,
    pub shell_ep_min: f64,
    pub shell_ep_stats: Stats,
    pub n_shell_pos_ep: usize,
    pub shell_pos_stats: Stats,
    pub n_shell_neg_ep: usize,
    pub shell_neg_stats: Stats,
}

pub fn calculate_shell_features(
    chain: &Chain,
    centroid: Point3D,
    surface_points: &[Point3D],
    ph: f64,
    n_points: usize,
) -> ShellFeatures {
    let surface_grid = build_surface_grid(surface_points, 2.0);
    let charged = charges_at_ph(chain, ph, ChargeMode::Formal);
    let direction_points = sunflower_sphere(centroid, 1.0, n_points);
    let mut direction_eps = Vec::new();
    for direction in &direction_points {
        let req_dis = required_distance(*direction, centroid, surface_points);
        let new_point = move_point(*direction, centroid, req_dis);
        let (a, b, c, d) = find_plane(new_point, centroid);
        let mut direction_ep = 0.0;
        for (atom, charge) in &charged {
            let (x, y, z) = atom.pos();
            let pos = Point3D {x: x, y: y, z: z};
            let projected = project_point(a, b, c, d, pos);  
            if let Some(exit) = find_exit(pos, projected, &surface_grid, 2.0) {
                direction_ep += map_ep_to_plane(pos, *charge, projected, exit);
            }
            
        }
        direction_eps.push(direction_ep);
    }
    let shell_ep_min = direction_eps.iter()
        .cloned()
        .fold(f64::MAX, f64::min);
    
    let shell_ep_max = direction_eps.iter()
        .cloned()
        .fold(f64::MIN, f64::max);
    let n_shell_pos_ep: Vec<f64> = direction_eps.iter().cloned().filter(|x| *x > 0.0).collect();
    let n_shell_neg_ep: Vec<f64> = direction_eps.iter().cloned().filter(|x| *x < 0.0).collect();
     return ShellFeatures {
        shell_ep_max: shell_ep_max,
        shell_ep_min: shell_ep_min,
        shell_ep_stats: standard_features(&direction_eps), 
        n_shell_pos_ep: n_shell_pos_ep.len(), 
        shell_pos_stats: standard_features(&n_shell_pos_ep),
        n_shell_neg_ep: n_shell_neg_ep.len(),
        shell_neg_stats: standard_features(&n_shell_neg_ep),
    };
}

pub fn standard_features(values:&[f64]) -> Stats {
    if values.is_empty() {
        return Stats { mean: 0.0, trimean: 0.0, median: 0.0, sum: 0.0, std: 0.0};
    }
    let sum = values.iter().sum();
    let mean = sum / values.len() as f64;
    let median = median(values).unwrap();
    let trimean = trimean(values);
    let std = (values.iter().map(|x| ((x - mean) as f64).powi(2)).sum::<f64>() / values.len() as f64).sqrt();
    return Stats {mean: mean, trimean: trimean, median: median, sum: sum, std: std};
}

pub fn find_plane(point_on_plane: Point3D, normal_vector_point: Point3D) -> (f64, f64, f64, f64) {
    let p = normal_vector_point - point_on_plane;
    let (a, b, c) = (p.x, p.y, p.z);
    let d = dot(Point3D{x: a, y: b, z: c}, point_on_plane);
    (a, b, c, d)
}

pub fn project_point(a: f64, b: f64, c: f64, d: f64, point: Point3D) -> Point3D {
    let k  = (d - a * point.x - b * point.y - c * point.z)
        / (a.powi(2) + b.powi(2) + c.powi(2));
    Point3D {
        x: point.x + a * k,
        y: point.y + b * k,
        z: point.z + c * k
    }
}

pub fn sunflower_sphere(center: Point3D, radius: f64, n_points: usize) -> Vec<Point3D> {
    let mut res: Vec<Point3D> = Vec::new();
    for i in 0..n_points {
        let index = i as f64 + 0.5;
        let phi = (1.0 - 2.0 * index / n_points as f64).acos();
        let theta = PI * (1.0 + 5.0_f64.sqrt()) * index;
        let x = theta.cos() * phi.sin() * radius + center.x;
        let y = theta.sin() * phi.sin() * radius + center.y;
        let z = phi.cos() * radius + center.z;
        res.push(Point3D { x: x, y: y, z: z}); 
    }
    res
}

pub fn move_point(point: Point3D, origin: Point3D, magnitude: f64) -> Point3D {
    let unit_vector = normalize(point - origin);
    origin + unit_vector * magnitude
}

pub fn maximal_distance(normal_vector: Point3D, vector_on_plane: Point3D, points: &[Point3D]) -> f64 {
    let unit_vector = normalize(normal_vector);
    let mut max: f64 = 0.0;
    for point in points {
        let d = dot(unit_vector, *point - vector_on_plane);
        if d > max {
            max = d;
        }    
    }
    max
}

pub fn required_distance(point_for_plane: Point3D, structure_center: Point3D, surface_points: &[Point3D]) -> f64 {
    let normal_vector = point_for_plane - structure_center;
    maximal_distance(normal_vector, structure_center, surface_points) + 1.0
}

pub fn build_surface_grid(points: &[Point3D], cell_size: f64) -> HashMap<(i64, i64, i64), Vec<Point3D>> {
    let mut map: HashMap<(i64, i64, i64), Vec<Point3D>> = HashMap::new();
    for point in points {
        let (xc, yc, zc) = point_to_cell(point.clone(), cell_size);
        map.entry((xc, yc, zc)).or_default().push(*point);
    }
    map
}

pub fn find_exit(
    point_vector: Point3D,
    projected_point_vector: Point3D,
    grid: &HashMap<(i64, i64, i64), Vec<Point3D>>,
    cell_size: f64,
) -> Option<Point3D> {
    let direction = normalize(projected_point_vector - point_vector);
    let total_distance = magnitude (projected_point_vector - point_vector);

    let mut visited: HashSet<(i64, i64, i64)> = HashSet::new();

    let steps = (total_distance.ceil() as usize) * 2;

    for i in 0..steps{
        let sample_point = point_vector + direction * (i as f64 / 2.0);
        let sample_cell = point_to_cell(sample_point, cell_size);
        let neighbors = neighbor_cells(sample_cell);
        visited.insert(sample_cell);
        for cell in neighbors {
            visited.insert(cell);
        }
    }

    let mut points: Vec<Point3D> = Vec::new();
    for cell in visited {
        if let Some(pts) = grid.get(&cell) {
            points.extend(pts);
        }
    }
    let mut max: f64 = 0.0;
    let mut best_exit: Option<Point3D> = None;
    for point in points {
        let vector = point - point_vector;
        let dot_prod = dot(direction, vector);
        if dot_prod <= 0.0 {
            continue;
        }
        let potential_exit = point_vector + direction * dot_prod;
        let perp_distance = magnitude(potential_exit - point);
        if (perp_distance * 10.0).round() / 10.0 <= 1.0 {
            if dot_prod > max {
                max = dot_prod;
                best_exit = Some(potential_exit);
            }
        }
    }
    best_exit
}

pub fn map_ep_to_plane(atom_position: Point3D, charge: f64, projected_point: Point3D, surface_exit: Point3D) -> f64 {
    let total_distance = magnitude(projected_point - atom_position);
    let protein_distance = magnitude(surface_exit - atom_position);
    let solvent_segment = (total_distance - protein_distance) * 1e-10;
    let protein_segment = protein_distance * 1e-10;
    let atom_charge = atom_charge_coulomb(charge);
    potential_multiple_media(atom_charge, vec![(solvent_segment, 80.0), (protein_segment, 4.0)])
}

pub fn median(values: &[f64]) -> Option<f64> {
    let mut valid_numbers: Vec<f64> = values.iter()
        .copied()
        .filter(|x| !x.is_nan())
        .collect();
    if valid_numbers.is_empty() {
        return None;
    }

    valid_numbers.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let mid = valid_numbers.len() / 2;
    if valid_numbers.len() %2 == 0 {
        Some((valid_numbers[mid-1] + valid_numbers[mid])/2.0)
    } else {
        Some(valid_numbers[mid])
    }
}

pub fn trimean(values: &[f64]) -> f64 {
    let q2 = median(values).unwrap();
    let lesser: Vec<f64> = values.iter().filter(|&&x| x < q2).copied().collect();
    let greater: Vec<f64> = values.iter().filter(|&&x| x > q2).copied().collect();
    let q1 = median(&lesser).unwrap_or(f64::NAN);
    let q3 = median(&greater).unwrap_or(f64::NAN);
    (q2*2.0 + q1 + q3)/4.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sunflower_sphere_points_sit_at_the_given_radius() {
        let center = Point3D { x: 1.0, y: -2.0, z: 3.0 };
        let radius = 7.5;
        let n_points = 120;

        let points = sunflower_sphere(center, radius, n_points);
        assert_eq!(points.len(), n_points, "should generate exactly n_points");

        for point in points {
            let dist = magnitude(point - center);
            assert!(
                (dist - radius).abs() < 1e-9,
                "point should sit at exactly `radius` from center, got distance {}",
                dist
            );
        }
    }
}
