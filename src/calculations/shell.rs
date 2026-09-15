use crate::calculations::geometry::*;
use std::f64::consts::PI;
use std::collections::HashMap;

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
        map.entry((xc, yc, zc)).or_default().push(point.clone());
    }
    map
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
