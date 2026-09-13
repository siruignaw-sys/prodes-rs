use std::f64::consts::PI;
use std::collections::HashSet;

#[derive(Debug, Clone, Copy)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

pub fn distance(p1: &Point3D, p2: &Point3D) -> f64 {
    ((p1.x-p2.x).powi(2) + (p1.y-p2.y).powi(2) + (p1.z-p2.z).powi(2)).sqrt()
}

pub fn point_to_cell(point: Point3D, cell_size: f64) -> (i64, i64, i64) {
    let x_c = (point.x / cell_size).floor() as i64;
    let y_c = (point.y / cell_size).floor() as i64;
    let z_c = (point.z / cell_size).floor() as i64;
    (x_c, y_c, z_c)
}

pub fn property_points_on_surface(surface_points: &[Point3D], cell_size: f64) -> Vec<Point3D> {
    let mut set: HashSet<(i64, i64, i64)> = HashSet::new();
    let mut out = Vec::new();
    for point in surface_points.iter(){
        let (x_g, y_g, z_g) = point_to_cell(*point, cell_size);
        set.insert((x_g, y_g, z_g));
    }
    for center in set.iter() {
        let center_point = Point3D { 
            x: center.0 as f64 * cell_size + cell_size / 2.0, 
            y: center.1 as f64 * cell_size + cell_size / 2.0, 
            z: center.2 as f64 * cell_size + cell_size / 2.0 
        };
        out.push(center_point);
    }
    out
}

pub fn fibonacci_sphere(center: Point3D, radius: f64, n_points: usize) -> Vec<Point3D>{
    let mut points = Vec::new(); 
    let golden_angle = PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..n_points {
        let i_f = i as f64;
        
        let y = 1.0 - (i_f / (n_points as f64 - 1.0)) * 2.0;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_count_matches_n_points_exactly() {
        let center = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let points = fibonacci_sphere(center, 5.0, 250);
        assert_eq!(points.len(), 250, "should generate exactly n_points, not a fixed constant");
    }

    #[test]
    fn first_and_last_points_sit_at_the_poles() {
        let center = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let radius = 5.0;
        let n_points = 250;
        let points = fibonacci_sphere(center, radius, n_points);

        let first = &points[0];
        let last = &points[n_points - 1];

        // i=0 should sit at y = +radius (north pole), i=n_points-1 at y = -radius (south pole).
        // This is the exact invariant the old hardcoded `98.0` denominator violated: with a
        // variable n_points, only `n_points - 1` keeps both poles on the sphere.
        assert!((first.y - radius).abs() < 1e-9, "first point should be at the north pole, got y={}", first.y);
        assert!((last.y - (-radius)).abs() < 1e-9, "last point should be at the south pole, got y={}", last.y);
    }

    #[test]
    fn poles_hold_for_a_different_n_points_too() {
        // Same check as above but with a different point count, to confirm the fix is
        // relative to n_points and not just correct for whatever value was tested above.
        let center = Point3D { x: 1.0, y: 2.0, z: 3.0 };
        let radius = 2.5;
        let n_points = 50;
        let points = fibonacci_sphere(center, radius, n_points);

        let first = &points[0];
        let last = &points[n_points - 1];

        assert!((first.y - (radius + center.y)).abs() < 1e-9, "first point should be at the north pole");
        assert!((last.y - (-radius + center.y)).abs() < 1e-9, "last point should be at the south pole");
    }
}
