use std::f64::consts::PI;

#[derive(Debug, Clone, Copy)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

pub fn distance(p1: &Point3D, p2: &Point3D) -> f64 {
    return ((p1.x-p2.x).powi(2) + (p1.y-p2.y).powi(2) + (p1.z-p2.z).powi(2)).sqrt()
}

pub fn point_to_cell(point: Point3D, cell_size: f64) -> (i64, i64, i64) {
    let x_c = (point.x / cell_size).floor() as i64;
    let y_c = (point.y / cell_size).floor() as i64;
    let z_c = (point.z / cell_size).floor() as i64;
    return (x_c, y_c, z_c);
}

pub fn fibonacci_sphere(center: Point3D, radius: f64) -> Vec<Point3D>{
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
