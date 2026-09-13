use std::f64::consts::PI;

pub fn atom_charge_coulomb(charge: f64) -> f64{
    charge * 1.6e-19
}

pub fn charge_simple(charge: f64, distance: f64, dielectric_constant: f64) -> f64 {
    let absolute_permittivity = 8.854e-12;
    let permittivity = dielectric_constant * absolute_permittivity;

    charge / (permittivity * distance * 4.0 * PI)
}


// not yet used, will be used in shell electrostatics
#[allow(dead_code)]
pub fn potential_multiple_media(charge: f64, distance_dielectrics: Vec<(f64, f64)>) -> f64 {
    let absolute_permittivity = 8.854e-12;
    let mut denominator = 0.0;

    for (distance, dielectric_constant) in distance_dielectrics.iter() {
        let permittivity = dielectric_constant * absolute_permittivity;
        denominator += permittivity * distance;
    }
    charge / (denominator * 4.0 * PI)
    
}
