use crate::calculations::geometry::*;
use crate::calculations::distance_functions::*;
use pdbtbx::*;

pub fn electrostatic_potential(point: &Point3D, charged_atoms: &[(&Atom, f64)]) -> f64 {
    let mut sum = 0.0;
    for (atom, charge) in charged_atoms {
        let (ax, ay, az) = atom.pos();
        let position = Point3D {x: ax, y: ay, z: az};
        let d = distance(&position, point) * 1e-10;
        let charge_coulombs = atom_charge_coulomb(*charge);
        sum += charge_simple(charge_coulombs, d, 4.0);
    }
    sum
}

pub fn surf_ep_neg_sum(property_points: &[Point3D], charged_atoms: &[(&Atom, f64)]) -> f64 {
    let mut sum = 0.0;
    for point in property_points {
        let pot = electrostatic_potential(point, charged_atoms);
        sum = if pot < 0.0 { sum + pot } else { sum };
    }
    sum
}

#[cfg(test)]
mod checks {
    use super::*;
    use std::io::{BufReader, Cursor};

    const ONE_ATOM_PDB: &str = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nTER\nEND\n";

    #[test]
    fn matches_hand_computed_potential() {
        let mut options = ReadOptions::new();
        options.set_level(StrictnessLevel::Loose);
        options.set_format(Format::Pdb);
        let cursor = Cursor::new(ONE_ATOM_PDB.as_bytes());
        let (pdb, _warnings) = options.read_raw(BufReader::new(cursor)).unwrap();
        let atom = pdb.atoms().next().unwrap();

        // atom at origin, charge = 0.5, evaluated at a point 10.0 Angstrom away along z
        let charge = 0.5;
        let charged_atoms: Vec<(&Atom, f64)> = vec![(atom, charge)];
        let eval_point = Point3D { x: 0.0, y: 0.0, z: 10.0 };

        let got = electrostatic_potential(&eval_point, &charged_atoms);

        // hand-computed expected value using the exact same formulas, independently
        let charge_coulombs = charge * 1.6e-19;
        let dist_m = 10.0 * 1e-10;
        let absolute_permittivity = 8.854e-12;
        let permittivity = 4.0 * absolute_permittivity;
        let expected = charge_coulombs / (permittivity * dist_m * 4.0 * std::f64::consts::PI);

        println!("got = {}, expected = {}", got, expected);
        assert!((got - expected).abs() < 1e-6, "got {}, expected {}", got, expected);
    }
}
