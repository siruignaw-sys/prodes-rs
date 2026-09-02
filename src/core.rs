pub fn pka_data(residue_name: &str) -> Option<(f64, bool, Vec<&'static str>)> {
    match residue_name {
        "ARG" => Some((13.8, true, vec!["NE", "NH1", "NH2"])),
        "ASP" => Some((3.86, false, vec!["OD1", "OD2"])),
        "CYS" => Some((8.33, false, vec!["SG"])),
        "GLU" => Some((4.25, false, vec!["OE1", "OE2"])),
        "HIS" => Some((6.0, true, vec!["ND1", "NE2"])),
        "LYS" => Some((10.5, true, vec!["NZ"])),
        "TYR" => Some((10.0, false, vec!["OH"])),
        _ => None
    }
}

pub fn termini_check(atom_name: &str, first: bool, last: bool) -> Option<(f64, bool)> {
    if atom_name == "N" && first {
        Some((9.69, true))
    }
    else if atom_name == "C" && last {
        Some((2.34, false))
    }
    else {
        None
    }
}
