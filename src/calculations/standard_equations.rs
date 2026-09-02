pub fn pos_charge(pka: f64, ph: f64) -> f64 {
    1.0 / (1.0 + 10f64.powf(ph - pka))
}

pub fn neg_charge(pka: f64, ph: f64) -> f64 {
    -1.0 / (1.0 + 10f64.powf(pka-ph))
}

