use std::collections::HashMap;

use ptable::Element;

pub fn molecular_weight(atoms: HashMap<&str, u32>) -> Result<f32, &'static str> {
    let mut weight: f32 = 0 as f32;
    let mut err: bool = false;
    for (symbol, count) in atoms {
        match Element::from_symbol(symbol) {
            Some(el) => weight += el.get_atomic_mass() * count as f32,
            None => err = true,
        }
    }
    if !err {
        Ok(weight)
    } else {
        Err("Invalid molecular formula")
    }
}

#[cfg(test)]
mod tests {
    use crate::molecule::molecular_weight;
    use math::round::half_up;
    use std::collections::HashMap;

    fn round(weight: f32) -> f64 {
        half_up(weight as f64, 2)
    }

    #[test]
    fn ethane() {
        let molecule: HashMap<&str, u32> = [("C", 2), ("H", 6)].iter().cloned().collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(weight, 30.07);
    }

    #[test]
    fn cellulose() {
        let molecule: HashMap<&str, u32> =
            [("C", 6), ("H", 10), ("O", 5)].iter().cloned().collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(round(weight), 162.14);
    }

    #[test]
    fn vanadium_acetylacetonate() {
        let molecule: HashMap<&str, u32> = [("V", 1), ("C", 15), ("H", 21), ("O", 6)]
            .iter()
            .cloned()
            .collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(round(weight), 348.27);
    }
}
