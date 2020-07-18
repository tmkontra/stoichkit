use std::collections::HashMap;

use ptable::Element;


pub fn molecular_weight(atoms: HashMap<String, u32>) -> Result<f32, &'static str> {
    let mut weight: f32 = 0 as f32;
    let mut err: bool = false;
    for (symbol, count) in atoms {
        match Element::from_symbol(symbol.as_str()) {
            Some(el) => weight += el.get_atomic_mass() * count as f32,
            None => err = true
        }
    }
    if !err { Ok(weight) } else { Err("Invalid molecular formula") }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use crate::molecule::molecular_weight;
    use math::round::half_up;

    fn s(s: &str) -> String {
        String::from(s)
    }

    fn round(weight: f32) -> f64 {
        half_up(weight as f64, 2)
    }

    #[test]
    fn ethane() {
        let molecule: HashMap<String, u32> =
            [(s("C"), 2),
            (s("H"), 6)]
            .iter().cloned().collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(weight, 30.07);
    }

    #[test]
    fn cellulose() {
        let molecule: HashMap<String, u32> =
            [(s("C"), 6),
            (s("H"), 10),
            (s("O"), 5)]
            .iter().cloned().collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(round(weight), 162.14);
    }

    #[test]
    fn vanadium_acetylacetonate() {
        let molecule: HashMap<String, u32> =
            [
                (s("V"), 1),
                (s("C"), 15),
                (s("H"), 21),
                (s("O"), 6)
            ].iter().cloned().collect();
        let weight = molecular_weight(molecule).unwrap();
        assert_eq!(round(weight), 348.27);
    }
}
