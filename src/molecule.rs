use std::collections::HashMap;

use ptable::Element;

pub fn molecular_weight(atoms: HashMap<Element, u32>) -> Result<f32, String> {
    let mut weight: f32 = 0 as f32;
    for (element, count) in atoms {
        let mass = element.get_atomic_mass();
        trace!("Adding {:?} x {:?} for element {:?}", count, mass, element);
        weight += mass * count as f32
    }
    Ok(weight)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use math::round::half_up;

    use crate::molecule::molecular_weight;
    use crate::test_utils::e;

    fn round(weight: f32) -> f64 {
        half_up(weight as f64, 2)
    }

    #[test]
    fn ethane() {
        let molecule: HashMap<&str, u32> = [("C", 2), ("H", 6)].iter().cloned().collect();
        let weight = molecular_weight(e(molecule)).unwrap();
        assert_eq!(weight, 30.07);
    }

    #[test]
    fn cellulose() {
        let molecule: HashMap<&str, u32> =
            [("C", 6), ("H", 10), ("O", 5)].iter().cloned().collect();
        let weight = molecular_weight(e(molecule)).unwrap();
        assert_eq!(round(weight), 162.14);
    }

    #[test]
    fn vanadium_acetylacetonate() {
        let molecule: HashMap<&str, u32> = [("V", 1), ("C", 15), ("H", 21), ("O", 6)]
            .iter()
            .cloned()
            .collect();
        let weight = molecular_weight(e(molecule)).unwrap();
        assert_eq!(round(weight), 348.27);
    }
}
