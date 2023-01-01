use std::collections::HashMap;

use crate::model::Element;
use crate::parse;

pub type ElementCounts = HashMap<Element, usize>;

#[derive(Clone, Debug)]
pub struct Compound {
    pub formula: String,
    pub atoms: ElementCounts,
    pub molar_mass: f32,
}

impl Compound {
    pub fn from_formula(formula: &str) -> Result<Compound, String> {
        Compound::new(formula)
    }

    pub fn new(formula: &str) -> Result<Compound, String> {
        let atoms: HashMap<Element, u64> = parse::parse_formula_v2(formula)?;
        let atoms = atoms
            .into_iter()
            .map(|(k, v)| (k, v as usize))
            .collect();
        let molecular_weight: f32 = Compound::molecular_weight(&atoms);
        Ok(Compound {
            formula: formula.to_string(),
            atoms,
            molar_mass: molecular_weight,
        })
    }

    fn molecular_weight(atoms: &ElementCounts) -> f32 {
        atoms.iter().fold(0 as f32, |acc, (e, count)| {
            acc + e.get_atomic_mass() * count.to_owned() as f32
        })
    }

    pub fn all_elements(&self) -> Vec<&Element> {
        return self.atoms.keys().collect();
    }
}

#[cfg(test)]
mod tests {
    use math::round::half_up;

    use crate::model::compound::Compound;

    fn round(weight: f32) -> f64 {
        half_up(weight as f64, 2)
    }

    #[test]
    fn ethane() {
        let compound = Compound::from_formula("C2H6").unwrap();
        let weight = compound.molar_mass;
        assert_eq!(weight, 30.07);
    }

    #[test]
    fn cellulose() {
        let compound = Compound::from_formula("C6H10O5").unwrap();
        let weight = compound.molar_mass;
        assert_eq!(round(weight), 162.14);
    }

    #[test]
    fn vanadium_acetylacetonate() {
        let compound = Compound::from_formula("V1C15H21O6").unwrap();
        let weight = compound.molar_mass;
        assert_eq!(round(weight), 348.27);
    }
}
