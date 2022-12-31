use crate::model::{Compound, Element};
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Reactant {
    pub compound: Compound,
    pub molar_coefficient: usize,
}

impl Reactant {
    pub fn of_compound(compound: Compound, coefficient: usize) -> Self {
        Reactant {
            compound,
            molar_coefficient: coefficient,
        }
    }

    pub fn from_formula(
        formula: &str,
        coefficient: usize,
    ) -> Result<Self, String> {
        let cmp = Compound::from_formula(formula);
        Ok(Reactant {
            compound: cmp?,
            molar_coefficient: coefficient,
        })
    }

    pub fn fold_elements(reactants: &Vec<Reactant>) -> HashMap<Element, usize> {
        return reactants
            .iter()
            .map(|s| (&s.compound.atoms, s.molar_coefficient))
            .fold(HashMap::new(), |mut acc, (item, coeff)| {
                for (e, c) in item {
                    let counter = acc.entry(e.to_owned()).or_insert(0);
                    *counter += *c * coeff;
                }
                acc
            });
    }

    pub fn format(&self, explicit: bool) -> String {
        if self.molar_coefficient != 1 || explicit {
            format!("{}*{}", self.molar_coefficient, self.compound.formula)
        } else {
            self.compound.formula.to_string()
        }
    }
}
