use std::cmp::Ordering;
use std::collections::HashMap;

use ptable::Element;

use crate::molecule::molecular_weight;
use crate::parse::parse_formula;

#[derive(Clone, Debug)]
pub struct Substance {
    pub formula: String,
    mass: f32,
    pub atoms: HashMap<Element, u32>,
    molar_mass: f32,
    molar_coefficient: u32,
}

impl Substance {
    pub fn from_formula(formula: &str) -> Result<Substance, String> {
        Substance::new(formula, 0., None)
    }

    pub fn new(
        formula: &str,
        mass: f32,
        molar_coefficient: Option<u32>,
    ) -> Result<Substance, String> {
        let atoms = parse_formula(formula);
        let molecular_weight = atoms.clone().and_then(molecular_weight)?;
        Ok(Substance {
            formula: formula.to_string(),
            mass,
            atoms: atoms?,
            molar_mass: molecular_weight,
            molar_coefficient: molar_coefficient.unwrap_or(1),
        })
    }

    pub fn moles(self: &Self) -> f32 {
        self.mass / self.molar_mass
    }

    pub fn molrxn(self: &Self) -> f32 {
        self.moles() / self.molar_coefficient as f32
    }
}

pub struct Reaction {
    pub reagents: Vec<Substance>,
    pub product: Substance,
}

impl Reaction {
    pub fn new(reagents: Vec<Substance>, product: Substance) -> Reaction {
        Reaction { reagents, product }
    }

    pub fn limiting_reagent(self: &Self) -> &Substance {
        self.reagents
            .iter()
            .min_by(|l, r| {
                l.molrxn()
                    .partial_cmp(&r.molrxn())
                    .unwrap_or(Ordering::Equal)
            })
            .map(|s| {
                debug!("Limiting reagent is {}", s.formula);
                s
            })
            .unwrap()
    }

    pub fn theoretical_yield(self: &Self) -> f32 {
        let limiting = self.limiting_reagent();
        trace!("{} moles of limiting reagent", limiting.moles());
        let exp_moles = limiting.moles()
            * (self.product.molar_coefficient as f32 / limiting.molar_coefficient as f32);
        debug!("Theoretical moles of product: {}", exp_moles);
        let exp_grams = exp_moles * self.product.molar_mass;
        debug!("Theoretical yield of product (g): {}", exp_grams);
        exp_grams
    }

    pub fn percent_yield(self: &Self) -> f32 {
        self.product.mass / self.theoretical_yield()
    }
}
