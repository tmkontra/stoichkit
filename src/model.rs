use std::cmp::Ordering;
use std::collections::HashMap;

use ptable::Element;

use crate::parse::parse_formula;

#[derive(Clone, Debug)]
pub struct Reactant {
    pub substance: Substance,
    pub molar_coefficient: u32,
}

impl Reactant {
    pub fn new(
        formula: &str,
        molar_coefficient: u32,
    ) -> Result<Reactant, String> {
        Ok(Reactant {
            substance: Substance::new(formula)?,
            molar_coefficient,
        })
    }
}

#[derive(Clone, Debug)]
pub struct Reagent {
    pub substance: Substance,
    pub molar_coefficient: u32,
    pub grams: f32,
}

impl Reagent {
    pub fn new(
        formula: &str,
        molar_coefficient: u32,
        grams: f32,
    ) -> Result<Reagent, String> {
        Ok(Reagent {
            substance: Substance::new(formula)?,
            molar_coefficient,
            grams,
        })
    }

    pub fn moles(&self) -> f32 {
        self.grams / self.substance.molecular_weight()
    }

    pub fn molrxn(&self) -> f32 {
        self.moles() / self.molar_coefficient as f32
    }
}

#[derive(Clone, Debug)]
pub struct Substance {
    pub formula: String,
    pub atoms: HashMap<Element, u32>,
}

impl Substance {
    pub fn new(formula: &str) -> Result<Substance, String> {
        let atoms = parse_formula(formula)?;
        Ok(Substance {
            formula: formula.to_string(),
            atoms,
        })
    }

    pub fn molecular_weight(&self) -> f32 {
        let mut weight: f32 = 0 as f32;
        for (element, count) in &self.atoms {
            let mass = element.get_atomic_mass();
            trace!("Adding {:?} x {:?} for element {:?}", count, mass, element);
            weight += mass * *count as f32
        }
        weight
    }
}

#[derive(Debug)]
pub struct UnbalancedReaction {
    pub reactants: Vec<Substance>,
    pub products: Vec<Substance>,
}

pub struct YieldReaction {
    pub reagents: Vec<Reagent>,
    pub product: Reagent,
}

impl YieldReaction {
    pub fn new(reagents: Vec<Reagent>, product: Reagent) -> YieldReaction {
        YieldReaction { reagents, product }
    }

    pub fn limiting_reagent(self: &Self) -> &Reagent {
        self.reagents
            .iter()
            .min_by(|l, r| {
                l.molrxn()
                    .partial_cmp(&r.molrxn())
                    .unwrap_or(Ordering::Equal)
            })
            .map(|r| {
                debug!("Limiting reagent is {}", r.substance.formula);
                r
            })
            .unwrap()
    }

    pub fn theoretical_yield(self: &Self) -> f32 {
        let limiting = self.limiting_reagent();
        trace!("{} moles of limiting reagent", limiting.moles());
        let exp_moles = limiting.moles()
            * (self.product.molar_coefficient as f32
                / limiting.molar_coefficient as f32);
        debug!("Theoretical moles of product: {}", exp_moles);
        let exp_grams = exp_moles * self.product.substance.molecular_weight();
        debug!("Theoretical yield of product (g): {}", exp_grams);
        exp_grams
    }

    pub fn percent_yield(self: &Self) -> f32 {
        self.product.grams / self.theoretical_yield()
    }
}
