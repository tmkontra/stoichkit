use std::cmp::{Eq, Ordering};
use std::collections::HashMap;
use std::{hash::Hash, iter::FromIterator};

use periodic_table_on_an_enum::Element as PElement;

use crate::molecule::molecular_weight;
use crate::parse::parse_formula;

#[derive(Clone, Debug)]
pub struct Compound {
    pub formula: String,
    pub atoms: HashMap<Element, u32>,
    molar_mass: f32,
}

impl Compound {
    pub fn from_formula(formula: &str) -> Result<Compound, String> {
        Compound::new(formula)
    }

    pub fn new(formula: &str) -> Result<Compound, String> {
        let atoms = parse_formula(formula);
        let molecular_weight = atoms.clone().and_then(molecular_weight)?;
        Ok(Compound {
            formula: formula.to_string(),
            atoms: atoms?,
            molar_mass: molecular_weight,
        })
    }
}

#[derive(Clone, Debug)]
pub struct Reactant {
    pub compound: Compound,
    pub molar_coefficient: u32,
}

impl Reactant {
    pub fn of_compound(compound: Compound, coeff: u32) -> Self {
        Reactant {
            compound,
            molar_coefficient: coeff,
        }
    }

    pub fn from_formula(formula: &str, coeff: u32) -> Result<Self, String> {
        let cmp = Compound::from_formula(formula);
        Ok(Reactant {
            compound: cmp?,
            molar_coefficient: coeff,
        })
    }
}

#[derive(Clone, Debug)]
pub struct Substance {
    pub reactant: Reactant,
    pub mass: f32,
}

impl Substance {
    pub fn of_reactant(reactant: Reactant, mass: f32) -> Self {
        Substance { reactant, mass }
    }

    pub fn from_formula(
        formula: &str,
        mass: f32,
        molar_coefficient: u32,
    ) -> Result<Substance, String> {
        let rct = Reactant::from_formula(formula, molar_coefficient);
        Ok(Substance {
            reactant: rct?,
            mass,
        })
    }

    pub fn moles(self: &Self) -> f32 {
        self.mass / self.reactant.compound.molar_mass
    }

    pub fn molrxn(self: &Self) -> f32 {
        self.moles() / self.reactant.molar_coefficient as f32
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Element {
    element: PElement,
}

impl Element {
    pub fn from_symbol(sym: &str) -> Option<Element> {
        PElement::from_symbol(sym).map(Element::from_pt_element)
    }

    pub fn from_atomic_number(z: usize) -> Option<Element> {
        PElement::from_atomic_number(z).map(Element::from_pt_element)
    }

    fn from_pt_element(element: PElement) -> Element {
        Element { element }
    }

    pub fn get_atomic_mass(&self) -> f32 {
        self.element.get_atomic_mass()
    }
}

#[derive(Debug, Clone)]
pub struct BalancedReaction {
    pub reactants: Vec<Reactant>,
    pub products: Vec<Reactant>,
}

struct ReactantMap(HashMap<String, u32>);

impl FromIterator<Reactant> for ReactantMap {
    fn from_iter<R>(reactants: R) -> Self
    where
        R: IntoIterator<Item = Reactant>,
    {
        let mut m = HashMap::new();
        for r in reactants {
            m.entry(r.compound.formula).or_insert(r.molar_coefficient);
        }
        Self(m)
    }
}

impl PartialEq for BalancedReaction {
    fn eq(&self, other: &Self) -> bool {
        let ReactantMap(r) = self.reactants.clone().into_iter().collect();
        let ReactantMap(p) = self.products.clone().into_iter().collect();
        let ReactantMap(or) = other.reactants.clone().into_iter().collect();
        let ReactantMap(op) = other.products.clone().into_iter().collect();
        r == or && p == op
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

impl BalancedReaction {
    pub fn new(reactants: Vec<Reactant>, products: Vec<Reactant>) -> Self {
        BalancedReaction {
            reactants,
            products,
        }
    }

    fn reactants_display_string(&self) -> String {
        let balr: Vec<String> = self
            .reactants
            .iter()
            .map(|r| format!("{} {}", r.molar_coefficient, r.compound.formula))
            .collect();
        balr.join(" + ")
    }

    fn products_display_string(&self) -> String {
        let balp: Vec<String> = self
            .products
            .iter()
            .map(|r| format!("{} {}", r.molar_coefficient, r.compound.formula))
            .collect();
        balp.join(" + ")
    }

    pub fn display_string(&self) -> String {
        format!(
            "{} = {}",
            self.reactants_display_string(),
            self.products_display_string()
        )
    }
}

#[derive(Debug, Clone)]
pub struct YieldReaction {
    pub reagents: Vec<Substance>,
    pub product: Substance,
}

impl YieldReaction {
    pub fn new(reagents: Vec<Substance>, product: Substance) -> YieldReaction {
        YieldReaction { reagents, product }
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
                debug!("Limiting reagent is {}", s.reactant.compound.formula);
                s
            })
            .unwrap()
    }

    pub fn theoretical_yield(self: &Self) -> f32 {
        let limiting = self.limiting_reagent();
        trace!("{} moles of limiting reagent", limiting.moles());
        let exp_moles = limiting.moles()
            * (self.product.reactant.molar_coefficient as f32
                / limiting.reactant.molar_coefficient as f32);
        debug!("Theoretical moles of product: {}", exp_moles);
        let exp_grams = exp_moles * self.product.reactant.compound.molar_mass;
        debug!("Theoretical yield of product (g): {}", exp_grams);
        exp_grams
    }

    pub fn percent_yield(self: &Self) -> f32 {
        self.product.mass / self.theoretical_yield()
    }
}
