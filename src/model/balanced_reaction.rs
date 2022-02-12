use std::collections::HashMap;
use std::iter::FromIterator;

use crate::model::{Element, Reactant};

#[derive(Debug, Clone)]
pub struct BalancedReaction {
    pub reactants: Vec<Reactant>,
    pub products: Vec<Reactant>,
}

impl BalancedReaction {
    pub fn new(
        reactants: Vec<Reactant>,
        products: Vec<Reactant>,
    ) -> Result<BalancedReaction, String> {
        match BalancedReaction::check_balance(&reactants, &products) {
            Ok(_) => Ok(BalancedReaction {
                reactants,
                products,
            }),
            Err((reactants, products)) => Err(format!(
                "Equation is not be balanced\n{:?}\n{:?}",
                reactants, products
            )),
        }
    }

    fn check_balance(
        reactants: &Vec<Reactant>,
        products: &Vec<Reactant>,
    ) -> Result<(), (HashMap<Element, usize>, HashMap<Element, usize>)> {
        let react_elems: HashMap<Element, usize> =
            Reactant::fold_elements(&reactants);
        let prod_elems: HashMap<Element, usize> =
            Reactant::fold_elements(&products);
        debug!(
            "Checking balanced?: Reagent elements: {:?} === Product elements: {:?}",
            react_elems, prod_elems
        );
        match react_elems.eq(&prod_elems) {
            true => Ok(()),
            false => Err((react_elems, prod_elems)),
        }
    }

    fn reactants_display_string(&self, explicit: bool) -> String {
        let balanced: Vec<String> =
            self.reactants.iter().map(|r| r.format(explicit)).collect();
        balanced.join(" + ")
    }

    fn products_display_string(&self, explicit: bool) -> String {
        let balanced: Vec<String> =
            self.products.iter().map(|r| r.format(explicit)).collect();
        balanced.join(" + ")
    }

    pub fn display_string(&self, explicit: bool) -> String {
        format!(
            "{} = {}",
            self.reactants_display_string(explicit),
            self.products_display_string(explicit)
        )
    }

    #[allow(dead_code)]
    pub(crate) fn all_coefficients(&self) -> Vec<usize> {
        self.reactants
            .iter()
            .chain(self.products.iter())
            .map(|r| r.molar_coefficient)
            .collect()
    }
}

struct ReactantMap(HashMap<String, usize>);

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
