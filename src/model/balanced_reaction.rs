use std::collections::HashMap;
use std::slice::Iter;

use crate::model::{ElementCounts, Reactant};

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
        reactants: &[Reactant],
        products: &[Reactant],
    ) -> Result<(), (ElementCounts, ElementCounts)> {
        let react_elems: ElementCounts = Reactant::element_counts(reactants);
        let prod_elems: ElementCounts = Reactant::element_counts(products);
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

impl From<std::slice::Iter<'_, Reactant>> for ReactantMap {
    fn from(reactants: Iter<'_, Reactant>) -> Self {
        let mut m = HashMap::new();
        for r in reactants {
            m.entry(r.compound.formula.clone())
                .or_insert(r.molar_coefficient);
        }
        Self(m)
    }
}

impl PartialEq for BalancedReaction {
    fn eq(&self, other: &Self) -> bool {
        let ReactantMap(r) = &self.reactants.iter().into();
        let ReactantMap(p) = &self.products.iter().into();
        let ReactantMap(or) = &other.reactants.iter().into();
        let ReactantMap(op) = &other.products.iter().into();
        r == or && p == op
    }
}
