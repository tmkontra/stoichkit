use std::collections::HashSet;

use crate::model::{BalancedReaction, Compound, Element, Reactant};
use crate::solve;

#[derive(Debug, Clone)]
pub struct Reaction {
    pub reactants: Vec<Compound>,
    pub products: Vec<Compound>,
}

impl Reaction {
    pub fn new(
        reactants: Vec<Compound>,
        products: Vec<Compound>,
    ) -> Result<Reaction, String> {
        let rxn = Reaction {
            reactants,
            products,
        };
        rxn.check_elements()?;
        Ok(rxn)
    }

    pub fn len(&self) -> usize {
        self.reactants.len() + self.products.len()
    }

    fn check_elements(&self) -> Result<(), String> {
        let reagent_atoms: HashSet<&Element> =
            Reaction::elements_from(&self.reactants);
        let product_atoms: HashSet<&Element> =
            Reaction::elements_from(&self.products);
        return if !&reagent_atoms.eq(&product_atoms) {
            let missing_products: HashSet<_> =
                reagent_atoms.difference(&product_atoms).collect();
            let missing_reagents: HashSet<_> =
                product_atoms.difference(&reagent_atoms).collect();
            Err(format!(
                "Equation cannot be balanced. Reagent elements that are not in products = {:?}. Product elements that are not in products = {:?}",
                missing_products, missing_reagents)
            )
        } else {
            Ok(())
        };
    }

    pub fn balance(&self) -> Result<BalancedReaction, String> {
        let mx = solve::build_matrix(
            &self.all_elements(),
            self.all_compounds(),
            self.len(),
        );
        let coefficients: Vec<f64> = solve::solve_system(mx, self.len() - 1)?;
        debug!("Got solution coefficients: {:?}", &coefficients);
        trace!("Converting to rationals");
        let scaled_coefficients: Vec<usize> =
            solve::normalize_coefficients(coefficients)?;
        let result: Vec<Reactant> = self
            .all_compounds()
            .into_iter()
            .zip(&mut scaled_coefficients.iter().map(|c| c.to_owned()))
            .map(|(c, coefficient)| {
                Reactant::of_compound(c.to_owned(), coefficient)
            })
            .collect();
        Reaction::check_all_nonzero(&result)?;
        let (reagents_result, products_result) =
            result.split_at(self.reactants.len());
        BalancedReaction::new(
            reagents_result.to_vec(),
            products_result.to_vec(),
        )
    }

    fn check_all_nonzero(reactants: &Vec<Reactant>) -> Result<(), String> {
        let zeroes = reactants
            .iter()
            .filter(|r| r.molar_coefficient == 0)
            .map(|c| c.compound.formula.as_str())
            .collect::<Vec<&str>>();
        if !zeroes.is_empty() {
            let formulas = zeroes.join(", ");
            let err_msg = format!(
                "0 coefficient is not a valid solution! Got 0 for: {}",
                formulas
            );
            return Err(err_msg);
        }
        Ok(())
    }

    pub fn all_elements(&self) -> Vec<&Element> {
        Reaction::elements_from(&self.reactants)
            .iter()
            .chain(Reaction::elements_from(&self.products).iter())
            .cloned()
            .collect()
    }

    pub fn all_compounds(&self) -> Vec<&Compound> {
        self.reactants.iter().chain(self.products.iter()).collect()
    }

    fn elements_from(compounds: &Vec<Compound>) -> HashSet<&Element> {
        compounds.iter().flat_map(|c| c.all_elements()).collect()
    }
}
