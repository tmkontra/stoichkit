use core::result::Result;
use core::result::Result::{Err, Ok};

use itertools::Itertools;

use crate::model::*;
use crate::model::{TheoreticalReaction, YieldReaction};

pub struct ReactionList {
    substances: Vec<String>
}

impl ReactionList {
    pub fn new(substances: Vec<String>) -> ReactionList {
        ReactionList { substances }
    }

    pub fn yield_reaction(&self) -> Result<YieldReaction, String> {
        let (reagent_input, product_input) = self.split_reagents_products();
        let reagents =
            ReactionList::mass_pairs_to_substances(reagent_input)?;
        let products = ReactionList::mass_pairs_to_substances(product_input)?;
        let product = match products.len() {
            0 => Err("Must specify a product!"),
            1 => Ok(products.first().unwrap()),
            _ => Err("Must specify only one product!")
        }?;
        Ok(YieldReaction::new(reagents.to_vec(), product.to_owned()))
    }

    pub fn theoretical_reaction(&self) -> Result<TheoreticalReaction, String> {
        let (reagent_input, product_input) = self.split_reagents_products();
        match (reagent_input.clone().len(), product_input.clone().len()) {
            (x, y) if x >= 1 as usize && y >= 1 as usize => Ok(()),
            _ => Err("Must provide at least 1 reactant and 1 product"),
        }?;
        let reactant_samples: Vec<Sample> =
            ReactionList::mass_pairs_to_substances(reagent_input)?;
        let products: Result<Vec<Reactant>, String> =
            product_input.into_iter()
                .enumerate()
                .map(|(i, p)| ReactionList::str_to_reactant(p, i))
                .collect();
        let reactants = reactant_samples
            .clone()
            .into_iter()
            .map(|r| r.reactant.to_owned())
            .collect();
        let reaction = BalancedReaction::new(reactants, products?)?;
        Ok(TheoreticalReaction::new(reaction, reactant_samples))
    }

    pub fn reaction(&self) -> Result<Reaction, String> {
        let (reagents, products) = self.split_reagents_products();
        let reagents: Result<Vec<Compound>, String> =
            reagents.into_iter().map(|f| Compound::from_formula(f.as_str())).collect();
        let products: Result<Vec<Compound>, String> =
            products.into_iter().map(|f| Compound::from_formula(f.as_str())).collect();
        Reaction::new(reagents?, products?)
    }

    fn split_reagents_products(&self) -> (Vec<String>, Vec<String>) {
        let reagent_input: Vec<String> = self
            .substances
            .clone()
            .into_iter()
            .take_while(|a| a.as_str() != "=")
            .collect();
        let rx_len = reagent_input.clone().len() + 1;
        let product_input: Vec<String> = self
            .substances
            .clone()
            .into_iter()
            .dropping(rx_len)
            .collect();
        (reagent_input, product_input)
    }

    pub fn substance_list(&self) -> Result<Vec<Sample>, String> {
        ReactionList::mass_pairs_to_substances(self.substances.clone())
    }

    fn str_to_reactant(
        formula: String,
        index: usize,
    ) -> Result<Reactant, String> {
        let stoich: Vec<&str> = formula.as_str().split('*').collect();
        let (coeff, formula): (usize, &str) = match stoich.len() {
            1 => (1, formula.as_str()),
            2 => (
                stoich.first().unwrap().parse::<usize>().map_err(|_| {
                    format!(
                        "Invalid coefficient {} for substance {}",
                        stoich.first().unwrap(),
                        formula
                    )
                })?,
                stoich.iter().cloned().last().unwrap(),
            ),
            _ => {
                return Err(format!(
                    "Invalid formula {} at position {}",
                    formula,
                    index + 1
                ))
            }
        };
        Reactant::from_formula(formula, coeff)
    }

    fn collect_mass_pairs(
        substance_strings: Vec<String>,
    ) -> Result<Vec<(String, f32)>, String> {
        substance_strings
            .chunks(2)
            .map(|c| c.to_vec().to_owned())
            .map(|pair| {
                if pair.len() < 2 {
                    return Err(format!(
                        "Got substance with no mass: {}",
                        pair[0]
                    ));
                } else {
                    match pair[1].parse::<f32>() {
                        Ok(mass) => Ok((pair[0].to_owned(), mass)),
                        Err(e) => Err(format!("Expected mass for '{}', but could not parse '{}' as float.", pair[0], pair[1]))
                    }
                }
            })
            .collect()
    }

    fn mass_pairs_to_substances(
        substance_strings: Vec<String>,
    ) -> Result<Vec<Sample>, String> {
        let mut substances = vec![];
        for (i, pair) in ReactionList::collect_mass_pairs(substance_strings)?
            .into_iter()
            .enumerate()
        {
            let reactant = ReactionList::str_to_reactant(pair.0.to_owned(), i)?;
            let substance = Sample::of_reactant(
                reactant, pair.1,
            );
            substances.push(substance);
        }
        Ok(substances)
    }
}
