extern crate clap;
#[macro_use]
extern crate log;

use std::fs::read_to_string;

use clap::Clap;
use itertools::Itertools;

use stoichkit::ext::parse_chemdraw_reaction;
use stoichkit::model::Reaction;
use stoichkit::model::Sample;
use stoichkit::model::{BalancedReaction, Compound, Reactant};
use stoichkit::model::{TheoreticalReaction, YieldReaction};

#[derive(Clap)]
#[clap(version = "0.2.0")]
struct Cli {
    #[clap(subcommand)]
    command: Subcommand,
}

#[derive(Clap)]
struct ReactionList {
    #[clap(
        about = "Fully balanced chemical reaction list: [<coeff>*]<formula> <grams>, coeff defaults to 1"
    )]
    substances: Vec<String>,
}

#[derive(Clap)]
enum Subcommand {
    TheoreticalYield(ReactionList),
    Yield(ReactionList),
    Balance(Unbal),
}

#[derive(Clap)]
struct Unbal {
    #[clap(about = "Chemical equation [...reactants] = [...products]")]
    substances: Vec<String>,
    #[clap(short, conflicts_with = "substances")]
    chemdraw_file: Option<String>,
}

impl Unbal {
    pub fn balance(&self) -> Result<String, String> {
        let (reagents, products) = match &self.chemdraw_file {
            None => {
                let reagent_input: Result<Vec<Compound>, _> = self
                    .substances
                    .iter()
                    .take_while(|a| a.as_str() != "=")
                    .map(|f| Compound::from_formula(f.as_str()))
                    .collect();
                let rx_len = reagent_input.clone()?.len() + 1;
                let product_input: Result<Vec<Compound>, _> = self
                    .substances
                    .iter()
                    .dropping(rx_len)
                    .map(|f| Compound::from_formula(f.as_str()))
                    .collect();
                match (
                    reagent_input.clone()?.len(),
                    product_input.clone()?.len(),
                ) {
                    (x, y) if x >= 1 as usize && y >= 1 as usize => Ok(()),
                    _ => Err("Must provide at least 1 reactant and 1 product"),
                }?;
                (reagent_input?, product_input?)
            }
            Some(f) => {
                let s = read_to_string(f)
                    .map_err(|_e| format!("Could not read file {:?}", f))?;
                let result = parse_chemdraw_reaction(s.as_str())?;
                info!(
                    "Parsed reaction {:?} = {:?}",
                    result.reactants, result.products
                );
                (result.reactants, result.products)
            }
        };
        let rxn = Reaction::new(reagents, products)?;
        let balanced = rxn.balance()?;
        Ok(balanced.display_string())
    }
}

impl ReactionList {
    pub fn yield_reaction(&self) -> Result<YieldReaction, String> {
        let substances =
            ReactionList::mass_pairs_to_substances(self.substances.clone())?;
        let (product, reagents) = substances
            .split_last()
            .ok_or_else(|| "Invalid substance list!".to_string())?;
        Ok(YieldReaction::new(reagents.to_vec(), product.to_owned()))
    }

    pub fn theoretical_reaction(&self) -> Result<TheoreticalReaction, String> {
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
    ) -> Result<Vec<(String, String)>, String> {
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
                    Ok((pair[0].to_owned(), pair[1].to_owned()))
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
                reactant,
                pair.1.clone().parse::<f32>().map_err(|_| {
                    format!("Invalid mass {} for substance {}", pair.1, pair.0)
                })?,
            );
            substances.push(substance);
        }
        Ok(substances)
    }
}

fn main() {
    env_logger::Builder::from_env("STOICHKIT_LOG").init();
    let opts: Cli = Cli::parse();
    match opts.command {
        Subcommand::TheoreticalYield(r) => {
            match r.theoretical_reaction().map(|r| r.yields()) {
                Ok(yields) => yields.iter().for_each(|y| println!("{}", y)),
                Err(msg) => println!("ERROR: {:?}", msg),
            }
        }
        Subcommand::Yield(r) => {
            match r.yield_reaction().map(|r| r.percent_yield()) {
                Ok(yld) => println!("Yield: {:?}", yld),
                Err(msg) => println!("ERROR: {:?}", msg),
            }
        }
        Subcommand::Balance(u) => {
            let result = u.balance().unwrap_or_else(|e| e);
            println!("{}", result);
        }
    }
}
