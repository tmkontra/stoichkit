mod molecule;
mod parse;
mod test_utils;

extern crate clap;
use clap::{App, Arg, ArgMatches, SubCommand};

use crate::molecule::molecular_weight;
use crate::parse::parse_formula;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

fn percent_yield(reaction: &Reaction) -> f32 {
    let rmoles = reaction.reagent.moles();
    let pmoles = reaction.product.moles();
    rmoles / pmoles
}

struct Substance {
    formula: String,
    mass: f32,
    atoms: HashMap<String, u32>,
    molecular_weight: f32,
}

impl Substance {
    pub fn new(formula: &str, mass: f32) -> Result<Substance, String> {
        let atoms = parse_formula(formula);
        let molecular_weight = atoms
            .clone()
            .and_then(|atoms| molecular_weight(atoms));
        match molecular_weight {
            Ok(wt) => Ok(Substance {
                formula: formula.to_string(),
                mass,
                atoms: atoms?,
                molecular_weight: wt,
            }),
            Err(_) => Err(String::from("Could not parse formula")),
        }
    }

    pub fn moles(self: &Self) -> f32 {
        self.mass / self.molecular_weight
    }
}

struct Reaction {
    pub reagent: Substance,
    pub product: Substance,
}

fn parse_reaction(args: ArgMatches) -> Option<Reaction> {
    let in_formula = args.value_of("reagent_formula");
    let in_grams_arg = args
        .value_of("reagent_weight")
        .and_then(|g| g.parse::<f32>().ok());
    let out_formula = args.value_of("product_formula");
    let out_grams_arg = args
        .value_of("product_weight")
        .and_then(|g| g.parse::<f32>().ok());
    let rg = Substance::new(in_formula?, in_grams_arg?).ok()?;
    let pd = Substance::new(out_formula?, out_grams_arg?).ok()?;
    Some(Reaction {
        reagent: rg,
        product: pd,
    })
}

fn main() {
    let matches: ArgMatches = App::new("StoichKit")
        .version("0.1")
        .arg(Arg::with_name("reagent_formula").required(true).index(1))
        .arg(Arg::with_name("reagent_weight").required(true).index(2))
        .arg(Arg::with_name("product_formula").required(true).index(3))
        .arg(Arg::with_name("product_weight").required(true).index(4))
        .get_matches();
    let reaction = parse_reaction(matches);
    let yld = reaction.as_ref().map(|r| percent_yield(r));
    println!("Yield: {:?}", yld)
}
