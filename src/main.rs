mod model;
mod molecule;
mod parse;
mod test_utils;

#[macro_use]
extern crate log;
extern crate env_logger;

extern crate clap;
use clap::{App, Arg, ArgMatches, SubCommand};

use crate::model::*;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

fn parse_reaction(args: ArgMatches) -> Result<Reaction, String> {
    let in_formula = args
        .value_of("reagent_formula")
        .ok_or("Invalid reagent formula");
    let in_grams_arg = args
        .value_of("reagent_weight")
        .ok_or(format!("Missing reagent weight"))
        .and_then(|g| g.parse::<f32>().or(Err(format!("Invalid reagent weight {}", g))));
    let out_formula = args.value_of("product_formula")
        .ok_or("Invalid product formula");
    let out_grams_arg = args
        .value_of("product_weight")
        .ok_or(format!("Missing product weight"))
        .and_then(|g| g.parse::<f32>().or(Err(format!("Invalid product weight {}", g))));
    let rg = Substance::new(in_formula?, in_grams_arg?)?;
    let pd = Substance::new(out_formula?, out_grams_arg?)?;
    Ok(Reaction::new(rg, pd))
}

fn main() {
    env_logger::Builder::from_env("STOICHKIT_LOG").init();
    let cli: ArgMatches = App::new("StoichKit")
        .version("0.1")
        .arg(Arg::with_name("reagent_formula").required(true).index(1))
        .arg(Arg::with_name("reagent_weight").required(true).index(2))
        .arg(Arg::with_name("product_formula").required(true).index(3))
        .arg(Arg::with_name("product_weight").required(true).index(4))
        .get_matches();
    match parse_reaction(cli).map(|r| r.percent_yield()) {
        Ok(yld) => println!("Yield: {:?}", yld),
        Err(msg) => println!("ERROR: {:?}", msg),
    }
}
