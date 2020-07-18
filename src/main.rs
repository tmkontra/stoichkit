mod molecule;
mod parse;

extern crate clap;
use clap::{App, Arg, SubCommand};

use crate::molecule::molecular_weight;
use crate::parse::parse_formula;

fn main() {
    let matches = App::new("StoichKit")
        .version("0.1")
        .arg(Arg::with_name("reagent_formula").required(true).index(1))
        .arg(Arg::with_name("reagent_weight").required(true).index(2))
        .arg(Arg::with_name("product_formula").required(true).index(3))
        .arg(Arg::with_name("product_weight").required(true).index(4))
        .get_matches();
    let in_formula = matches.value_of("reagent_formula").unwrap();
    let in_grams_arg = matches
        .value_of("reagent_weight")
        .and_then(|g| g.parse::<f32>().ok());
    let out_formula = matches.value_of("product_formula").unwrap();
    let parsed = parse_formula(in_formula).and_then(|atoms| molecular_weight(atoms));
    let parsed_out = parse_formula(out_formula).and_then(|atoms| molecular_weight(atoms));
    let out_grams_arg = matches
        .value_of("product_weight")
        .and_then(|g| g.parse::<f32>().ok());
    match parsed {
        Ok(molwt) => {
            println!("In MolWt: {:?}", molwt);
            match parsed_out {
                Ok(pw) => {
                    println!("Out MolWt: {:?}", pw);
                    println!("Ratio: {:?}", pw / molwt);
                    match (in_grams_arg, out_grams_arg) {
                        (Some(r), Some(p)) => {
                            let yld = (r / molwt) / (p / pw);
                            println!("Yield: {:?}", yld)
                        }
                        _ => println!("Could not calculate"),
                    }
                }
                Err(msg) => println!("ERROR: {:?}", msg),
            }
        }
        Err(msg) => println!("ERROR: {:?}", msg),
    }
}
