mod model;
mod molecule;
mod parse;
mod test_utils;

#[macro_use]
extern crate log;
extern crate env_logger;

extern crate clap;
use clap::{App, Arg, ArgMatches, Clap};

use crate::model::*;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

#[derive(Clap)]
#[clap(version = "0.1.0")]
struct ReactionList {
    substances: Vec<String>,
}

impl ReactionList {
    pub fn reaction(&self) -> Result<Reaction, String> {
        let mut substances: Vec<Substance> = vec![];
        for (i, pair) in self.substances.chunks(2).map(|c| c.to_vec()).enumerate() {
            if pair.len() < 2 {
                Err(format!("Got substance with no mass at position {}", i + 1))?;
            }
            let stoich: Vec<&str> = pair[0].as_str().split('*').collect();
            let (coeff, formula) = match stoich.len() {
                1 => (None, pair[0].as_str()),
                2 => (
                    Some(stoich.first().unwrap().parse::<u32>().map_err(|_| {
                        format!(
                            "Invalid coeffcieint {} for substance {}",
                            stoich.first().unwrap(),
                            pair[0]
                        )
                    })?),
                    stoich.last().unwrap().clone(),
                ),
                _ => Err(format!("Invalid formula {} at position {}", pair[0], i + 1).to_string())?,
            };
            let substance = Substance::new(
                formula,
                pair[1]
                    .clone()
                    .parse::<f32>()
                    .map_err(|_| format!("Invalid mass {} for substance {}", pair[1], pair[0]))?,
                coeff,
            )?;
            substances.push(substance);
        }
        let (product, reagents) = substances
            .split_last()
            .ok_or(format!("Invalid substance list!"))?;
        Ok(Reaction::new(reagents.to_vec(), product.to_owned()))
    }
}

fn main() {
    env_logger::Builder::from_env("STOICHKIT_LOG").init();
    let opts: ReactionList = ReactionList::parse();
    match opts.reaction().map(|r| r.percent_yield()) {
        Ok(yld) => println!("Yield: {:?}", yld),
        Err(msg) => println!("ERROR: {:?}", msg),
    }
}
