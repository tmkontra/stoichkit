extern crate clap;

use clap::Clap;
use itertools::Itertools;

use stoichkit::model::{Reaction, Substance};
use stoichkit::solve::balance;

#[derive(Clap)]
struct ReactionList {
    #[clap(
        about = "Fully balanced chemical reaction list: [<coeff>*]<formula> <grams>, coeff defaults to 1"
    )]
    substances: Vec<String>,
}

#[derive(Clap)]
enum Subcommand {
    Yield(ReactionList),
    Balance(Unbal),
}

#[derive(Clap)]
struct Unbal {
    #[clap(about = "Chemical equation [...reactants] = [...products]")]
    substances: Vec<String>,
}

impl Unbal {
    pub fn balance(&self) -> Result<String, String> {
        let r: Result<Vec<Substance>, _> = self
            .substances
            .iter()
            .take_while(|a| a.as_str() != "=")
            .map(|f| Substance::from_formula(f.as_str()))
            .collect();
        let rx_len = r.clone()?.len() + 1;
        let p: Result<Vec<Substance>, _> = self
            .substances
            .iter()
            .dropping(rx_len)
            .map(|f| Substance::from_formula(f.as_str()))
            .collect();
        match (r.clone()?.len(), p.clone()?.len()) {
            (x, y) if x >= 1 as usize && y >= 1 as usize => Ok(()),
            _ => Err("Must provide at least 1 reactant and 1 product"),
        }?;
        let balanced = balance(r?, p?)?;
        let balr: Vec<String> = balanced
            .iter()
            .take(rx_len - 1)
            .map(|(e, c)| format!("{} {}", c, e))
            .collect();
        let balp: Vec<String> = balanced
            .iter()
            .dropping(rx_len - 1)
            .map(|(e, c)| format!("{} {}", c, e))
            .collect();
        let result = format!("{} = {}", balr.join(" + "), balp.join(" + "));
        Ok(result)
    }
}

#[derive(Clap)]
#[clap(version = "0.2.3")]
struct Cli {
    #[clap(subcommand)]
    subcmd: Subcommand,
    //
    // #[clap(subcommand)]
    // bal: Unbal
}

impl ReactionList {
    pub fn reaction(&self) -> Result<Reaction, String> {
        let mut substances: Vec<Substance> = vec![];
        for (i, pair) in
            self.substances.chunks(2).map(|c| c.to_vec()).enumerate()
        {
            if pair.len() < 2 {
                return Err(format!(
                    "Got substance with no mass at position {}",
                    i + 1
                ));
            }
            let stoich: Vec<&str> = pair[0].as_str().split('*').collect();
            let (coeff, formula): (Option<u32>, &str) = match stoich.len() {
                1 => (None, pair[0].as_str()),
                2 => (
                    Some(stoich.first().unwrap().parse::<u32>().map_err(
                        |_| {
                            format!(
                                "Invalid coeffcieint {} for substance {}",
                                stoich.first().unwrap(),
                                pair[0]
                            )
                        },
                    )?),
                    stoich.iter().cloned().last().unwrap(),
                ),
                _ => {
                    return Err(format!(
                        "Invalid formula {} at position {}",
                        pair[0],
                        i + 1
                    ))
                }
            };
            let substance = Substance::new(
                formula,
                pair[1].clone().parse::<f32>().map_err(|_| {
                    format!(
                        "Invalid mass {} for substance {}",
                        pair[1], pair[0]
                    )
                })?,
                coeff,
            )?;
            substances.push(substance);
        }
        let (product, reagents) = substances
            .split_last()
            .ok_or_else(|| "Invalid substance list!".to_string())?;
        Ok(Reaction::new(reagents.to_vec(), product.to_owned()))
    }
}

fn main() {
    env_logger::Builder::from_env("STOICHKIT_LOG").init();
    let opts: Cli = Cli::parse();
    match opts.subcmd {
        Subcommand::Yield(r) => match r.reaction().map(|r| r.percent_yield()) {
            Ok(yld) => println!("Yield: {:?}", yld),
            Err(msg) => println!("ERROR: {:?}", msg),
        },
        Subcommand::Balance(u) => {
            let result = u.balance().unwrap_or_else(|e| e);
            println!("{}", result);
        }
    }
}
