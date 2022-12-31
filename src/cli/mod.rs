use clap::{Args, Parser, Subcommand};

use crate::ext::chemdraw;
use crate::model::{BalancedReaction, Compound, Units, YieldUnits};
use crate::model::ReactionList;

#[derive(Parser)]
#[clap(name = "stoichkit")]
#[clap(about = "A stoichiometry toolkit.", long_about = None, version = "0.6.0")]
pub struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

impl Cli {
    pub fn run(self) {
        let result = match self.command {
            Commands::TheoreticalYield(TheoreticalArgs {
                substances,
                units,
            }) => {
                let units = units.unwrap_or(YieldUnits::Mass);
                ReactionList::new(substances)
                    .theoretical_reaction()
                    .map(|r| {
                        ResultList::print(
                            r.yields(&units)
                                .into_iter()
                                .map(|(r, amt)| (r.compound, amt))
                                .collect(),
                            units.into(),
                        )
                    })
            }
            Commands::Yield(YieldArgs { substances }) => {
                ReactionList::new(substances).yield_reaction().map(|yld| {
                    ResultList::print(
                        vec![(
                            yld.product.reactant.compound.to_owned(),
                            yld.percent_yield(),
                        )],
                        Units::Percent,
                    )
                })
            }
            Commands::Balance(args) => args.balance().map(|balanced| {
                println!("{}", balanced.display_string(args.explicit))
            }),
            Commands::Moles { substances } => {
                ReactionList::new(substances).substance_list().map(|subs| {
                    subs.iter().for_each(|s| {
                        println!(
                            "{} mol",
                            s.mass / s.reactant.compound.molar_mass
                        )
                    })
                })
            }
        };
        match result {
            Ok(_) => (),
            Err(err) => println!("ERROR: {}", err),
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    TheoreticalYield(TheoreticalArgs),
    Yield(YieldArgs),
    Balance(UnbalancedEquationArgs),
    Moles { substances: Vec<String> }, // Pvnrt(GasArgs) TODO,
}

#[derive(Args)]
struct TheoreticalArgs {
    substances: Vec<String>,
    #[clap(short, long, arg_enum)]
    units: Option<YieldUnits>,
}

#[derive(Args)]
struct YieldArgs {
    substances: Vec<String>,
}

struct ResultList {}

impl ResultList {
    pub(crate) fn print(list: Vec<(Compound, f32)>, units: Units) {
        list.iter().for_each(|(product, yld)| {
            println!("{} {} {}", product.formula, yld, units)
        })
    }
}

#[derive(Args)]
struct UnbalancedEquationArgs {
    #[clap(help = "Chemical equation [...reactants] = [...products]")]
    substances: Vec<String>,
    #[clap(short, conflicts_with = "substances")]
    chemdraw_file: Option<String>,
    #[clap(short = 'x', long)]
    explicit: bool,
}

impl UnbalancedEquationArgs {
    pub fn balance(&self) -> Result<BalancedReaction, String> {
        let rxn = match &self.chemdraw_file {
            None => ReactionList::new(self.substances.to_owned()).reaction(),
            Some(file) => chemdraw::parse_chemdraw_file(file),
        }?;
        rxn.balance()
    }
}
