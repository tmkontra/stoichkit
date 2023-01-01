use std::fmt::Error;
use clap::{Args, Parser, Subcommand};

use crate::ext::chemdraw;
use crate::model::ReactionList;
use crate::model::{Compound, Units, YieldUnits};

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
            Commands::TheoreticalYield(
                TheoreticalYieldArgs { reaction_list, units}
            ) => Cli::run_theoretical_yield_command(reaction_list, units),
            Commands::Yield(
                YieldArgs { reaction_list }
            ) => Cli::run_yield_command(reaction_list),
            Commands::Balance(BalanceEquationArgs{ reaction_list, chemdraw_file, explicit }) =>
                Cli::run_balance_command(reaction_list, chemdraw_file, explicit),
            Commands::Moles(MolesArgs { reaction_list }) =>
                Cli::run_moles_command(reaction_list),
        };
        match result {
            Ok(_) => (),
            Err(err) => println!("ERROR: {}", err),
        }
    }

    fn run_theoretical_yield_command(
        reaction_list: ReactionList,
        units: Option<YieldUnits>,
    ) -> Result<(), String> {
        let units = units.unwrap_or(YieldUnits::Mass);
        reaction_list
            .parse_theoretical_reaction()
            .map(|r| {
                print_result_list(
                    r.yields(&units)
                        .iter()
                        .map(|(r, amt)| (&r.compound, *amt))
                        .collect(),
                    units.into(),
                )
            })
    }

    fn run_yield_command(reaction_list: ReactionList) -> Result<(), String> {
        reaction_list.parse_yield_reaction().map(|yld| {
            print_result_list(
                vec![(&yld.product.reactant.compound, yld.percent_yield())],
                Units::Percent,
            )
        })
    }

    fn run_balance_command(reaction_list: ReactionList,
                           chemdraw_file: Option<String>,
                           explicit: bool) -> Result<(), String> {
        let rxn = match chemdraw_file {
            Some(file) => chemdraw::parse_chemdraw_file(file.as_ref()),
            None => reaction_list.parse_reaction(),
        }?;
        let balanced_rxn = rxn.balance();
        balanced_rxn.map(|balanced| {
            println!("{}", balanced.display_string(explicit))
        })
    }

    fn run_moles_command(reaction_list: ReactionList) -> Result<(), String> {
        reaction_list.to_samples().map(|subs| {
            subs.iter().for_each(|s| {
                println!("{} mol", s.mass / s.reactant.compound.molar_mass)
            })
        })
    }
}

#[derive(Subcommand)]
enum Commands {
    TheoreticalYield(TheoreticalYieldArgs),
    Yield(YieldArgs),
    Balance(BalanceEquationArgs),
    Moles(MolesArgs), // Pvnrt(GasArgs) TODO,
}

#[derive(Args)]
struct TheoreticalYieldArgs {
    #[clap(parse(try_from_str = parse_reaction_list))]
    reaction_list: ReactionList,
    #[clap(short, long, arg_enum)]
    units: Option<YieldUnits>,
}

fn parse_reaction_list(arg: &str) -> Result<ReactionList, Error> {
    let args = arg.split(" ").map(String::from).collect();
    Ok(ReactionList::new(args))
}

#[derive(Args)]
struct YieldArgs {
    #[clap(parse(try_from_str = parse_reaction_list))]
    reaction_list: ReactionList,
}

#[derive(Args)]
struct BalanceEquationArgs {
    #[clap(help = "Chemical equation [...reactants] = [...products]")]
    #[clap(parse(try_from_str = parse_reaction_list))]
    reaction_list: ReactionList,
    #[clap(short, conflicts_with = "substances")]
    chemdraw_file: Option<String>,
    #[clap(short = 'x', long)]
    explicit: bool,
}

#[derive(Args)]
struct MolesArgs {
    #[clap(parse(try_from_str = parse_reaction_list))]
    reaction_list: ReactionList
}

fn print_result_list(list: Vec<(&Compound, f32)>, units: Units) {
    list.iter().for_each(|(product, yld)| {
        println!("{} {} {}", product.formula, yld, units)
    })
}
