extern crate clap;
extern crate core;
#[macro_use]
extern crate log;

use crate::clap::Parser;
use stoichkit::cli::Cli;

fn main() {
    env_logger::Builder::from_env("STOICHKIT_LOG").init();
    let opts: Cli = Cli::parse();
    opts.run();
}
