StoichKit
===

A toolkit for stoichiometry.

### Features
- `balance`: Balances a chemical equation
- `yield`: Calculates percent yield
  - given a fully balanced chemical reaction, and respective masses (in grams)

### Usage

```$xslt
$ stoichkit --help
stoichkit 0.2.0

USAGE:
    stoichkit <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    balance    
    help       Prints this message or the help of the given subcommand(s)
    yield 
```

```
$ stoichkit yield "2*Al" 2.8 "3*Cl2" 4.25 "2*AlCl3" 4.889
Yield: 0.91756254
```
```
$ stoichkit balance Al Cl2 = AlCl3                     
2 Al + 3 Cl2 = 2 AlCl3
```

### Installation

1. Clone this repo
2. With [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) installed, run `cargo build --release`
3. Run `./target/release/stoichkit` or copy that binary to a bin folder.

`stoichkit` equation balancer uses the `nalgebra-linalg` solver, which requires a BLAS installation.

On macOS BLAS can be installed via `brew install openblas`.


### Roadmap

- [x] Accept all reagents and determine limiting reagent.
- [ ] Add desktop GUI (maybe)
- [x] Use StoichKit to power a web UI (`stoichkitweb` available [here](https://github.com/ttymck/stoichkitweb))
- [x] Implement chemical equation balancer
- [ ] Test linear algebra suite on other platforms
 
