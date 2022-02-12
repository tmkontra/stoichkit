StoichKit
===

A toolkit for stoichiometry.

### Features
- `balance`: Balances a chemical equation
- `yield`: Calculates percent yield
  - given a fully balanced chemical reaction, and respective masses (in grams)
- `moles` calculates moles given formula and mass (grams)
- `theoretical-yield`: Calculate theoretical yield of all products 
  - given fully balanced chemical equation and reagent masses

### Usage

```$xslt
stoichkit 0.6.0
A stoichiometry toolkit.

USAGE:
    stoichkit <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    balance              
    help                 Print this message or the help of the given subcommand(s)
    moles                
    theoretical-yield    
    yield
```

### Examples

#### Balance
```
$ stoichkit balance H2O O2 = H2O2
2*H2O + O2 = 2*H2O2

$ stoichkit balance -x H2O O2 = H2O2
2*H2O + 1*O2 = 2*H2O2
```

#### Moles
```
$ stoichkit moles C4H6 0.7254
0.013410485 mol
```

#### Theoretical Yield
```
$ stoichkit theoretical-yield "2*H2O2" 4.0 = "2*H2O" O2 
H2O 2.1185393 g
O2 1.8814605 g

$ stoichkit theoretical-yield --units moles "2*H2O2" 4.0 = "2*H2O" O2
H2O 0.11759864 mol
O2 0.05879932 mol
```

#### Yield

```
$ stoichkit yield '2*H2O2' 4.0 = '2*H2O' 2.1184621`
H2O 0.9999635 %
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
- [ ] Ideal Gas calculations
 
