StoichKit
===

A toolkit for stoichiometry.

### Features
- Calculates percent yield
  - given limiting reagent, product and respective masses (in grams)

### Usage

```$xslt
$ stoichkit --help
StoichKit 0.1

USAGE:
    stoichkit <reagent_formula> <reagent_grams> <product_formula> <product_grams>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

ARGS:
    <reagent_formula>    
    <reagent_grams>      
    <product_formula>    
    <product_grams>

$ stoichkit C2H6 30.00213 "(CH3)2" 26.98
Yield: 1.1120138
```

### Installation

1. Clone this repo
2. With [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) installed, run `cargo build --release`
3. Run `./target/release/stoichkit` or copy that binary to a bin folder.


### Roadmap

1. Accept all reagents and determine limiting reagent.
2. Add desktop GUI
3. Use StoichKit to power a web UI.
 
