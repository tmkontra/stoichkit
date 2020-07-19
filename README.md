StoichKit
===

A toolkit for stoichiometry.

### Features
- Calculates percent yield
  - given a fully balanced chemical reaction, and respective masses (in grams)

### Usage

```$xslt
$ stoichkit --help
stoichkit 0.1.0

USAGE:
    stoichkit [substances]...

ARGS:
    <substances>...    Fully balanced chemical reaction list: [<coeff>*]<formula> <grams>, coeff defaults to 1

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
```

```
$ stoichkit "2*Al" 2.8 "3*Cl2" 4.25 "2*AlCl3" 4.889
Yield: 0.91756254
```

### Installation

1. Clone this repo
2. With [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) installed, run `cargo build --release`
3. Run `./target/release/stoichkit` or copy that binary to a bin folder.


### Roadmap

- [x] Accept all reagents and determine limiting reagent.
- [ ] Add desktop GUI
- [ ] Use StoichKit to power a web UI.
 
