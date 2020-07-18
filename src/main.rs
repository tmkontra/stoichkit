mod parse;
mod molecule;

use crate::parse::parse_formula;
use crate::molecule::molecular_weight;

fn main() {
    let formula = std::env::args().nth(1).expect("Formula Required");
    let grams_arg = std::env::args()
        .nth(2)
        .ok_or("No grams")
        .and_then(|g| g.parse::<f32>().map_err(|_| "No grams"));
    let parsed = parse_formula(formula.as_str());
    match parsed
        .and_then(|atoms| molecular_weight(atoms)) {
        Ok(molwt) => {
            println!("MolWt: {:?}", molwt);
            match grams_arg {
                Ok(grams) => println!("Moles: {:?}", grams / molwt),
                Err(_) => ()
            }
        },
        Err(msg) => println!("ERROR: {:?}", msg)
    }
}
