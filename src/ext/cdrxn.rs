use crate::model::{Substance, UnbalancedReaction};
use serde::Deserialize;

#[derive(Deserialize, Debug)]
struct ChemdrawReaction {
    #[serde(rename(deserialize = "STEPS"))]
    steps: Vec<Step>,
}

impl ChemdrawReaction {
    pub fn reactants(&self) -> Result<Vec<Substance>, String> {
        self.get_molecules(|s| {
            s.reactants
                .iter()
                .chain(s.reagents.iter())
                .cloned()
                .collect()
        })
    }

    pub fn products(&self) -> Result<Vec<Substance>, String> {
        self.get_molecules(|s| s.products.clone())
    }

    fn get_molecules<F: Fn(&Step) -> Vec<Molecule>>(
        &self,
        f: F,
    ) -> Result<Vec<Substance>, String> {
        self.steps
            .first()
            .map(f)
            .unwrap_or(Vec::new())
            .iter()
            .map(|m: &Molecule| m.to_substance())
            .collect::<Result<Vec<Substance>, String>>()
    }
}

#[derive(Deserialize, Debug)]
struct Step {
    #[serde(rename(deserialize = "REACTANTS"))]
    reactants: Vec<Molecule>,
    #[serde(rename(deserialize = "REAGENTS"))]
    reagents: Vec<Molecule>,
    #[serde(rename(deserialize = "PRODUCTS"))]
    products: Vec<Molecule>,
}

#[derive(Clone, Deserialize, Debug)]
struct Molecule {
    #[serde(rename(deserialize = "NAME"))]
    name: String,
    #[serde(rename(deserialize = "FORMULA"))]
    raw_formula: String,
    #[serde(rename(deserialize = "SMILES"))]
    smiles: String,
}

impl Molecule {
    pub fn to_substance(&self) -> Result<Substance, String> {
        let formula =
            self.raw_formula.replace("<sub>", "").replace("</sub>", "");
        Substance::new(formula.as_str())
    }
}

pub fn parse_chemdraw_reaction(
    document: &str,
) -> Result<UnbalancedReaction, String> {
    let parsed: Vec<ChemdrawReaction> = serde_json::from_str(document)
        .map_err(|e| format!("Could not parse: {:?}", e))?;
    let rxn = parsed.first().ok_or("No reactions!".to_string())?;
    Ok(UnbalancedReaction {
        reactants: rxn.reactants()?,
        products: rxn.products()?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        let document = r#"[{
"STEPS":[{
"REACTANTS":[{"ID":"54","NAME":"dimethylamine","FORMULA":"C<sub>2</sub>H<sub>7</sub>N","MW":"45.0850000005 g/mol","EM":45.0578492299,"INCHI":"InChI=1S/C2H7N/c1-3-2/h3H,1-2H3","SMILES":"CNC"}],
"REAGENTS":[],
"PRODUCTS":[{"ID":"59","NAME":"methane","FORMULA":"CH<sub>4</sub>","MW":"16.04300000025 g/mol","EM":16.031300128399998,"INCHI":"InChI=1S/CH4/h1H4","SMILES":"C"}]}]}]
"#;
        let result = parse_chemdraw_reaction(document).unwrap();
        assert_eq!(result.reactants.len(), 1);
        assert_eq!(result.products.len(), 1);
        assert_eq!(result.products.first().unwrap().formula, "CH4");
    }
}
