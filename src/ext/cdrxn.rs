use crate::model::Substance;
use serde::Deserialize;

#[derive(Deserialize, Debug)]
struct Reaction {
    #[serde(rename(deserialize = "STEPS"))]
    steps: Vec<Step>,
}

impl Reaction {
    pub fn reactants(&self) -> Vec<String> {
        self.steps
            .first()
            .map(|s| s.reactants.iter().chain(s.reagents.iter()).collect())
            .unwrap_or(Vec::new())
            .iter()
            .map(|m| m.formula())
            .collect()
    }

    pub fn products(&self) -> Vec<String> {
        self.steps
            .first()
            .map(|s| &s.products)
            .unwrap_or(&Vec::new())
            .iter()
            .map(|m| m.formula())
            .collect()
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

#[derive(Deserialize, Debug)]
struct Molecule {
    #[serde(rename(deserialize = "NAME"))]
    name: String,
    #[serde(rename(deserialize = "FORMULA"))]
    raw_formula: String,
    #[serde(rename(deserialize = "SMILES"))]
    smiles: String,
}

impl Molecule {
    pub fn formula(&self) -> String {
        self.raw_formula.replace("<sub>", "").replace("</sub>", "")
    }
}

#[derive(Debug)]
pub struct ParsedReaction {
    pub reactants: Vec<Substance>,
    pub products: Vec<Substance>,
}

pub fn parse_chemdraw_reaction(
    document: &str,
) -> Result<ParsedReaction, String> {
    let parsed: Vec<Reaction> = serde_json::from_str(document)
        .map_err(|e| format!("Could not parse: {:?}", e))?;
    let rxn = parsed.first().ok_or("No reactions!".to_string())?;
    let reactants: Vec<Substance> = rxn
        .reactants()
        .iter()
        .map(|r| Substance::new(r))
        .collect::<Result<Vec<Substance>, String>>()?;
    let products: Vec<Substance> = rxn
        .products()
        .iter()
        .map(|r| Substance::new(r))
        .collect::<Result<Vec<Substance>, String>>()?;
    Ok(ParsedReaction {
        reactants,
        products,
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
