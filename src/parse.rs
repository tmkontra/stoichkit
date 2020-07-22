use std::collections::hash_map::RandomState;
use std::collections::HashMap;
use std::panic;

use ptable::Element;

fn get_element(symbol: &str) -> Result<Element, String> {
    match symbol.chars().all(|c| c.is_ascii_alphabetic()) {
        true => {
            let e: Result<Option<Element>, _> =
                panic::catch_unwind(|| Element::from_symbol(symbol));
            e.unwrap_or(None)
                .ok_or(format!("Invalid symbol {}", symbol))
        }
        false => Err(format!("Invalid symbol {}", symbol)),
    }
}

// translated from https://leetcode.com/articles/number-of-atoms/#
pub fn parse_formula(formula: &str) -> Result<HashMap<Element, u32, RandomState>, String> {
    let mut stack: Vec<HashMap<Element, u32>> = vec![HashMap::new()];
    let mut i: usize = 0;
    let N: usize = formula.len();
    let mut broken: bool = false;
    trace!("Parsing formula {:?}", formula);
    while i < N && !broken {
        match formula.chars().nth(i).unwrap() {
            '(' | '[' | '{' => {
                trace!("Start of group at {:?}", i);
                stack.push(HashMap::new());
                i += 1;
            }
            ')' | ']' | '}' => {
                match stack.pop() {
                    Some(top) => {
                        i += 1;
                        let i_start = i;
                        while i < N && formula.chars().nth(i).unwrap().is_digit(10) {
                            i += 1;
                        }
                        let mult: u32 = match i_start == i {
                            true => 1,
                            false => {
                                let multiplicity_str = formula.get(i_start..i).unwrap();
                                multiplicity_str.parse::<u32>().unwrap()
                            }
                        };
                        trace!("Got multiplicity {:?} for group {:?}", mult, top);
                        for (elem, v) in top {
                            match stack.last() {
                                Some(last) => {
                                    let curr = match last.get(&elem) {
                                        Some(val) => val.clone(),
                                        None => {
                                            stack.last_mut().unwrap().insert(elem, 0).unwrap_or(0)
                                        }
                                    };
                                    let addt = v.clone() * mult;
                                    let new = curr.clone() + addt;
                                    stack.last_mut().unwrap().insert(elem, new);
                                }
                                None => broken = true,
                            }
                        }
                    }
                    None => broken = true,
                };
            }
            token if token.is_ascii_alphanumeric() => {
                let mut i_start = i;
                i += 1;
                while i < N && formula.chars().nth(i).unwrap().is_lowercase() {
                    i += 1;
                }
                let name = formula.get(i_start..i).unwrap();
                i_start = i;
                trace!("Captured symbol {:?}", name);
                let elem: Element = get_element(name)?;
                while i < N && formula.chars().nth(i).unwrap().is_digit(10) {
                    i += 1
                }
                let mult: u32 = match i_start == i {
                    true => 1,
                    false => {
                        let multiplicity_str = formula.get(i_start..i).unwrap();
                        multiplicity_str.parse::<u32>().unwrap()
                    }
                };
                trace!("Got multiplicity {:?} for element {:?}", mult, name);
                match stack.last() {
                    Some(last) => {
                        let curr = match last.get(&elem) {
                            Some(val) => val.clone(),
                            None => stack.last_mut().unwrap().insert(elem, 0).unwrap_or(0),
                        };
                        let new = curr.to_owned() + mult;
                        stack.last_mut().unwrap().insert(elem, new);
                    }
                    None => broken = true,
                }
            }
            invalid => {
                error!("Got invalid character {:?} in formula {}", invalid, formula);
                broken = true;
            }
        }
    }

    if broken {
        Err(String::from("Could not parse"))
    } else {
        let result = stack
            .last()
            .ok_or(format!("Error parsing {}", formula))?
            .clone()
            .to_owned();
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use crate::parse::parse_formula;
    use crate::test_utils::e;
    use ptable::Element;
    use std::collections::hash_map::RandomState;
    use std::collections::HashMap;

    #[test]
    fn ethane() {
        let formula = "C2H6";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("C", 2), ("H", 6)].iter().cloned().collect();
        assert_eq!(result, e(expected));
    }

    #[test]
    fn nickel_tert_butoxide() {
        let formula = "Ni[OC(CH3)3]2";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("Ni", 1), ("O", 2), ("C", 8), ("H", 18)]
            .iter()
            .cloned()
            .collect();
        assert_eq!(result, e(expected));
    }

    #[test]
    fn missing_open_bracket() {
        let formula = "(C2H3)3)";
        let result = parse_formula(formula);
        assert!(result.is_err(), "{:?}", result.ok());
    }

    #[test]
    #[should_panic] // todo: better enforce matching brakcets
    fn missing_closing_bracket() {
        let formula = "((C2H3)3";
        let result = parse_formula(formula);
        assert!(result.is_err(), "{:?}", result.ok());
    }

    #[test]
    fn mismatched_bracket_types() {
        let formula: &str = "{C2H6)12";
        let result: HashMap<Element, u32, RandomState> = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("C", 24), ("H", 72)].iter().cloned().collect();
        assert_eq!(result, e(expected));
    }

    #[test]
    fn invalid_symbol() {
        let formula: &str = "CkCoNv30";
        let result = parse_formula(formula);
        assert!(result.is_err());
    }

    #[test]
    fn invalid_string() {
        let formula: &str = "30k";
        let result = parse_formula(formula);
        assert!(result.is_err());
    }
}
