mod v2;

use std::collections::hash_map::RandomState;
use std::collections::HashMap;
use std::fmt::format;
use std::panic;

use crate::model::Element;

// translated from https://leetcode.com/articles/number-of-atoms/#
pub fn parse_formula(
    formula: &str,
) -> Result<HashMap<Element, usize, RandomState>, String> {
    let mut stack: Vec<HashMap<Element, usize>> = vec![HashMap::new()];
    let mut i: usize = 0;
    let formula_len: usize = formula.len();
    let mut broken: bool = false;
    trace!("Parsing formula {:?}", formula);
    while i < formula_len && !broken {
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
                        while i < formula_len
                            && formula.chars().nth(i).unwrap().is_digit(10)
                        {
                            i += 1;
                        }
                        let mult: usize = match i_start == i {
                            true => 1,
                            false => {
                                let multiplicity_str =
                                    formula.get(i_start..i).unwrap();
                                multiplicity_str.parse::<usize>().unwrap()
                            }
                        };
                        trace!(
                            "Got multiplicity {:?} for group {:?}",
                            mult,
                            top
                        );
                        for (elem, v) in top {
                            match stack.last() {
                                Some(last) => {
                                    let curr = match last.get(&elem) {
                                        Some(val) => *val,
                                        None => stack
                                            .last_mut()
                                            .unwrap()
                                            .insert(elem, 0)
                                            .unwrap_or(0),
                                    };
                                    let addt = v * mult;
                                    let new = curr + addt;
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
                while i < formula_len
                    && formula.chars().nth(i).unwrap().is_lowercase()
                {
                    i += 1;
                }
                let name = formula.get(i_start..i).unwrap();
                i_start = i;
                trace!("Captured symbol {:?}", name);
                let elem: Element = element_from_string(name)?;
                while i < formula_len
                    && formula.chars().nth(i).unwrap().is_digit(10)
                {
                    i += 1
                }
                let mult: usize = match i_start == i {
                    true => 1,
                    false => {
                        let multiplicity_str = formula.get(i_start..i).unwrap();
                        multiplicity_str.parse::<usize>().unwrap()
                    }
                };
                trace!("Got multiplicity {:?} for element {:?}", mult, name);
                match stack.last() {
                    Some(last) => {
                        let curr = match last.get(&elem) {
                            Some(val) => *val,
                            None => stack
                                .last_mut()
                                .unwrap()
                                .insert(elem, 0)
                                .unwrap_or(0),
                        };
                        let new = curr.to_owned() + mult;
                        stack.last_mut().unwrap().insert(elem, new);
                    }
                    None => broken = true,
                }
            }
            invalid => {
                error!(
                    "Got invalid character {:?} in formula {}",
                    invalid, formula
                );
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
            .to_owned();
        Ok(result)
    }
}

fn element_from_string(symbol: &str) -> Result<Element, String> {
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

pub fn get_element_by_id(id: u8) -> Result<Element, String> {
    let e: Result<Option<Element>, _> =
        panic::catch_unwind(|| Element::from_atomic_number(id.into()));
    e.unwrap_or(None)
        .ok_or(format!("Invalid atomic number {}", id))
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::RandomState;
    use std::collections::HashMap;

    use crate::model::Element;
    use crate::parse::parse_formula;
    use crate::test_utils::parse_elements;

    #[test]
    fn ethane() {
        let formula = "C2H6";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, usize> =
            [("C", 2), ("H", 6)].iter().cloned().collect();
        assert_eq!(result, parse_elements(expected));
    }

    #[test]
    fn nickel_tert_butoxide() {
        let formula = "Ni[OC(CH3)3]2";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, usize> =
            [("Ni", 1), ("O", 2), ("C", 8), ("H", 18)]
                .iter()
                .cloned()
                .collect();
        assert_eq!(result, parse_elements(expected));
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
        let result: HashMap<Element, usize, RandomState> =
            parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, usize> =
            [("C", 24), ("H", 72)].iter().cloned().collect();
        assert_eq!(result, parse_elements(expected));
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
