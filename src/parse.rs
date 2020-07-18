use std::collections::hash_map::RandomState;
use std::collections::HashMap;

// translated from https://leetcode.com/articles/number-of-atoms/#
pub fn parse_formula(formula: &str) -> Result<HashMap<String, u32, RandomState>, String> {
    let mut stack: Vec<HashMap<String, u32>> = vec![HashMap::new()];
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
                        for (name, v) in top {
                            match stack.last() {
                                Some(last) => {
                                    let curr = match last.get(name.as_str()) {
                                        Some(val) => val.clone(),
                                        None => stack
                                            .last_mut()
                                            .unwrap()
                                            .insert(name.to_string(), 0)
                                            .unwrap_or(0),
                                    };
                                    let addt = v.clone() * mult;
                                    let new = curr.clone() + addt;
                                    stack.last_mut().unwrap().insert(name, new);
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
                trace!("Found element {:?}", name);
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
                        let curr = match last.get(name) {
                            Some(val) => val.clone(),
                            None => stack
                                .last_mut()
                                .unwrap()
                                .insert(name.to_string(), 0)
                                .unwrap_or(0),
                        };
                        let new = curr.to_owned() + mult;
                        stack.last_mut().unwrap().insert(name.to_string(), new);
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
    let result = stack
        .last()
        .ok_or(format!("Error parsing {}", formula))?
        .clone();
    for (symbol, _) in result {
        match symbol.chars().all(|c| c.is_ascii_alphabetic()) {
            true => continue,
            false => broken = true,
        }
    }
    if broken {
        Err(String::from("Could not parse"))
    } else {
        Ok(stack.last().unwrap().clone())
    }
}

#[cfg(test)]
mod tests {
    use crate::parse::parse_formula;
    use std::collections::hash_map::RandomState;
    use std::collections::HashMap;

    use crate::test_utils::e;

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
    fn missing_close_bracket() {
        let formula = "((C2H3)3";
        let result = parse_formula(formula);
        assert!(result.is_err(), "{:?}", result.ok());
    }

    #[test]
    fn mismatched_bracket_types() {
        let formula: &str = "{C2H6)12";
        let result: HashMap<String, u32, RandomState> = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("C", 24), ("H", 72)].iter().cloned().collect();
        assert_eq!(result, e(expected));
    }
}
