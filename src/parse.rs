use std::collections::hash_map::RandomState;
use std::collections::HashMap;

// translated from https://leetcode.com/articles/number-of-atoms/#
pub fn parse_formula(formula: &str) -> Result<HashMap<&str, u32, RandomState>, String> {
    let mut stack: Vec<HashMap<&str, u32>> = vec![HashMap::new()];
    let mut i: usize = 0;
    let N: usize = formula.len();
    let mut broken: bool = false;
    while i < N && !broken {
        println!("{:?}", formula.chars().nth(i));
        println!("{:?}", formula.chars());
        match formula.chars().nth(i).unwrap() {
            '(' | '[' | '{' => {
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
                        for (name, v) in top {
                            match stack.last() {
                                Some(last) => {
                                    let curr = match last.get(name) {
                                        Some(val) => val.clone(),
                                        None => {
                                            stack.last_mut().unwrap().insert(name, 0).unwrap_or(0)
                                        }
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
            _ => {
                let mut i_start = i;
                i += 1;
                println!("character so i_start {:?} and i {:?}", i_start, i);
                while i < N && formula.chars().nth(i).unwrap().is_lowercase() {
                    println!("passing {:?} at {:?}", formula.chars().nth(i), i);
                    i += 1;
                }
                let name = formula.get(i_start..i).unwrap();
                i_start = i;
                println!(
                    "digit so name {:?} and i_start {:?} for mult",
                    name, i_start
                );
                while i < N && formula.chars().nth(i).unwrap().is_digit(10) {
                    println!("{:?} at {:?}", formula.chars().nth(i), i);
                    i += 1
                }
                let mult: u32 = match i_start == i {
                    true => 1,
                    false => {
                        let multiplicity_str = formula.get(i_start..i).unwrap();
                        multiplicity_str.parse::<u32>().unwrap()
                    }
                };
                match stack.last() {
                    Some(last) => {
                        let curr = match last.get(name) {
                            Some(val) => val.clone(),
                            None => stack.last_mut().unwrap().insert(name, 0).unwrap_or(0),
                        };
                        let new = curr.to_owned() + mult;
                        stack.last_mut().unwrap().insert(name, new);
                    }
                    None => broken = true,
                }
            }
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
    use std::collections::HashMap;

    #[test]
    fn ethane() {
        let formula = "C2H6";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("C", 2), ("H", 6)].iter().cloned().collect();
        assert_eq!(result, expected);
    }

    #[test]
    fn nickel_tert_butoxide() {
        let formula = "Ni[OC(CH3)3]2";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("Ni", 1), ("O", 2), ("C", 8), ("H", 18)]
            .iter()
            .cloned()
            .collect();
        assert_eq!(result, expected);
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
        let formula = "{C2H6)12";
        let result = parse_formula(formula).ok().unwrap();
        let expected: HashMap<&str, u32> = [("C", 24), ("H", 72)].iter().cloned().collect();
        assert_eq!(result, expected);
    }
}
