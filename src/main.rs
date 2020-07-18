use std::collections::HashMap;
use std::collections::hash_map::RandomState;

// translated from https://leetcode.com/articles/number-of-atoms/#
fn parse_formula(formula: &str) -> HashMap<String, u32, RandomState> {
    let mut stack: Vec<HashMap<String, u32>> = vec![HashMap::new()];
    let mut i = 0;
    let N = formula.len();
    while i < N {
        println!("{:?}", formula.chars().nth(i));
        println!("{:?}", formula.chars());
        match formula.chars().nth(i).unwrap() {
            '(' => {
                stack.push(HashMap::new());
                i += 1;
            },
            ')' => {
                let top = stack.pop().unwrap();
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
                    let curr = match stack
                        .last()
                        .unwrap()
                        .get(name.as_str()) {
                        Some(val) => val.clone(),
                        None => stack.last_mut().unwrap().insert(name.to_owned(), 0).unwrap_or(0)
                    };
                    let addt = v.clone() * mult;
                    let new = curr.clone() + addt;
                    stack.last_mut().unwrap().insert(name.to_owned(), new);
                }
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
                println!("digit so name {:?} and i_start {:?} for mult", name, i_start);
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
                let curr = match stack
                    .last()
                    .unwrap()
                    .get(name) {
                    Some(val) => val.clone(),
                    None => stack.last_mut().unwrap().insert(name.to_owned(), 0).unwrap_or(0)
                };
                let new = curr.to_owned() + mult;
                stack.last_mut().unwrap().insert(name.to_string(), new);
            }
        }
    }
    return stack.last().unwrap().clone();
}


fn main() {
    let formula = std::env::args().nth(1).expect("Formula Required");
    let result = parse_formula(formula.as_str());
    println!("{:?}", result)
}
