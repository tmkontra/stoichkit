use crate::model::{Substance};
use std::collections::HashSet;
use std::iter;
use counter::Counter;
use ptable::Element;
use std::convert::TryInto;
use rug::Rational;

use nalgebra::{DMatrix, Matrix, Vector, DVector};
use nalgebra::linalg::*;
use num::integer::lcm;
use std::cmp::max;
use itertools::Itertools;


fn balance(reagents: Vec<Substance>, products: Vec<Substance>) {
    let mut all_atoms: HashSet<&Element> = HashSet::new();
    for r in &reagents {
        for e in r.atoms.keys() {
            all_atoms.insert(e);
        };
    }
    for p in &products {
        for e in p.atoms.keys() {
            all_atoms.insert(e);
        }
    }
    let mut elements: Vec<&Element> = Vec::new();
    let mut matrix: Vec<f64> = Vec::new();
    for (i, element) in all_atoms.iter().enumerate() {
        let mut rge: Vec<f64> = reagents
            .iter()
            .map(|s| s.atoms.get(element).map(|c| *c).unwrap_or(0 as u32).clone())
            .map(|c| c as f64)
            .collect();
        let pde: Vec<f64> = products
            .iter()
            .map(|s| s.atoms.get(element).map(|c| *c).unwrap_or(0 as u32).clone())
            .map(|c| c as f64 * -1 as f64)
            .collect();
        rge.extend(&pde);
        elements.append(&mut vec![element.to_owned()]);
        matrix.append(&mut rge);
    }
    println!("{:?}", matrix);
    println!("{:?}", matrix.len());
    println!("{:?}", (elements.len()));
    println!("{:?}", matrix.len() - (elements.len()));
    let mx = DMatrix::from_row_slice(elements.len(), reagents.len() + products.len(), matrix.as_slice());
    println!("{:?}", mx);
    let c = reagents.len()+products.len()-1;
    let (a, b) = mx.columns_range_pair(0..c, c..);
    println!("{:?}", b);
    println!("{:?}", a);
    let x = a.svd(true, true);
    println!("{:?}", x);
    let solve = x.solve(&b, 0.0).unwrap();
    println!("{:?}", solve);
    let c: Vec<f32> = solve.column(0).iter().map(|c| *c as f32).collect();
    println!("{:?}", c);
    let fc: Vec<Rational> = c.iter().map(|c| Rational::from_f32(*c).unwrap()).collect();
    println!("{:?}", fc);
    let denoms: Vec<_> = fc.iter().map(|c| c.denom()).collect();
    println!("{:?}", denoms);
    let mult = denoms
        .iter()
        .combinations(2)
        .fold(1 as i32, |mult, cur| {
            println!("{:?}", cur);
            max(mult, lcm(cur[0].to_i32().unwrap(), cur[1].to_i32().unwrap()))
        });
    println!("LCM: {:?}", mult);
    let fin: Vec<_> = fc.iter().map(|f| f * Rational::from((mult, 1))).map(|f| f.numer().to_owned()).collect();
    println!("{:?}", fin);
}

#[cfg(test)]
mod tests {
    use crate::parse::parse_formula;
    use crate::test_utils::e;
    use ptable::Element;
    use std::collections::hash_map::RandomState;
    use std::collections::HashMap;
    use crate::model::*;
    use crate::solve::balance;

    #[test]
    fn test() {
        let rg = vec![Substance::new("Al", 3.0, None).unwrap(),
                      Substance::new("Cl2", 3.0, None).unwrap()];
        let pd = vec![Substance::new("AlCl3", 3.0, None).unwrap()];
        let result = balance(rg, pd);
        assert_eq!(result, ());
    }

    #[test]
    fn test_2() {
        let rg = vec![Substance::new("C2H2Cl4", 3.0, None).unwrap(),
                      Substance::new("Ca(OH)2", 3.0, None).unwrap()];
        let pd = vec![Substance::new("C2H2Cl3", 3.0, None).unwrap(),
                      Substance::new("CaCl2", 3.0, None).unwrap(),
                      Substance::new("H2O", 3.0, None).unwrap()];
        let result = balance(rg, pd);
        assert_eq!(result, ());
    }
}

