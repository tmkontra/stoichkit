use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::iter;

use counter::Counter;
use itertools::Itertools;
use nalgebra::linalg::*;
use nalgebra::{DMatrix, DVector, Matrix, Vector};
use ndarray::prelude::*;
use ndarray_linalg::Solve;
use num::abs;
use num::integer::lcm;
use ptable::Element;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use rug::Rational;

use crate::model::Substance;

pub fn balance(reagents: Vec<Substance>, products: Vec<Substance>) -> Vec<(String, i32)> {
    let mut all_atoms: HashSet<&Element> = HashSet::new();
    for r in &reagents {
        for e in r.atoms.keys() {
            all_atoms.insert(e);
        }
    }
    for p in &products {
        for e in p.atoms.keys() {
            all_atoms.insert(e);
        }
    }
    let mut elements: Vec<&Element> = all_atoms.iter().map(|c| c.clone()).collect();
    let mut subs: Vec<&Substance> = Vec::new();
    let mut matrix: Vec<f64> = Vec::new();
    for (i, element) in all_atoms.iter().enumerate() {
        for rge in &reagents {
            subs.push(rge);
            let rxe = rge
                .atoms
                .get(element)
                .map(|c| *c)
                .unwrap_or(0 as u32)
                .clone();
            matrix.push(rxe as f64);
        }
        for pde in &products {
            subs.push(pde);
            let rxe = pde
                .atoms
                .get(element)
                .map(|c| *c)
                .unwrap_or(0 as u32)
                .clone();
            matrix.push(rxe as f64);
        }
    }
    let mx = DMatrix::from_row_slice(
        elements.len(),
        reagents.len() + products.len(),
        matrix.as_slice(),
    );
    let c = reagents.len() + products.len() - 1;
    let (a, b) = mx.columns_range_pair(0..c, c..);
    let x = a.svd(true, true);
    let solve = x.solve(&b, 0.0).unwrap();
    let c: Vec<f32> = solve.column(0).iter().map(|c| *c as f32).collect();
    let fc: Vec<Rational> = c
        .iter()
        .map(|c| limit_denominator(&Rational::from_f32(*c).unwrap().abs(), 100))
        .collect();
    let denoms: Vec<_> = fc.iter().map(|c| c.denom()).collect();
    let mult = denoms.iter().combinations(2).fold(1 as i32, |mult, cur| {
        max(
            mult,
            lcm(cur[0].to_i32().unwrap(), cur[1].to_i32().unwrap()),
        )
    });
    let mut fin: Vec<_> = fc
        .iter()
        .map(|f| f * Rational::from((mult, 1)))
        .map(|f| f.numer().to_i32().unwrap())
        .collect();
    fin.push(mult);
    let result: Vec<(String, i32)> = subs
        .iter()
        .map(|s| s.to_owned().formula.clone())
        .zip(&mut fin.iter().map(|c| c.to_owned()))
        .collect();
    result
}

fn limit_denominator(f: &Rational, max_denom: u32) -> Rational {
    debug!("Limiting denom for {:?} to at most {:?}", f, max_denom);
    if max_denom < 1 || f.denom() <= &max_denom {
        f.clone()
    } else {
        let (mut p0, mut q0, mut p1, mut q1) = (0, 1, 1, 0);
        let (mut n, mut d) = (f.numer().to_i32().unwrap(), f.denom().to_i32().unwrap());
        let mut a: i32;
        let mut q2: i32;
        loop {
            a = n / d;
            q2 = q0 + (a * q1);
            if q2 > max_denom as i32 {
                break;
            }
            let p0new: i32 = p1.clone();
            let q0new: i32 = q1.clone();
            let p1new: i32 = (p0.clone() + (a.clone() * p1.clone()));
            let q1new: i32 = q2.clone();
            p0 = p0new;
            q0 = q0new;
            p1 = p1new;
            q1 = q1new;
            let nnew: i32 = d.clone();
            let dnew: i32 = (n.clone() - (a.clone() * d.clone()));
            n = nnew;
            d = dnew;
        }
        let k = (max_denom as i32 - q0) / q1;
        let bound1: Rational = Rational::from((p0 + k * p1, q0 + k * q1));
        let bound2: Rational = Rational::from((p1, q1));
        let d1f: Rational = (bound1.clone() - f);
        let d2f: Rational = bound2.clone() - f;
        return if d1f.abs() <= d2f.abs() {
            bound1.clone()
        } else {
            bound2.clone()
        };
    }
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::RandomState;
    use std::collections::HashMap;

    use ptable::Element;

    use crate::model::*;
    use crate::parse::parse_formula;
    use crate::solve::balance;
    use crate::test_utils::e;

    #[test]
    fn test() {
        let rg = vec![
            Substance::new("Al", 3.0, None).unwrap(),
            Substance::new("Cl2", 3.0, None).unwrap(),
        ];
        let pd = vec![Substance::new("AlCl3", 3.0, None).unwrap()];
        let bal = balance(rg, pd);
        let result: Vec<(&str, i32)> = bal
            .iter()
            .map(|(s, c)| (s.as_str(), c.to_owned()))
            .collect();
        assert_eq!(result, vec![("Al", 2), ("Cl2", 3), ("AlCl3", 2)]);
    }

    #[test]
    // C6H5COOH + O2 = CO2 + H2O
    fn test_2() {
        let rg = vec![
            Substance::new("C6H5COOH", 3.0, None).unwrap(),
            Substance::new("O2", 3.0, None).unwrap(),
        ];
        let pd = vec![
            Substance::new("CO2", 3.0, None).unwrap(),
            Substance::new("H2O", 3.0, None).unwrap(),
        ];
        let bal = balance(rg, pd);
        let result: Vec<(&str, i32)> = bal
            .iter()
            .map(|(s, c)| (s.as_str(), c.to_owned()))
            .collect();
        assert_eq!(
            result,
            vec![("C6H5COOH", 2), ("O2", 15), ("CO2", 14), ("H2O", 6)]
        );
    }

    #[test]
    // KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2
    fn test_3() {
        let rg = vec![
            Substance::new("KMnO4", 3.0, None).unwrap(),
            Substance::new("HCl", 3.0, None).unwrap(),
        ];
        let pd = vec![
            Substance::new("KCl", 3.0, None).unwrap(),
            Substance::new("MnCl2", 3.0, None).unwrap(),
            Substance::new("H2O", 3.0, None).unwrap(),
            Substance::new("Cl2", 3.0, None).unwrap(),
        ];
        let bal = balance(rg, pd);
        let result: Vec<(&str, i32)> = bal
            .iter()
            .map(|(s, c)| (s.as_str(), c.to_owned()))
            .collect();
        assert_eq!(
            result,
            vec![
                ("KMnO4", 2),
                ("HCl", 16),
                ("KCl", 2),
                ("MnCl2", 2),
                ("H2O", 8),
                ("Cl2", 5)
            ]
        );
    }
}
