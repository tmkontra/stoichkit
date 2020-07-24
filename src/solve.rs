use std::cmp::max;
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use nalgebra::DMatrix;
use num::integer::lcm;
use ptable::Element;
use rug::Rational;

use crate::model::Substance;

pub fn balance(
    reagents: Vec<Substance>,
    products: Vec<Substance>,
) -> Result<Vec<(String, i32)>, String> {
    let mut reagent_atoms: HashSet<&Element> = HashSet::new();
    let mut product_atoms: HashSet<&Element> = HashSet::new();
    for r in &reagents {
        for e in r.atoms.keys() {
            reagent_atoms.insert(e);
        }
    }
    for p in &products {
        for e in p.atoms.keys() {
            product_atoms.insert(e);
        }
    }
    if !&reagent_atoms.eq(&product_atoms) {
        let missing_products: HashSet<_> = reagent_atoms.difference(&product_atoms).collect();
        let missing_reagents: HashSet<_> = product_atoms.difference(&reagent_atoms).collect();
        return Err(format!(
            "Equation cannot be balanced. Reagent elements that are not in products = {:?}. Product elements that are not in products = {:?}",
            missing_products, missing_reagents)
        );
    }
    let all_atoms: Vec<&Element> = reagent_atoms
        .iter()
        .chain(product_atoms.iter())
        .cloned()
        .collect();
    let elements: Vec<&Element> = all_atoms.iter().cloned().collect();
    let mut subs: Vec<&Substance> = Vec::new();
    let mut matrix: Vec<f64> = Vec::new();
    debug!("Building matrix");
    for element in all_atoms.iter() {
        debug!("Getting coefficients for {:?}", element);
        for rge in &reagents {
            subs.push(rge);
            let rxe = rge.atoms.get(element).cloned().unwrap_or(0 as u32);
            trace!("Pushing {:?}*{:?} from {:?}", rxe, element, rge);
            matrix.push(rxe as f64);
        }
        for pde in &products {
            subs.push(pde);
            let rxe = pde.atoms.get(element).cloned().unwrap_or(0 as u32);
            trace!("Pushing {:?}*{:?} from {:?}", rxe, element, pde);
            matrix.push(rxe as f64);
        }
    }
    debug!("Constructing matrix");
    let mx = DMatrix::from_row_slice(
        elements.len(),
        reagents.len() + products.len(),
        matrix.as_slice(),
    );
    let c = reagents.len() + products.len() - 1;
    let (a, b) = mx.columns_range_pair(0..c, c..);
    let x = a.svd(true, true);
    debug!("Solving equation system");
    let solve = x.solve(&b, 0.0).unwrap();
    let c: Vec<f32> = solve.column(0).iter().map(|c| *c as f32).collect();
    let fc: Vec<Rational> = c
        .iter()
        .cloned()
        .map(|c| Rational::from_f32(c).ok_or(format!("Could not rational: {:?}", c)))
        .collect::<Result<Vec<Rational>, String>>()?
        .iter()
        .cloned()
        .map(|r| limit_denominator(&r.abs(), 100))
        .collect::<Result<Vec<Rational>, String>>()?;
    let denoms: Vec<_> = fc.iter().map(|c| c.denom()).collect();
    let mult = denoms
        .iter()
        .cloned()
        .map(|i| i.to_i32().ok_or_else(|| "Could not balance!".to_string()))
        .collect::<Result<Vec<i32>, String>>()?
        .iter()
        .cloned()
        .combinations(2)
        .fold(1 as i32, |mult, cur| max(mult, lcm(cur[0], cur[1])));
    let mut fin: Vec<i32> = fc
        .iter()
        .map(|f| f * Rational::from((mult, 1)))
        .map(|f| {
            f.numer()
                .to_i32()
                .ok_or_else(|| "Could not i32".to_string())
        })
        .collect::<Result<Vec<i32>, String>>()?;
    fin.push(mult);
    let result: Vec<(String, i32)> = subs
        .iter()
        .map(|s| s.to_owned().formula.clone())
        .zip(&mut fin.iter().map(|c| c.to_owned()))
        .collect();
    let (reag, prod) = result.split_at(reagents.len());
    if check_balance(reag.to_vec(), prod.to_vec())? {
        Ok(result)
    } else {
        Err(format!("Equation could not be balanced!"))
    }
}

fn check_balance(
    reactants: Vec<(String, i32)>,
    products: Vec<(String, i32)>,
) -> Result<bool, String> {
    let react_subs: Vec<Substance> = reactants
        .iter()
        .map(|(f, c)| Substance::new(f.as_str(), 0.0, Some(*c as u32)))
        .collect::<Result<Vec<Substance>, String>>()?;
    let prod_subs: Vec<Substance> = products
        .iter()
        .map(|(f, c)| Substance::new(f.as_str(), 0.0, Some(*c as u32)))
        .collect::<Result<Vec<Substance>, String>>()?;
    let react_elems: HashMap<Element, u32> = react_subs
        .iter()
        .map(|s| (&s.atoms, s.molar_coefficient))
        .fold(HashMap::new(), |mut acc, (item, coeff)| {
            for (e, c) in item {
                let counter = acc.entry(e.to_owned()).or_insert(0);
                *counter += c * coeff;
            }
            acc
        });
    let prod_elems: HashMap<Element, u32> = prod_subs
        .iter()
        .map(|s| (&s.atoms, s.molar_coefficient))
        .fold(HashMap::new(), |mut acc, (item, coeff)| {
            for (e, c) in item {
                let counter = acc.entry(e.to_owned()).or_insert(0);
                *counter += c * coeff;
            }
            acc
        });
    debug!(
        "Checking balanced?: Reagent elements: {:?} === Product elements: {:?}",
        react_elems, prod_elems
    );
    Ok(react_elems.eq(&prod_elems))
}

fn limit_denominator(given: &Rational, max_denominator: u32) -> Result<Rational, String> {
    debug!(
        "Limiting denom for {:?} to at most {:?}",
        given, max_denominator
    );
    if max_denominator < 1 || given.denom() <= &max_denominator {
        Ok(given.to_owned())
    } else {
        let (mut p0, mut q0, mut p1, mut q1) = (0, 1, 1, 0);
        let (mut n, mut d) = (
            given
                .numer()
                .to_i32()
                .ok_or_else(|| format!("No numerator for: {:?}", given))?,
            given
                .denom()
                .to_i32()
                .ok_or_else(|| format!("No denominator for: {:?}", given))?,
        );
        let mut a: i32;
        let mut q2: i32;
        loop {
            a = n / d;
            q2 = q0 + (a * q1);
            if q2 > max_denominator as i32 {
                break;
            }
            let p0_new: i32 = p1;
            let q0_new: i32 = q1;
            let p1_new: i32 = p0 + (a * p1);
            let q1_new: i32 = q2;
            p0 = p0_new;
            q0 = q0_new;
            p1 = p1_new;
            q1 = q1_new;
            let n_new: i32 = d;
            let d_new: i32 = n - (a * d);
            n = n_new;
            d = d_new;
        }
        let k = (max_denominator as i32 - q0) / q1;
        let bound1: Rational = Rational::from((p0 + k * p1, q0 + k * q1));
        let bound2: Rational = Rational::from((p1, q1));
        let d1f: Rational = bound1.clone() - given;
        let d2f: Rational = bound2.clone() - given;
        if d1f.abs() <= d2f.abs() {
            Ok(bound1)
        } else {
            Ok(bound2)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::model::*;
    use crate::solve::balance;

    fn _formulas_to_substances(formulas: Vec<&str>) -> Vec<Substance> {
        formulas
            .iter()
            .cloned()
            .map(|f| Substance::new(f, 1.0, None).unwrap())
            .collect()
    }

    fn _expect_solution(reagents: Vec<&str>, products: Vec<&str>, expected: Vec<(&str, u32)>) {
        let solution = balance(
            _formulas_to_substances(reagents),
            _formulas_to_substances(products),
        )
        .unwrap();
        let result: Vec<(&str, u32)> = solution
            .iter()
            .map(|(s, c)| (s.as_str(), *c as u32))
            .collect();
        assert_eq!(result, expected)
    }

    #[test]
    fn test_AlCl3() {
        let rg = vec!["Al", "Cl2"];
        let pd = vec!["AlCl3"];
        let expected = vec![("Al", 2), ("Cl2", 3), ("AlCl3", 2)];
        _expect_solution(rg, pd, expected)
    }

    #[test]
    // C6H5COOH + O2 = CO2 + H2O
    fn test_CO2() {
        let rg = vec!["C6H5COOH", "O2"];
        let pd = vec!["CO2", "H2O"];
        let expected = vec![("C6H5COOH", 2), ("O2", 15), ("CO2", 14), ("H2O", 6)];
        _expect_solution(rg, pd, expected)
    }

    #[test]
    // KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2
    fn test_KMnO4() {
        let rg = vec!["KMnO4", "HCl"];
        let pd = vec!["KCl", "MnCl2", "H2O", "Cl2"];
        let expected = vec![
            ("KMnO4", 2),
            ("HCl", 16),
            ("KCl", 2),
            ("MnCl2", 2),
            ("H2O", 8),
            ("Cl2", 5),
        ];
        _expect_solution(rg, pd, expected)
    }

    #[test]
    fn test_H2O() {
        let rg = vec!["H2", "O2"];
        let pd = vec!["H2O"];
        let expected = vec![("H2", 2), ("O2", 1), ("H2O", 2)];
        _expect_solution(rg, pd, expected)
    }

    #[test]
    fn test_missing_products() {
        //Fe3 + Cl5 = Cl2Fe5H2O
        let rg = vec!["Fe3", "Cl5"];
        let pd = vec!["Cl2Fe5H2O"];
        let result = balance(_formulas_to_substances(rg), _formulas_to_substances(pd));
        assert!(
            result.is_err(),
            format!("Balance solution was not Err: {:?}", result),
        )
    }

    #[test]
    fn test_impossible_reaction() {
        // H2O + NO2 = HNO3
        let rg = vec!["H2O", "NO2"];
        let pd = vec!["HNO3"];
        let result = balance(_formulas_to_substances(rg), _formulas_to_substances(pd));
        assert!(
            result.is_err(),
            format!("Balance solution was not Err: {:?}", result),
        )
    }
}
