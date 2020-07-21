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

use ndarray::prelude::*;
use ndarray_linalg::Solve;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use num::abs;


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
    println!("Found elements: {:?}", elements);
    println!("Matrix Values: {:?}", &matrix);
    let mx = DMatrix::from_row_slice(elements.len(), reagents.len() + products.len(), matrix.as_slice());
    let c = reagents.len() + products.len() - 1;
    let (a, b) = mx.columns_range_pair(0..c, c..);
    println!("A =  {:?}", a);
    println!("b = {:?}", b);
    let x = a.svd(true, true);
    println!("SVD: {:?}", x);
    let solve = x.solve(&b, 0.0).unwrap();
    println!("x = {:?}", solve);
    let c: Vec<f32> = solve.column(0).iter().map(|c| *c as f32).collect();
    println!("Coeffs {:?}", c);
    let fc: Vec<Rational> = c.iter().map(|c| Rational::from_f32(*c).unwrap()).collect();
    println!("Rational Coeffs: {:?}", fc);
    let denoms: Vec<_> = fc.iter().map(|c| c.denom()).collect();
    println!("Denoms: {:?}", denoms);
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
    // let rounded = roundr(c[0]);
    // println!("ROUNDED: {:?}", rounded);
    let p = limit_denominator(&fc[0], 10);
    println!("LimDenom: {:?}", p);
    let py = pymain(c[0].clone());
    println!("PyResult: {:?}", py.unwrap());
}

fn limit_denominator(f: &Rational, max_denom: u32) -> Rational {
    if max_denom < 1 || f.denom() < &max_denom {
        f.clone()
    } else {
        let (mut p0, mut q0, mut p1, mut q1) = (0, 1, 1, 0);
        let (mut n, mut d) = (f.numer().to_i32().unwrap(), f.denom().to_i32().unwrap());
        let mut a: i32;
        let mut q2: i32;
        loop {
            a = n / d;
            q2 = q0+a*q1;
            if q2 > max_denom as i32 {
                break
            }
            let p0new = p1.clone();
            let q0new = q1.clone();
            let p1new = (p0.clone()+a.clone()*p1.clone());
            let q1new = q2.clone();
            p0 = p0new;
            q0 = q0new;
            p1 = p1new;
            q1 = q1new;
            let nnew = d.clone();
            let dnew = (n.clone() - a.clone() * d.clone());
            n = nnew;
            d = dnew;
        }
        let k = (max_denom as i32 - q0) / q1;
        let bound1: Rational = Rational::from((p0+k*p1, q0+k*q1));
        let bound2: Rational = Rational::from((p1, q1));
        println!("Bound1: {:?} or Bound2: {:?}", &bound1, &bound2);
        let d1f: Rational = (bound1.clone() - f);
        let d2f: Rational = bound2.clone() - f;
        println!("D1: {:?} -- D2: {:?}", d1f.clone().abs(), d2f.clone().abs());
        return if d1f.abs() <= d2f.abs() {
            bound1.clone()
        } else {
            bound2.clone()
        }
    }
}

fn pymain(v: f32) -> Result<(i32, i32), ()> {
    let gil = Python::acquire_gil();
    let py = gil.python();
    pymain_(py, v).map_err(|e| {
        // We can't display Python exceptions via std::fmt::Display,
        // so print the error here manually.
        e.print_and_set_sys_last_vars(py);
    })
}

fn pymain_(py: Python, v: f32) -> PyResult<(i32, i32)> {
    let globals = [
        ("fractions", py.import("fractions")?),
    ].into_py_dict(py);
    let locals = [
        ("v", v),
    ].into_py_dict(py);
    let code = "(fractions.Fraction(v).limit_denominator(100).numerator, fractions.Fraction(v).limit_denominator(100).denominator)";
    let ret: (i32, i32) = py.eval(code, Some(&globals), Some(&locals))?.extract()?;
    println!("Got {} / {}", ret.0, ret.1);
    Ok(ret)
}


fn roundr(f: f32) -> Result<(i32, i32), PyErr>{
    let gil = Python::acquire_gil();
    let py = gil.python();
    let activators = PyModule::from_code(py, r#"
    def nd(x):
        from fractions import Fraction
        f = Fraction(x).limit_denominator(64)
        return (f.numerator, f.denominator)
    "#, "activators.py", "activators").unwrap();
    let total: (i32, i32) = activators.call1("nd", (f,))?.extract()?;
    println!("{:?}", total);
    Ok(total)
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

