use std::cmp::max;

use itertools::Itertools;
use nalgebra::DMatrix;
use num::integer::lcm;
use rug::Rational;

use crate::model::Compound;
use crate::model::Element;

pub fn build_matrix(
    all_elements: &Vec<&Element>,
    all_compounds: Vec<&Compound>,
    ncols: usize,
) -> DMatrix<f64> {
    let mut matrix: Vec<f64> = Vec::new();
    debug!("Building matrix");
    for element in all_elements {
        debug!("Getting coefficients for {:?}", element);
        for compound in all_compounds.clone() {
            let coefficient =
                compound.atoms.get(element).cloned().unwrap_or(0 as usize);
            trace!(
                "Pushing {:?}*{:?} from {:?}",
                coefficient,
                element,
                compound
            );
            matrix.push(coefficient as f64);
        }
    }
    debug!("Constructing matrix");
    DMatrix::from_row_slice(all_elements.len(), ncols, matrix.as_slice())
}

pub fn solve_system(
    mx: DMatrix<f64>,
    ncols: usize,
) -> Result<Vec<f64>, String> {
    let (a, b) = mx.columns_range_pair(0..ncols, ncols..);
    let solution = if !a.is_square() {
        debug!("Solving non-square matrix by SVD");
        let x = a.svd(true, true);
        debug!("Solving equation system");
        x.solve(&b, 0.0)?
    } else {
        debug!("Solving square matrix by LU");
        let x = a.lu();
        debug!("Solving equation system");
        x.solve(&b)
            .expect(format!("Failed to solve matrix! {:?}", x).as_str())
    };
    let coefficients =
        solution.column(0).iter().map(|c| c.to_owned()).collect();
    Ok(coefficients)
}

pub fn normalize_coefficients(
    coefficients: Vec<f64>,
) -> Result<Vec<usize>, String> {
    let rational_coefficients: Vec<Rational> = coefficients
        .iter()
        .map(|c| {
            trace!("Constructing rational from {:?}", &c);
            Rational::from_f64(c.to_owned()).ok_or_else(|| {
                format!("Could not construct rational from: {:?}", &c)
            })
        })
        .collect::<Result<Vec<Rational>, String>>()?;
    trace!("Got rational coefficients: {:?}", &rational_coefficients);
    trace!("Limiting denominators");
    let rational_limited: Vec<Rational> = rational_coefficients
        .iter()
        .map(|r| limit_denominator(&r.clone().abs(), 100))
        .collect::<Result<Vec<Rational>, String>>()?;
    trace!("Limited rational coefficients: {:?}", rational_limited);
    let denominators: Vec<usize> = rational_limited
        .iter()
        .map(|c| {
            c.denom().to_usize().ok_or_else(|| {
                format!(
                    "Could not convert denominator {:?} to u64!",
                    &c.denom()
                )
            })
        })
        .collect::<Result<Vec<usize>, String>>()?;
    trace!("Got denominators: {:?}", &denominators);
    let scale: usize = denominators
        .iter()
        .cloned()
        .combinations(2)
        .fold(1, |acc, cur| {
            max(acc, lcm(cur[0] as usize, cur[1] as usize))
        });
    debug!("Scaling coefficients by: {}", scale);
    let mut scaled_coefficients: Vec<usize> = rational_limited
        .iter()
        .map(|f| f * Rational::from((scale, 1)))
        .map(|f| {
            f.numer().to_usize().ok_or_else(|| {
                format!("Could not convert scaled {:?} to u64", &f.numer())
            })
        })
        .collect::<Result<Vec<usize>, String>>()?;
    scaled_coefficients.push(scale);
    trace!("Got scaled coefficients: {:?}", scaled_coefficients);
    Ok(scaled_coefficients)
}

pub fn limit_denominator(
    given: &Rational,
    max_denominator: u64,
) -> Result<Rational, String> {
    debug!(
        "Limiting denominator for {:?} to at most {:?}",
        given, max_denominator
    );
    if max_denominator < 1 || given.denom() <= &max_denominator {
        Ok(given.to_owned())
    } else {
        let (mut p0, mut q0, mut p1, mut q1) = (0u64, 1u64, 1u64, 0u64);
        let (mut n, mut d) = (
            given
                .numer()
                .to_u64()
                .ok_or_else(|| format!("No numerator for: {:?}", given))?,
            given
                .denom()
                .to_u64()
                .ok_or_else(|| format!("No denominator for: {:?}", given))?,
        );
        let mut a: u64;
        let mut q2: u64;
        loop {
            a = n / d;
            q2 = q0 + (a * q1);
            if q2 > max_denominator {
                break;
            }
            let p0_new: u64 = p1;
            let q0_new: u64 = q1;
            let p1_new: u64 = p0 + (a * p1);
            let q1_new: u64 = q2;
            p0 = p0_new;
            q0 = q0_new;
            p1 = p1_new;
            q1 = q1_new;
            let n_new: u64 = d;
            let d_new: u64 = n - (a * d);
            n = n_new;
            d = d_new;
        }
        let k = (max_denominator - q0) / q1;
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
#[allow(non_snake_case)]
mod tests {
    use crate::model::Reaction;
    use crate::model::*;

    macro_rules! parse_balanced_reagent {
        (($subst:tt, $coef: tt)) => {
            Reactant::from_formula(stringify!($subst), $coef).unwrap()
        };
        ($subst:tt) => {
            Reactant::from_formula(stringify!($subst), 1).unwrap()
        };
    }

    macro_rules! new_reaction {
        ($reactSubst:tt $( + $reactSubstTail:tt)* = $prodSubst:tt $( + $prodSubstTail:tt)*) => {{
            // Al + Cl2 = AlCl3
            let reagents = vec![stringify!($reactSubst) $(,stringify!($reactSubstTail))*];
            let products = vec![stringify!($prodSubst) $(,stringify!($prodSubstTail))*];
            let r = reagents.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            let p = products.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            Reaction::new(r, p)
        }};
    }

    /**
       e.x. expect_balanced!(H2O = H2 + O2 => (H2O, 2) = (H2, 2) + O2)
    */
    macro_rules! expect_balanced {
        ($reactSubst:tt $( + $reactSubstTail:tt)* = $prodSubst:tt $( + $prodSubstTail:tt)* => $expectedReactant:tt $( + $expectedReactantTail:tt)* = $expectedProduct:tt $( + $expectedProductTail:tt)*) => {
            // Al + Cl2 = AlCl3
            let reagents = vec![stringify!($reactSubst) $(,stringify!($reactSubstTail))*];
            let products = vec![stringify!($prodSubst) $(,stringify!($prodSubstTail))*];
            let exp_reag = vec![parse_balanced_reagent!($expectedReactant) $(, parse_balanced_reagent!($expectedReactantTail))*];
            let exp_prod = vec![parse_balanced_reagent!($expectedProduct) $(, parse_balanced_reagent!($expectedProductTail))*];
            let r = reagents.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            let p = products.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            let reaction = Reaction::new(r, p).unwrap();
            let solution = reaction.balance().unwrap();
            assert_eq!(solution, BalancedReaction::new(exp_reag, exp_prod).unwrap())
        };
    }

    macro_rules! expect_coefficients {
        ($reactSubst:tt $( + $reactSubstTail:tt)* = $prodSubst:tt $( + $prodSubstTail:tt)* => $expectedReactant:literal $( + $expectedReactantTail:literal)*) => {
            // Al + Cl2 = AlCl3
            let reagents = vec![stringify!($reactSubst) $(,stringify!($reactSubstTail))*];
            let products = vec![stringify!($prodSubst) $(,stringify!($prodSubstTail))*];
            let exp_reag = vec![($expectedReactant) $(, ($expectedReactantTail))*];
            let r = reagents.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            let p = products.into_iter().map(|r| Compound::from_formula(r).unwrap()).collect();
            let reaction = Reaction::new(r, p).unwrap();
            let solution = reaction.balance().unwrap();
            assert_eq!(solution.all_coefficients(), exp_reag);
        };
    }

    fn _formulas_to_compounds(formulas: Vec<&str>) -> Vec<Compound> {
        formulas
            .iter()
            .cloned()
            .map(|f| Compound::from_formula(f).unwrap())
            .collect()
    }

    #[test]
    fn test_AlCl3() {
        expect_balanced!(
            Al + Cl2 = AlCl3 =>
                (Al, 2) + (Cl2, 3) = (AlCl3, 2)
        );
    }

    #[test]
    // C6H5COOH + O2 = CO2 + H2O
    fn test_CO2() {
        expect_balanced!(
            C6H5COOH + O2 = CO2 + H2O =>
                (C6H5COOH, 2) + (O2, 15) = (CO2, 14) + (H2O, 6)
        );
    }

    #[test]
    // KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2
    fn test_KMnO4() {
        expect_balanced!(
            KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2 =>
                (KMnO4, 2) + (HCl, 16) = (KCl, 2) + (MnCl2, 2) + (H2O, 8) + (Cl2, 5)
        );
    }

    #[test]
    fn test_H2O() {
        expect_balanced!(
            H2 + O2 = H2O =>
                (H2, 2) + (O2, 1) = (H2O, 2)
        );
    }

    #[test]
    fn test_missing_products() {
        // Fe3 + Cl5 = Cl2Fe5H2O
        let rxn: Result<Reaction, String> =
            new_reaction!(Fe3 + Cl5 = Cl2Fe5H2O);
        assert!(
            rxn.is_err(),
            format!("Balance solution was not Err: {:?}", rxn),
        )
    }

    #[test]
    fn test_impossible_reaction() {
        // H2O + NO2 = HNO3
        let rxn = new_reaction!(H2O + NO2 = HNO3).unwrap();
        let result = rxn.balance();
        assert!(
            result.is_err(),
            format!("Balance solution was not Err: {:?}", result),
        )
    }

    #[test]
    fn test_impossible_reaction_2() {
        // some coefficients from matrix solution will be 0
        let rxn = new_reaction!(Cu + HNO3 = CuN2O6 + NO2 + H2O2).unwrap();
        let result = rxn.balance();
        assert!(
            result.is_err(),
            format!("Balance solution was not Err: {:?}", result)
        )
    }

    #[test]
    fn test_batch() {
        expect_coefficients!(CO2 + H2O = C6H12O6 + O2 => 6 + 6 + 1 + 6);
        expect_coefficients!(SiCl4 + H2O = H4SiO4 + HCl => 1 + 4 + 1 + 4);
        expect_coefficients!(Al + HCl = AlCl3 + H2 => 2 + 6 + 2 + 3);
        expect_coefficients!(Na2CO3 + HCl = NaCl + H2O + CO2 => 1 + 2 + 2 + 1 + 1);
        expect_coefficients!(C7H6O2 + O2 = CO2 + H2O => 2 + 15 + 14 + 6);
        expect_coefficients!(Fe2S3O12 + KOH = K2SO4 + FeO3H3 => 1 + 6 + 3 + 2);
        expect_coefficients!(Ca3P2O8 + SiO2 = P4O10 + CaSiO3 => 2 + 6 + 1 + 6);
        expect_coefficients!(KClO3 = KClO4 + KCl => 4 + 3 + 1);
        expect_coefficients!(Al2S3O12 + CaO2H2 = AlO3H3 + CaSO4 => 1 + 3 + 2 + 3);
        expect_coefficients!(H2SO4 + HI = H2S + I2 + H2O => 1 + 8 + 1 + 4 + 4);
    }
}
