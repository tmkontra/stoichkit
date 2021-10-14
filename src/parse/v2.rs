use std::collections::hash_map::RandomState;
use std::collections::HashMap;

use itertools::fold;
use nom::branch::alt;
use nom::bytes::complete::{take_until, take_while, take_while1};
use nom::character::complete::{digit1, one_of};
use nom::combinator::map;
use nom::multi::many0;
use nom::{
    bytes::complete::{tag, take_while_m_n},
    character::complete::char,
    combinator::{map_res, opt},
    multi::many1,
    number::complete::be_u64,
    sequence::delimited,
    sequence::pair,
    sequence::tuple,
    AsChar, IResult, Parser,
};
use periodic_table_on_an_enum::Element::Hydrogen;

use crate::model::Element;
use crate::parse::v1::element_from_string;

/// Parse an elemental symbol, e.g. H, He, Na, S, Co
fn symbol(sym: &str) -> IResult<&str, Element> {
    map_res(
        map(
            // one capital letter and possibly a second lower-case letter
            pair(
                take_while_m_n(1, 1, |c: char| {
                    c.to_owned().is_ascii_uppercase()
                }),
                opt(take_while1(|c: char| c.to_owned().is_ascii_lowercase())),
            ),
            |(one, two): (&str, Option<&str>)| {
                format!("{}{}", one, two.unwrap_or(""))
            },
        ),
        |s| element_from_string(&s),
    )(sym)
}

/// Parse the multiplier following a symbol, must be an integer (1 or more digits)
fn multiplier(mult: &str) -> IResult<&str, u64> {
    map_res(digit1, |s: &str| s.parse::<u64>())(mult)
}

/// Parse an element with its multiplier
fn element(elem: &str) -> IResult<&str, (Element, u64)> {
    // no multiplier means 1
    pair(symbol, opt(multiplier).map(|m| m.unwrap_or(1)))(elem)
}

/// Sums a vector of elements and their multipliers,
/// rolling them into a map
fn sum_elements(vec: Vec<(Element, u64)>) -> HashMap<Element, u64> {
    vec.into_iter()
        .fold(HashMap::new(), |mut acc, (element, num)| {
            let entry: &mut u64 = acc.entry(element.to_owned()).or_insert(0);
            *entry += num;
            acc
        })
}

/// Parses a group, which is a sequence of elements, each optionally
/// followed by a multiplier, e.g. COOH or CO2H or NaCl2
fn group(group: &str) -> IResult<&str, HashMap<Element, u64>> {
    map(many1(element), sum_elements)(group)
}

/// Parses a multi-group, which is a parenthetical group
/// optionally followed by a multiplier, e.g. (SO4) or (SO4)2
fn multi_group(multi_group: &str) -> IResult<&str, HashMap<Element, u64>> {
    map(
        pair(
            delimited(one_of("([{"), group, one_of(")]}")),
            opt(multiplier).map(|m| m.unwrap_or(1)),
        ),
        |(group, multiplier)| {
            group.to_owned().iter().fold(
                HashMap::new(),
                |mut acc, (el, num)| {
                    let entry: &mut u64 = acc.entry(el.to_owned()).or_insert(0);
                    *entry += num * multiplier;
                    acc
                },
            )
        },
    )(multi_group)
}

/// Parses an asterisk followed by "H2O"
fn hydrate(hydrate: &str) -> IResult<&str, HashMap<Element, u64>> {
    let (hydrate, _) = tag("*")(hydrate)?;
    map(
        pair(opt(multiplier).map(|m| m.unwrap_or(1)), tag("H2O")),
        |(mx, _hydrate)| {
            let mut map = HashMap::new();
            map.insert(Element::from_symbol("H").unwrap(), 2 * mx);
            map.insert(Element::from_symbol("O").unwrap(), mx);
            map
        },
    )(hydrate)
}

/// Sums a vector of group results
fn sum_groups(vec: Vec<HashMap<Element, u64>>) -> HashMap<Element, u64> {
    vec.into_iter().fold(HashMap::new(), |mut acc, m| {
        for (el, mul) in m {
            let c = acc.entry(el).or_insert(0);
            *c += mul;
        }
        acc
    })
}

/// Parses the full formula
/// Must contain at least one group or multi-group, e.g. H2O or (SO4)2
/// May also have a "hydrate" suffix, e.g. *6H2O
fn formula_parser(formula: &str) -> IResult<&str, HashMap<Element, u64>> {
    map(
        map(
            pair(many1(alt((multi_group, group))), opt(hydrate)),
            |(groups, maybe_hydrate)| match &maybe_hydrate {
                None => groups,
                Some(hydrate) => {
                    let mut hydrate_group = groups.to_owned();
                    let mut h = vec![hydrate.to_owned()].to_owned();
                    hydrate_group.append(h.as_mut());
                    hydrate_group
                }
            },
        ),
        sum_groups,
    )(formula)
}

/// Parse a compound formula string.
/// # Examples
/// parse_formula_v2("H2O")
/// parse_formula_v2("H2(SO4)2")
pub fn parse_formula_v2(
    formula: &str,
) -> Result<HashMap<Element, u64>, String> {
    let (rem, elems) =
        formula_parser(formula).map_err(|_e| "failed!".to_string())?;
    Ok(elems)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::model::Element;
    use crate::parse::v2::parse_formula_v2;

    #[test]
    fn test_H2O() {
        let formula = "H2O";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("H").unwrap(), 2);
        exp.insert(Element::from_symbol("O").unwrap(), 1);
        assert_eq!(map, exp);
    }

    #[test]
    fn test_multi_group() {
        let formula = "H2(SO4)2";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("H").unwrap(), 2);
        exp.insert(Element::from_symbol("S").unwrap(), 2);
        exp.insert(Element::from_symbol("O").unwrap(), 8);
        assert_eq!(map, exp);
    }

    #[test]
    fn test_hydrate() {
        let formula = "H2(SO4)2*6H2O";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("H").unwrap(), 14);
        exp.insert(Element::from_symbol("S").unwrap(), 2);
        exp.insert(Element::from_symbol("O").unwrap(), 14);
        assert_eq!(map, exp);
    }

    #[test]
    fn test_organic() {
        let formula = "C6H5COOH";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("C").unwrap(), 7);
        exp.insert(Element::from_symbol("H").unwrap(), 6);
        exp.insert(Element::from_symbol("O").unwrap(), 2);
        assert_eq!(map, exp);
    }

    #[test]
    fn test_brackets() {
        let formula = "H2[SO4]2";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("H").unwrap(), 2);
        exp.insert(Element::from_symbol("S").unwrap(), 2);
        exp.insert(Element::from_symbol("O").unwrap(), 8);
        assert_eq!(map, exp);
    }

    #[test]
    fn test_curly_brace() {
        let formula = "H2{SO4}2";
        let map = parse_formula_v2(formula).unwrap();
        let mut exp = HashMap::new();
        exp.insert(Element::from_symbol("H").unwrap(), 2);
        exp.insert(Element::from_symbol("S").unwrap(), 2);
        exp.insert(Element::from_symbol("O").unwrap(), 8);
        assert_eq!(map, exp);
    }
}
