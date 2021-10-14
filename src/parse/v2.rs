use crate::model::Element;
use itertools::fold;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

use crate::parse::element_from_string;
use nom::branch::alt;
use nom::bytes::complete::{take_until, take_while, take_while1};
use nom::character::complete::digit1;
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

fn sym2(sym: &str) -> IResult<&str, String> {
    map(
        pair(
            take_while_m_n(1, 1, |c: char| c.to_owned().is_ascii_uppercase()),
            opt(take_while1(|c: char| c.to_owned().is_ascii_lowercase())),
        ),
        |(one, two): (&str, Option<&str>)| {
            format!("{}{}", one, two.unwrap_or(""))
        },
    )(sym)
}

fn element_from_String(string: String) -> Result<Element, String> {
    element_from_string(&string)
}

fn sym(sym: &str) -> IResult<&str, Element> {
    map_res(sym2, element_from_String)(sym)
}

fn mult(mult: &str) -> IResult<&str, u64> {
    map_res(digit1, |s: &str| s.parse::<u64>())(mult)
}

fn elem(elem: &str) -> IResult<&str, (Element, Option<u64>)> {
    pair(sym, opt(mult))(elem)
}

fn elems_to_map(vec: Vec<(Element, Option<u64>)>) -> HashMap<Element, u64> {
    vec.into_iter()
        .map(|(el, maybeMult)| (el, maybeMult.unwrap_or(1)))
        .collect()
}

fn group(group: &str) -> IResult<&str, HashMap<Element, u64>> {
    map(many1(elem), elems_to_map)(group)
}

fn multi_group(multi_group: &str) -> IResult<&str, HashMap<Element, u64>> {
    map(
        pair(delimited(tag("("), group, tag(")")), opt(mult)),
        |(g, m)| {
            let mx: u64 = m.unwrap_or(1);
            g.to_owned()
                .iter()
                .fold(HashMap::new(), |mut acc, (el, num)| {
                    let entry: &mut u64 = acc.entry(el.to_owned()).or_insert(0);
                    *entry += num * mx;
                    acc
                })
        },
    )(multi_group)
}

fn hydrate(hydrate: &str) -> IResult<&str, HashMap<Element, u64>> {
    let (hydrate, _) = tag("*")(hydrate)?;
    map(pair(opt(mult), tag("H2O")), |(maybeMult, _H2O)| {
        let mx = maybeMult.unwrap_or(1);
        let mut map = HashMap::new();
        map.insert(Element::from_symbol("H").unwrap(), 2 * mx);
        map.insert(Element::from_symbol("O").unwrap(), mx);
        map
    })(hydrate)
}

fn fold_maps(vec: Vec<HashMap<Element, u64>>) -> HashMap<Element, u64> {
    vec.into_iter().fold(HashMap::new(), |mut acc, m| {
        for (el, mul) in m {
            let c = acc.entry(el).or_insert(0);
            *c += mul;
        }
        acc
    })
}

fn formula_parser(formula: &str) -> IResult<&str, Vec<HashMap<Element, u64>>> {
    map(
        pair(many1(alt((multi_group, group))), opt(hydrate)),
        |(groups, maybeHydrate)| match &maybeHydrate {
            None => groups,
            Some(hydrate) => {
                let mut g = groups.to_owned();
                let mut h = vec![hydrate.to_owned()].to_owned();
                g.append(h.as_mut());
                g
            }
        },
    )(formula)
}

pub fn parse_formula_v2(
    formula: &str,
) -> Result<HashMap<Element, u64>, String> {
    let (rem, elems) =
        formula_parser(formula).map_err(|_e| "failed!".to_string())?;
    Ok(fold_maps(elems))
}

#[cfg(test)]
mod tests {
    use crate::model::Element;
    use crate::parse::v2::parse_formula_v2;
    use std::collections::HashMap;

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
}
