use std::collections::HashMap;

use crate::model::Element;

#[allow(dead_code)]
pub fn parse_elements(
    expected: HashMap<&str, usize>,
) -> HashMap<Element, usize> {
    expected
        .iter()
        .map(|p| (Element::from_symbol(p.0).unwrap(), p.1.to_owned()))
        .collect()
}
