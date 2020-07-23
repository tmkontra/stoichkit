use std::collections::HashMap;

use ptable::Element;

#[allow(dead_code)]
pub fn e(expected: HashMap<&str, u32>) -> HashMap<Element, u32> {
    expected
        .iter()
        .map(|p| (Element::from_symbol(p.0).unwrap(), p.1.to_owned()))
        .collect()
}
