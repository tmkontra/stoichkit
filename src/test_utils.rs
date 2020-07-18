use std::collections::HashMap;


pub fn e(expected: HashMap<&str, u32>) -> HashMap<String, u32> {
    expected.iter().map(|p| (p.0.to_string(), p.1.to_owned())).collect()
}
