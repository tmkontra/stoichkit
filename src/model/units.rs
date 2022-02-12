use std::fmt::{Display, Formatter};

pub enum Units {
    Grams,
    Percent,
    Moles,
}

impl Display for Units {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Units::Grams => "g",
            Units::Percent => "%",
            Units::Moles => "mol",
        };
        write!(f, "{}", str)
    }
}
