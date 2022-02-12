use crate::model::{
    yield_reaction, BalancedReaction, Reactant, Sample, Units, YieldReaction,
};
use clap::ArgEnum;
use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub struct TheoreticalReaction {
    pub reaction: BalancedReaction,
    pub reactants: Vec<Sample>,
}

#[derive(ArgEnum, Clone)]
pub enum YieldUnits {
    Mass,
    Moles,
}

impl Into<Units> for YieldUnits {
    fn into(self) -> Units {
        match self {
            YieldUnits::Mass => Units::Grams,
            YieldUnits::Moles => Units::Moles,
        }
    }
}

impl TheoreticalReaction {
    pub fn yields(&self, units: &YieldUnits) -> Vec<(Reactant, f32)> {
        let limiting =
            yield_reaction::limiting_reagent(self.reactants.to_owned());
        self.reaction
            .products
            .iter()
            .map(|p| yield_reaction::theoretical_yield(&limiting, p))
            .zip(self.reaction.products.to_owned())
            .map(|(moles, product)| match units {
                YieldUnits::Mass => (
                    product.to_owned(),
                    moles * product.compound.molar_mass.to_owned(),
                ),
                YieldUnits::Moles => (product, moles),
            })
            .collect()
    }
}

impl TheoreticalReaction {
    pub fn new(
        reaction: BalancedReaction,
        reactants: Vec<Sample>,
    ) -> TheoreticalReaction {
        TheoreticalReaction {
            reaction,
            reactants,
        }
    }
}
