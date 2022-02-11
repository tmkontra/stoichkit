use std::cmp::Ordering;

use crate::model::{
    yield_reaction, BalancedReaction, Reactant, Sample, YieldReaction,
};

#[derive(Debug, Clone)]
pub struct TheoreticalReaction {
    pub reaction: BalancedReaction,
    pub reactants: Vec<Sample>,
}

impl TheoreticalReaction {
    pub fn yields(&self) -> Vec<f32> {
        let limiting =
            yield_reaction::limiting_reagent(self.reactants.to_owned());
        self.reaction
            .products
            .iter()
            .map(|p| yield_reaction::theoretical_yield(&limiting, p))
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
