use crate::model::{
    yield_reaction, BalancedReaction, Reactant, Sample, Units,
};
use clap::ArgEnum;

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

impl From<YieldUnits> for Units {
    fn from(val: YieldUnits) -> Self {
        match val {
            YieldUnits::Mass => Units::Grams,
            YieldUnits::Moles => Units::Moles,
        }
    }
}

impl TheoreticalReaction {
    pub fn yields(&self, units: &YieldUnits) -> Vec<(&Reactant, f32)> {
        let limiting =
            yield_reaction::limiting_reagent(&self.reactants);
        self.reaction
            .products
            .iter()
            .map(|p| yield_reaction::theoretical_yield(&limiting, p))
            .zip(&self.reaction.products)
            .map(|(moles, product)| match units {
                YieldUnits::Mass => (
                    product,
                    moles * product.compound.molar_mass,
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
