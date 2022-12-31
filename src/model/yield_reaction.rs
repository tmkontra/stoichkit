use std::cmp::Ordering;

use crate::model::{Reactant, Sample};

#[derive(Debug, Clone)]
pub struct YieldReaction {
    pub reagents: Vec<Sample>,
    pub product: Sample,
}

pub fn limiting_reagent(reagents: Vec<Sample>) -> Sample {
    reagents
        .into_iter()
        .min_by(|l, r| {
            l.molrxn()
                .partial_cmp(&r.molrxn())
                .unwrap_or(Ordering::Equal)
        })
        .map(|s| {
            debug!("Limiting reagent is {}", s.reactant.compound.formula);
            s
        })
        .unwrap()
}

pub fn theoretical_yield(limiting: &Sample, product: &Reactant) -> f32 {
    trace!("{} moles of limiting reagent", limiting.moles());
    let exp_moles = limiting.moles()
        * (product.molar_coefficient as f32
            / limiting.reactant.molar_coefficient as f32);
    debug!("Theoretical moles of product: {}", exp_moles);
    exp_moles
}

impl YieldReaction {
    pub fn new(reagents: Vec<Sample>, product: Sample) -> YieldReaction {
        YieldReaction { reagents, product }
    }

    pub fn limiting_reagent(self: &Self) -> Sample {
        limiting_reagent(self.reagents.to_owned())
    }

    pub fn theoretical_yield(self: &Self) -> f32 {
        let limiting = self.limiting_reagent();
        theoretical_yield(&limiting, &self.product.reactant)
            * self.product.reactant.compound.molar_mass
    }

    pub fn percent_yield(self: &Self) -> f32 {
        self.product.mass / self.theoretical_yield()
    }
}
