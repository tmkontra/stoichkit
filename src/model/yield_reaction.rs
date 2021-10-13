use std::cmp::Ordering;

use crate::model::Sample;

#[derive(Debug, Clone)]
pub struct YieldReaction {
    pub reagents: Vec<Sample>,
    pub product: Sample,
}

impl YieldReaction {
    pub fn new(reagents: Vec<Sample>, product: Sample) -> YieldReaction {
        YieldReaction { reagents, product }
    }

    pub fn limiting_reagent(self: &Self) -> &Sample {
        self.reagents
            .iter()
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

    pub fn theoretical_yield(self: &Self) -> f32 {
        let limiting = self.limiting_reagent();
        trace!("{} moles of limiting reagent", limiting.moles());
        let exp_moles = limiting.moles()
            * (self.product.reactant.molar_coefficient as f32
            / limiting.reactant.molar_coefficient as f32);
        debug!("Theoretical moles of product: {}", exp_moles);
        let exp_grams = exp_moles * self.product.reactant.compound.molar_mass;
        debug!("Theoretical yield of product (g): {}", exp_grams);
        exp_grams
    }

    pub fn percent_yield(self: &Self) -> f32 {
        self.product.mass / self.theoretical_yield()
    }
}
