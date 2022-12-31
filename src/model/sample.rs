use crate::model::Reactant;

#[derive(Clone, Debug)]
pub struct Sample {
    pub reactant: Reactant,
    pub mass: f32,
}

impl Sample {
    pub fn of_reactant(reactant: Reactant, mass: f32) -> Self {
        Sample { reactant, mass }
    }

    pub fn from_formula(
        formula: &str,
        mass: f32,
        molar_coefficient: usize,
    ) -> Result<Sample, String> {
        let rct = Reactant::from_formula(formula, molar_coefficient);
        Ok(Sample {
            reactant: rct?,
            mass,
        })
    }

    pub fn moles(&self) -> f32 {
        self.mass / self.reactant.compound.molar_mass
    }

    pub fn molrxn(&self) -> f32 {
        self.moles() / self.reactant.molar_coefficient as f32
    }
}
