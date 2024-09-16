use anyhow::Result;
use fraction::{One, Zero};

use crate::{activity_key::ActivityKeyTranslator, ebi_traits::{ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage, ebi_trait_queriable_stochastic_language::EbiTraitQueriableStochasticLanguage}, follower_semantics::FollowerSemantics, math::fraction::Fraction};

pub fn uemsc(log1: Box<dyn EbiTraitFiniteStochasticLanguage>, mut language2: Box<dyn EbiTraitQueriableStochasticLanguage>) -> Result<Fraction> {
    let mut sum = Fraction::zero();

    for (trace, probability1) in log1.iter_trace_probability() {
        let translator = ActivityKeyTranslator::new(log1.get_activity_key(), language2.get_activity_key_mut());
        let probability2 = language2.get_probability(&FollowerSemantics::Trace(&translator.translate_trace(trace)))?;

        // if !probability2.is_zero() {
            log::debug!("trace {:?} probability in model {}", trace, probability2);
        // }

        if *probability1 > probability2 {
            sum += probability1;
            sum -= probability2;
        }
    }

    sum = sum.one_minus();
    return Ok(sum);
}