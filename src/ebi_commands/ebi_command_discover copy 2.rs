use crate::{
    ebi_framework::{
        activity_key::HasActivityKey, ebi_command::EbiCommand, ebi_input::EbiInputType, ebi_object::{EbiObject, EbiObjectType}, ebi_output::{EbiOutput, EbiOutputType}, ebi_trait::EbiTrait
    }, ebi_objects::{finite_stochastic_language::FiniteStochasticLanguage, labelled_petri_net::LabelledPetriNet, process_tree::ProcessTree, stochastic_transition_system::StochasticTransitionSystem}, ebi_traits::{ebi_trait_event_log::EbiTraitEventLog, ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage}, math::fraction, techniques::{
        alignment_stochastic_miner::AlignmentMiner, occurrences_stochastic_miner::{
            OccurrencesStochasticMinerLPN, OccurrencesStochasticMinerTree,
        }, sdfa_discovery::convert_log_to_stochastic_transition_system, uniform_stochastic_miner::{UniformStochasticMinerLPN, UniformStochasticMinerTree}
    }
};
use std::thread;
use std::collections::HashMap;

pub const EBI_DISCOVER: EbiCommand = EbiCommand::Group {
    name_short: "disc",
    name_long: Some("discover"),
    explanation_short: "Discover a stochastic process model.",
    explanation_long: None,
    children: &[
        &EBI_DISCOVER_SDFA,
    ],
};

pub const EBI_DISCOVER_SDFA: EbiCommand = EbiCommand::Command { 
    name_short: "disc_sdfa", 
    name_long: Some("discover_sdfa"), 
    explanation_short: "Discover an sdfa to target size from the event log.", 
    explanation_long: None, 
    latex_link: None, 
    cli_command: None, 
    exact_arithmetic: false, 
    input_types: &[ 
        &[&EbiInputType::Trait(EbiTrait::EventLog)], 
        &[&EbiInputType::Usize],
    ], 
    input_names: &["EVENT_LOG_FILE", "TARGET_SIZE"], 
    input_helps: &[ "An event log.",  "The target size of the sdfa."], 
    execute: |mut inputs, _| {
        // get the event log
        let mut log: Box<dyn EbiTraitEventLog + 'static> = inputs.remove(0).to_type::<dyn EbiTraitEventLog>()?;

        // get the stochastic language of the log
        let trace_num = log.len();
        let mut multiset = log.to_multiset();

        let mut multiset_map = HashMap::new();
        // for each key in multiset, divide the frequency by the trace number
        for (key, value) in multiset.iter_mut() {
            let mut trace_prob = fraction::Fraction::from(value.clone());
            trace_prob /= fraction::Fraction::from(trace_num);
            multiset_map.insert(key.clone(), trace_prob);
        }

        let activity_key = log.get_activity_key().clone();
        let new_activity_key = log.get_activity_key().clone();
        let new_multiset_map = multiset_map.clone();
        let mut stochastic_language = FiniteStochasticLanguage::from((multiset_map.clone(), activity_key.clone()));

        log.translate_using_activity_key(stochastic_language.get_activity_key_mut());

        // get the target model size
        let target_size = *inputs.remove(0).to_type::<usize>()?;

        let mut sts = StochasticTransitionSystem::new();
        //create transition system
        for trace_index in 0..log.len() {
            let trace = log.get_trace(trace_index).unwrap(); 
            let mut state = sts.get_and_add_initial_state();

            for activity in trace {
                state = sts.take_or_add_transition(state, *activity, 1);
            }
        }
        for i in 0..sts.get_sources().len() {
            //add the source states
            let state = sts.get_sources()[i];
            let freq = sts.get_transition_frequencies()[i];
            sts.update_terminating_frequency(state, freq);
        }

        let size_before_filtering = sts.get_size();
        // filter the transition system by fourty percent
        sts.filter_by_percentage(40.0);
        sts.set_activity_key(stochastic_language.get_activity_key());
        let sdfa = convert_log_to_stochastic_transition_system(Box::new(stochastic_language), sts, target_size, size_before_filtering)?;
        Ok(EbiOutput::Object(EbiObject::StochasticDeterministicFiniteAutomaton(sdfa)))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticDeterministicFiniteAutomaton)
};