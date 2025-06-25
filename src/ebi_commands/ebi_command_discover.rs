use crate::{
    ebi_framework::{
        activity_key::HasActivityKey, ebi_command::EbiCommand, ebi_input::EbiInputType, ebi_object::{EbiObject, EbiObjectType}, ebi_output::{EbiOutput, EbiOutputType}, ebi_trait::EbiTrait
    }, ebi_objects::{finite_stochastic_language::FiniteStochasticLanguage, labelled_petri_net::LabelledPetriNet, process_tree::ProcessTree, stochastic_transition_system::StochasticTransitionSystem}, ebi_traits::{ebi_trait_event_log::EbiTraitEventLog, ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage}, math::fraction, techniques::{
        alignment_stochastic_miner::AlignmentMiner, occurrences_stochastic_miner::{
            OccurrencesStochasticMinerLPN, OccurrencesStochasticMinerTree,
        }, sdfa_discovery::convert_log_to_stochastic_transition_system, uniform_stochastic_miner::{UniformStochasticMinerLPN, UniformStochasticMinerTree}
    }
};
use std::collections::HashMap;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rand::Rng;

pub const EBI_DISCOVER: EbiCommand = EbiCommand::Group {
    name_short: "disc",
    name_long: Some("discover"),
    explanation_short: "Discover a stochastic process model.",
    explanation_long: None,
    children: &[
        &EBI_DISCOVER_ALIGNMENTS,
        &EBI_DISCOVER_OCCURRENCE,
        &EBI_DISCOVER_UNIFORM,
        &EBI_DISCOVER_SDFA,
        &EBI_DISCOVER_SDFA_PARALLEL,
    ],
};

pub const EBI_DISCOVER_ALIGNMENTS: EbiCommand = EbiCommand::Command {
    name_short: "ali",
    name_long: Some("alignments"),
    explanation_short: "Give each transition a weight that matches the aligned occurrences of its label. The model must be livelock-free.",
    explanation_long: None,
    latex_link: Some("~\\cite{DBLP:conf/icpm/BurkeLW20}"),
    cli_command: None,
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Object(EbiObjectType::LabelledPetriNet)],
    ],
    input_names: &["FILE_1", "FILE_2"],
    input_helps: &[
        "A finite stochastic language (log) to get the occurrences from.",
        "A labelled Petri net with the control flow.",
    ],
    execute: |mut inputs, _| {
        let language = inputs
            .remove(0)
            .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let lpn = inputs.remove(0).to_type::<LabelledPetriNet>()?;
        Ok(EbiOutput::Object(EbiObject::StochasticLabelledPetriNet(
            lpn.mine_stochastic_alignment(language)?,
        )))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticLabelledPetriNet),
};

pub const EBI_DISCOVER_OCCURRENCE: EbiCommand = EbiCommand::Group {
    name_short: "occ",
    name_long: Some("occurrence"),
    explanation_short: "Give each transition a weight that matches the occurrences of its label; silent transitions get a weight of 1.",
    explanation_long: None,
    children: &[&EBI_DISCOVER_OCCURRENCE_LPN, &EBI_DISCOVER_OCCURRENCE_PTREE],
};

pub const EBI_DISCOVER_OCCURRENCE_LPN: EbiCommand = EbiCommand::Command {
    name_short: "lpn",
    name_long: Some("labelled-petri-net"),
    explanation_short: "Give each transition a weight that matches the occurrences of its label; silent transitions get a weight of 1.",
    explanation_long: None,
    latex_link: Some("~\\cite{DBLP:conf/icpm/BurkeLW20}"),
    cli_command: None,
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Object(EbiObjectType::LabelledPetriNet)],
    ],
    input_names: &["FILE_1", "FILE_2"],
    input_helps: &[
        "A finite stochastic language (log) to get the occurrences from.",
        "A labelled Petri net with the control flow.",
    ],
    execute: |mut inputs, _| {
        let language = inputs
            .remove(0)
            .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let lpn = inputs.remove(0).to_type::<LabelledPetriNet>()?;
        Ok(EbiOutput::Object(EbiObject::StochasticLabelledPetriNet(
            lpn.mine_occurrences_stochastic_lpn(language),
        )))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticLabelledPetriNet),
};

pub const EBI_DISCOVER_OCCURRENCE_PTREE: EbiCommand = EbiCommand::Command {
    name_short: "ptree",
    name_long: Some("process-tree"),
    explanation_short: "Give each leaf a weight that matches the occurrences of its label; silent leaves get a weight of 1.",
    explanation_long: None,
    latex_link: Some("~\\cite{DBLP:conf/icpm/BurkeLW20}"),
    cli_command: None,
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Object(EbiObjectType::ProcessTree)],
    ],
    input_names: &["LANG", "TREE"],
    input_helps: &[
        "A finite stochastic language (log) to get the occurrences from.",
        "A process tree with the control flow.",
    ],
    execute: |mut inputs, _| {
        let language = inputs
            .remove(0)
            .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let lpn = inputs.remove(0).to_type::<ProcessTree>()?;
        Ok(EbiOutput::Object(EbiObject::StochasticProcessTree(
            lpn.mine_occurrences_stochastic_tree(language),
        )))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticProcessTree),
};

pub const EBI_DISCOVER_UNIFORM: EbiCommand = EbiCommand::Group {
    name_short: "uni",
    name_long: Some("uniform"),
    explanation_short: "Give each transition a weight of 1.",
    explanation_long: None,
    children: &[&EBI_DISCOVER_UNIFORM_LPN, &EBI_DISCOVER_UNIFORM_PTREE],
};

pub const EBI_DISCOVER_UNIFORM_LPN: EbiCommand = EbiCommand::Command {
    name_short: "lpn",
    name_long: Some("labelled-petri-net"),
    explanation_short: "Give each transition a weight of 1 in a labelled Petri net.",
    explanation_long: None,
    latex_link: None,
    cli_command: None,
    exact_arithmetic: true,
    input_types: &[&[&EbiInputType::Object(EbiObjectType::LabelledPetriNet)]],
    input_names: &["MODEL"],
    input_helps: &["A labelled Petri net."],
    execute: |mut inputs, _| {
        let lpn = inputs.remove(0).to_type::<LabelledPetriNet>()?;
        Ok(EbiOutput::Object(EbiObject::StochasticLabelledPetriNet(
            lpn.mine_uniform_stochastic_lpn(),
        )))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticLabelledPetriNet),
};

pub const EBI_DISCOVER_UNIFORM_PTREE: EbiCommand = EbiCommand::Command {
    name_short: "ptree",
    name_long: Some("process-tree"),
    explanation_short: "Give each leaf a weight of 1 in a process tree.",
    explanation_long: None,
    latex_link: None,
    cli_command: None,
    exact_arithmetic: true,
    input_types: &[&[&EbiInputType::Object(EbiObjectType::ProcessTree)]],
    input_names: &["TREE"],
    input_helps: &["A process tree."],
    execute: |mut inputs, _| {
        let tree = inputs.remove(0).to_type::<ProcessTree>()?;
        Ok(EbiOutput::Object(EbiObject::StochasticProcessTree(
            tree.mine_uniform_stochastic_tree(),
        )))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticProcessTree),
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
        println!("log has {} traces", trace_num);
        let mut multiset = log.to_multiset();

        let mut multiset_map = HashMap::new();
        // for each key in multiset, divide the frequency by the trace number
        for (key, value) in multiset.iter_mut() {
            let mut trace_prob = fraction::Fraction::from(value.clone());
            trace_prob /= fraction::Fraction::from(trace_num);
            multiset_map.insert(key.clone(), trace_prob);
        }

        let activity_key = log.get_activity_key().clone();
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

        // let size_before_filtering = sts.get_size();
        println!("initial traces: {:?} and initial size: {}", sts.get_trace_number(), sts.get_size());
        // println!("before filtering, {:?} and {:?} and {:?} and {:?}", sts.get_sources(), sts.get_targets(), sts.get_activities(), sts.get_transition_frequencies());
        sts.filter_by_percentage(26.0).unwrap();
        // println!("after filtering, {:?} and {:?} and {:?} and {:?}", sts.get_sources(), sts.get_targets(), sts.get_activities(), sts.get_transition_frequencies());
        println!("after traces: {:?} and initial size: {}", sts.get_trace_number(), sts.get_size());

        sts.set_activity_key(stochastic_language.get_activity_key());
        let sdfa = convert_log_to_stochastic_transition_system(Box::new(stochastic_language), sts, target_size)?;
        Ok(EbiOutput::Object(EbiObject::StochasticDeterministicFiniteAutomaton(sdfa)))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticDeterministicFiniteAutomaton)
};


pub const EBI_DISCOVER_SDFA_PARALLEL: EbiCommand = EbiCommand::Command { 
    name_short: "disc_sdfa_parallel", 
    name_long: Some("discover_sdfa_parallel"), 
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
        let mut stochastic_language = FiniteStochasticLanguage::from((multiset_map.clone(), activity_key.clone()));
        log.translate_using_activity_key(stochastic_language.get_activity_key_mut());

        // get the target model size
        let target_size = *inputs.remove(0).to_type::<usize>()?;

         // Create a new transition system for each thread
        let mut sts = StochasticTransitionSystem::new();
        
        // Create transition system
        for trace_index in 0..log.len() {
            let trace = log.get_trace(trace_index).unwrap(); 
            let mut state = sts.get_and_add_initial_state();

            for activity in trace {
                state = sts.take_or_add_transition(state, *activity, 1);
            }
        }
        
        for i in 0..sts.get_sources().len() {
            // Add the source states
            let state = sts.get_sources()[i];
            let freq = sts.get_transition_frequencies()[i];
            sts.update_terminating_frequency(state, freq);
        }

        // let size_before_filtering = sts.get_size();
        let filter_percentages = vec![1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0];
        let sts_list: Vec<StochasticTransitionSystem> = sts.filter_by_multiple_percentages(&filter_percentages)?;
        // Process each threshold in parallel
        let results: Vec<_> = sts_list
            .into_par_iter()
            .map(|mut each_sts| {
                each_sts.set_activity_key(stochastic_language.get_activity_key());
                println!("start computation for sts of size: {} and target size: {}", each_sts.get_size(), target_size);
                
                // Convert to SDFA
                let stochastic_language_clone = FiniteStochasticLanguage::from((multiset_map.clone(), activity_key.clone()));
                convert_log_to_stochastic_transition_system(
                    Box::new(stochastic_language_clone), 
                    each_sts, 
                    target_size, 
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Randomly select one of the results
        let mut rng = rand::thread_rng();
        let random_index = rng.gen_range(0..results.len());
        let selected_sdfa = results.into_iter().nth(random_index).unwrap();

        Ok(EbiOutput::Object(EbiObject::StochasticDeterministicFiniteAutomaton(selected_sdfa)))
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::StochasticDeterministicFiniteAutomaton)
};