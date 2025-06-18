use std::{cmp::Ordering, collections::{BinaryHeap, HashMap, HashSet}, io::Cursor};
use anyhow::Result;
use crate::{ebi_objects::{stochastic_deterministic_finite_automaton::StochasticDeterministicFiniteAutomaton, stochastic_transition_system::StochasticTransitionSystem}, ebi_traits::{ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage, ebi_trait_queriable_stochastic_language::EbiTraitQueriableStochasticLanguage}, math::{fraction::Fraction, traits::Zero}, techniques::unit_earth_movers_stochastic_conformance::UnitEarthMoversStochasticConformance};

use super::entropic_relevance::EntropicRelvance;



// Create a wrapper struct for our tuple
#[derive(Debug)]
struct TupleData {
    idx_to_sts: usize,
    red_states: HashSet<usize>,
    heuristic: Fraction,
}

// Implement PartialEq manually
impl PartialEq for TupleData {
    fn eq(&self, other: &Self) -> bool {
        self.idx_to_sts == other.idx_to_sts && 
        self.red_states == other.red_states &&
        self.heuristic == other.heuristic
    }
}

// For BinaryHeap, we'll focus on making a valid implementation that works for our needs
impl Eq for TupleData {}

// The key part: implement ordering based on the f64 priority value
impl Ord for TupleData {
    fn cmp(&self, other: &Self) -> Ordering {
        // Use partial_cmp for f64 values, defaulting to Equal in case of NaN
        // This reverses the order to make it a max-heap by priority
        self.heuristic.partial_cmp(&other.heuristic)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for TupleData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


pub fn convert_log_to_stochastic_transition_system(
    stochastic_language_from_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    sts_from_log: StochasticTransitionSystem,
    target_size: usize,
    size_before_filtering: usize
) -> Result<StochasticDeterministicFiniteAutomaton> {
    // merge the stochastic transition system to the target size
    Ok(merge_stochastic_transition_systems(stochastic_language_from_log, sts_from_log, target_size,size_before_filtering))
}


fn merge_stochastic_transition_systems(
    stochastic_language_from_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    sts_from_log: StochasticTransitionSystem,
    target_size: usize,
    size_before_filtering: usize
) -> StochasticDeterministicFiniteAutomaton {
    // start from the initial stochastic transition system
    let mut curr_sts: StochasticTransitionSystem = sts_from_log.clone(); 
    let new_sdfa = curr_sts.to_stochastic_deterministic_finite_automaton().unwrap();

    // get the initial entropic relevance and uemsc
    let init_uemsc = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
        Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    let init_er = stochastic_language_from_log.entropic_relevance(
        Box::new(new_sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    // println!("initial source states: {:?}", curr_sts.get_sources());
    // println!("initial target states: {:?}", curr_sts.get_targets());
    // println!("initial activities: {:?}", curr_sts.get_activities());
    // println!("initial frequencies: {:?}", curr_sts.get_transition_frequencies());
    // println!("initial size: {:?}", curr_sts.get_size());
    let init_size = curr_sts.get_size();
    // Create a binary heap to store tuples of (sts, red_state, blue_state, red_states, blue_states, heuristic)
    let mut heap: BinaryHeap<TupleData> = BinaryHeap::new();

    // map the sts to an index 
    let mut idx = 0;
    let mut idx_to_sts: HashMap<usize, StochasticTransitionSystem> = HashMap::new();
    idx_to_sts.insert(idx, curr_sts.clone());

    // helper variable to keep track of the current state
    let mut count: i32 =0;
    let mut curr_idx=0;
    let mut red_states = HashSet::new();
    red_states.insert(curr_sts.get_initial_state().clone());
    let mut blue_states = curr_sts.get_blue_states(&red_states);
    // loop until the target size    

    let init_size = curr_sts.get_size();

    let mut another_size = curr_sts.get_size();
    // while curr_sts.get_size() > 100
    while another_size > target_size
    {   
         // initialize the set of red states and blue states

        // println!("red states: {:?}", red_states);
        // println!("blue states: {:?}", blue_states);

        // get the current er and model size before any trial merging
        let sdfa_before: StochasticDeterministicFiniteAutomaton = curr_sts.to_stochastic_deterministic_finite_automaton().unwrap();

        let entropic_rel = stochastic_language_from_log.entropic_relevance(
            Box::new(sdfa_before) as Box<dyn EbiTraitQueriableStochasticLanguage>
        ).unwrap();
        let size = curr_sts.get_size();
        
        // trial merge each pair of red and blue states
        for each_red_state in red_states.iter() {
            for each_blue_state in blue_states.iter() {
                let mut new_sts: StochasticTransitionSystem = curr_sts.clone();
                // try to merge the blue into the red state
                let merge_result = new_sts.trial_merge(
                    *each_red_state,
                    *each_blue_state,
                    &mut red_states.clone(),
                    &mut blue_states.clone());
                
                match merge_result {
                    Ok(()) => {
                        // println!("merged red state: {:?} with blue state: {:?}", each_red_state, each_blue_state);

                        // println!("source states: {:?}", new_sts.get_sources());
                        // println!("target states: {:?}", new_sts.get_targets());
                        let size_reduction = size - new_sts.get_size();
                        let sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
                        let new_entropic_rel = stochastic_language_from_log.entropic_relevance(
                            Box::new(sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
                        ).unwrap();
                        let mut entropic_loss = new_entropic_rel.approximate().unwrap();
                        entropic_loss  -= entropic_rel.clone().approximate().unwrap();
                        let mut new_heuristic = Fraction::from(size_reduction);
                        new_heuristic /= entropic_loss.clone();

                        // we should not merge the states, add the blue 
                        if curr_sts.get_size() - new_sts.get_size()<=1 || new_heuristic <= Fraction::zero(){
                            // println!("skipping merge");
                            // update the red states and blue states after the trial merge
                            let mut updated_red_states = red_states.clone();
                            updated_red_states.insert(*each_blue_state);
                            let updated_blue_states = new_sts.get_blue_states(&updated_red_states);

                            // map the index to the sts
                            idx += 1;
                            let tuple = TupleData {
                                idx_to_sts: idx,
                                red_states: updated_red_states.clone(),
                                heuristic: Fraction::zero(),
                            };
                            idx_to_sts.insert(idx, curr_sts.clone());
                            heap.push(tuple);
                        }
                        else{
                            // update the red states and blue states after the trial merge
                            let updated_red_states = red_states.clone();

                            // map the index to the sts
                            idx += 1;
                            let tuple = TupleData {
                                idx_to_sts: idx,
                                red_states: updated_red_states.clone(),
                                heuristic: new_heuristic.clone(),
                            };
                            idx_to_sts.insert(idx, new_sts.clone());
                            heap.push(tuple);
                            // println!("pushing to heap");

                        }
                    },
                        Err(e) => {
                            continue; // skip this iteration if the merge fails
                        }
                    }
            }
        }
        // get the stochastic transition tuple with the max heuristic
        if heap.len() == 0 {
            println!("The heap is empty");
            break; // exit the loop if there are no more tuples to process
        }
        // get the sts with the max heuristic
        let max_tuple = heap.pop().unwrap();

        let new_sts_with_max_heuristic = idx_to_sts.get(&max_tuple.idx_to_sts).unwrap().clone();

        let new_sdfa = new_sts_with_max_heuristic.to_stochastic_deterministic_finite_automaton().unwrap();
        let uemsc = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
            Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
        ).unwrap();
        println!("\ncount:{:?}, heap len: {}, max heuristic: {:.3}, uEMSC: {:.3}, size reduction: {}", count, heap.len(), max_tuple.heuristic, uemsc, init_size - new_sts_with_max_heuristic.get_size());
        red_states = max_tuple.red_states.clone();
        blue_states = curr_sts.get_blue_states(&red_states);


        // println!("red states: {:?}", red_states);
        // println!("blue states: {:?}", blue_states);
        // println!("source states: {:?}", new_sts_with_max_heuristic.get_sources());
        // println!("target states: {:?}", new_sts_with_max_heuristic.get_targets());
        curr_sts = new_sts_with_max_heuristic.clone();
        let mut another_sts = curr_sts.clone();
        another_sts.minimize_hopcroft().unwrap();
        another_size = another_sts.get_size();
        count += 1;
    }

    let mut new_sts = curr_sts.clone();   
    let new_sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
    let uemsc_before_merging = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
        Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    let er_before_merging = stochastic_language_from_log.entropic_relevance(
        Box::new(new_sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    println!("start merging equivalent states");
    let size_before_merging = new_sts.get_size();
    new_sts.minimize_hopcroft().unwrap();
    println!("finished merging equivalent states");

    let new_sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
    // println!("after minimization, source states: {:?}", new_sdfa.get_sources());
    // println!("after minimization, target states: {:?}", new_sdfa.get_targets());
    // println!("after minimization, activities: {:?}", new_sdfa.get_activities());
    let uemsc = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
        Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    let er = stochastic_language_from_log.entropic_relevance(
        Box::new(new_sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    println!("num of merging: {:?}", count);
    println!("initial uemsc:{:?}", init_uemsc);
    println!("final uemsc before merging: {:?}", uemsc_before_merging);
    println!("final uemsc: {:?}", uemsc);
    println!("initial entropic relevance: {:?}", init_er);
    println!("final entropic relevance before merging: {:?}", er_before_merging);

    println!("final entropic relevance: {:?}", er);
    println!("initial size: {:?}", size_before_filtering);
    println!("size before collapsing states: {:?}", size_before_merging);
    println!("final size: {:?}", new_sts.get_size());
    return new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
}