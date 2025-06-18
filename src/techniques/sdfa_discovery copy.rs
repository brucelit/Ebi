use std::{cmp::Ordering, collections::{BinaryHeap, HashMap, HashSet}, io::Cursor};
use anyhow::Result;
use crate::{ebi_objects::{stochastic_deterministic_finite_automaton::StochasticDeterministicFiniteAutomaton, stochastic_transition_system::StochasticTransitionSystem}, ebi_traits::{ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage, ebi_trait_queriable_stochastic_language::EbiTraitQueriableStochasticLanguage}, math::fraction::Fraction, techniques::unit_earth_movers_stochastic_conformance::UnitEarthMoversStochasticConformance};

use super::entropic_relevance::EntropicRelvance;



// Create a wrapper struct for our tuple
#[derive(Debug)]
struct TupleData {
    idx_to_sts: usize,
    red_state: usize,
    blue_state: usize,
    red_states: HashSet<usize>,
    blue_states: HashSet<usize>,
    heuristic: Fraction,
}

// Implement PartialEq manually
impl PartialEq for TupleData {
    fn eq(&self, other: &Self) -> bool {
        self.idx_to_sts == other.idx_to_sts && 
        self.red_state == other.red_state && 
        self.red_states == other.red_states &&
        self.blue_states == other.blue_states &&
        self.blue_state == other.blue_state && 
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
) -> Result<StochasticDeterministicFiniteAutomaton> {
    // merge the stochastic transition system to the target size
    Ok(merge_stochastic_transition_systems(stochastic_language_from_log, sts_from_log, target_size))
}


fn merge_stochastic_transition_systems(
    stochastic_language_from_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    sts_from_log: StochasticTransitionSystem,
    target_size: usize
) -> StochasticDeterministicFiniteAutomaton {
    // start from the initial stochastic transition system
    let mut curr_sts: StochasticTransitionSystem = sts_from_log.clone(); 
    println!("initial size: {:?}", curr_sts.get_size());
    println!("source states: {:?}", curr_sts.get_sources());
    println!("target states: {:?}", curr_sts.get_targets());
    println!("activities: {:?}", curr_sts.get_activities());    

    let new_sdfa = curr_sts.to_stochastic_deterministic_finite_automaton().unwrap();

    // get the initial entropic relevance and uemsc
    let uemsc = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
        Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    let er = stochastic_language_from_log.entropic_relevance(
        Box::new(new_sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    println!("initial uemsc: {:?}", uemsc);
    println!("initial entropic relevance: {:?}", er);

    // initialize the set of red states and blue states
    let mut red_states = HashSet::new();
    red_states.insert(sts_from_log.get_initial_state().clone());
    let mut blue_states = sts_from_log.get_blue_states(&red_states);
    println!("red states: {:?}", red_states);
    println!("blue states: {:?}", blue_states);


    // Create a binary heap to store tuples of (sts, red_state, blue_state, red_states, blue_states, heuristic)
    let mut heap: BinaryHeap<TupleData> = BinaryHeap::new();

    // map the sts to an index 
    let mut idx = 0;
    let mut idx_to_sts: HashMap<usize, StochasticTransitionSystem> = HashMap::new();
    idx_to_sts.insert(idx, curr_sts.clone());

    // helper variable to keep track of the current state
    let mut count: i32 =0;
    let mut curr_idx=0;

    // loop until the target size    
    while curr_sts.get_size() > target_size && count < 7
    {
        println!("\ncount: {} and len of heap: {}", count, heap.len());

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
                new_sts.trial_merge(
                    *each_red_state,
                    *each_blue_state,
                    &mut red_states.clone(),
                    &mut blue_states.clone());

                // get the heuristic
                let size_reduction = size - new_sts.get_size();

                if size_reduction<=1{
                    // we should not merge the states, add the blue state to the read states, and put the updated sts back to the heap
                    trial_merge_with_blue_state(
                        &mut new_sts,
                        size,
                        *each_red_state,
                        *each_blue_state,
                        &mut red_states.clone(),
                        &mut blue_states.clone(),
                    );
                }

                else{
                    let sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
                    let new_entropic_rel = stochastic_language_from_log.entropic_relevance(
                        Box::new(sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
                    ).unwrap();
                    let mut entropic_loss = new_entropic_rel.approximate().unwrap();
                    entropic_loss  -= entropic_rel.clone().approximate().unwrap();
                    let mut new_heuristic = Fraction::from(size_reduction);
                    new_heuristic /= entropic_loss.clone();

                    // update the red states and blue states after the trial merge
                    let updated_red_states = red_states.clone();
                    let updated_blue_states = new_sts.get_blue_states(&red_states);

                    // map the index to the sts
                    idx += 1;
                    let tuple = TupleData {
                        idx_to_sts: idx,
                        red_state: *each_red_state,
                        blue_state: *each_blue_state,
                        red_states: updated_red_states.clone(),
                        blue_states: updated_blue_states.clone(),
                        heuristic: new_heuristic,
                    };
                    idx_to_sts.insert(idx, new_sts.clone());
                    heap.push(tuple);
                }
            }
        }

            // get the stochastic transition tuple with the max heuristic
            if heap.len() == 0 {
                println!("The heap is empty");
            }
            // get the sts with the max heuristic
            let max_tuple = heap.pop().unwrap();

            let mut new_sts_with_max_heuristic = idx_to_sts.get(&max_tuple.idx_to_sts).unwrap().clone();
            println!("in loop, red state: {:?}, blue state: {:?}, heuristic: {:?}", max_tuple.red_state, max_tuple.blue_state, max_tuple.heuristic);
            println!("in loop, red states: {:?}, blue states: {:?}", max_tuple.red_states, max_tuple.blue_states);

            // perform the actual merge
            new_sts_with_max_heuristic.finish_merge(
                max_tuple.red_state,
                max_tuple.blue_state,
                &mut max_tuple.red_states.clone(),
                &mut max_tuple.blue_states.clone());

            // If the merge only reduce one state, give up the merge
            // if new_sts_with_max_heuristic.get_size() >= size-1  {
            //     // we should not merge the states, add the blue state to the read states, and put the updated sts back to the heap
            //     let mut temp_red_states = red_states.clone();
            //     temp_red_states.insert(max_tuple.blue_state);
            //     let temp_blue_states = new_sts_with_max_heuristic.get_blue_states(&temp_red_states);

            //     for each_red_state in temp_red_states.iter() {
            //         for each_blue_state in temp_blue_states.iter() {
            //             let mut new_sts: StochasticTransitionSystem = curr_sts.clone();
    
            //             // try to merge the blue into the red state
            //             new_sts.trial_merge(
            //                 *each_red_state,
            //                 *each_blue_state,
            //                 &mut temp_red_states.clone(),
            //                 &mut temp_blue_states.clone(),
            //             );
            
            //             // get the heuristic
            //             let size_reduction = size - new_sts.get_size();
            //             let sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
            
            //             let new_entropic_rel = stochastic_language_from_log.entropic_relevance(
            //                 Box::new(sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
            //             ).unwrap();
            //             let mut entropic_loss = new_entropic_rel.approximate().unwrap();
            //             entropic_loss  -= entropic_rel.clone().approximate().unwrap();
            //             let mut new_heuristic = Fraction::from(size_reduction);
            //             new_heuristic /= entropic_loss.clone();
            
            //             let sts_tuple = TupleData {
            //                 idx_to_sts: curr_idx,
            //                 red_state: *each_red_state,
            //                 blue_state: *each_blue_state,
            //                 red_states: temp_red_states.clone(),
            //                 blue_states: temp_blue_states.clone(),
            //                 heuristic: new_heuristic,
            //             };
            //             heap.push(sts_tuple);
            //         }
            //     }
            //     continue;
            // }
            
            // // the merge is valid, so we can update the current sts
            // else{
                //for the current sts, get the red states set and update the blue states set
            red_states = max_tuple.red_states.clone();
            blue_states = curr_sts.get_blue_states(&red_states);
            curr_sts = new_sts_with_max_heuristic.clone();
            curr_idx = max_tuple.idx_to_sts;
            count += 1;
    }

    let new_sts = curr_sts.clone();
    println!("source states: {:?}", new_sts.get_sources());
    println!("target states: {:?}", new_sts.get_targets());
    println!("activities: {:?}", new_sts.get_activities());        
    let new_sdfa = new_sts.to_stochastic_deterministic_finite_automaton().unwrap();
    let uemsc = stochastic_language_from_log.unit_earth_movers_stochastic_conformance(
        Box::new(new_sdfa.clone()) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    let er = stochastic_language_from_log.entropic_relevance(
        Box::new(new_sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
    ).unwrap();
    println!("entropic relevance: {:?}", er);
    println!("final uemsc: {:?}", uemsc);
    println!("final size: {:?}", curr_sts.get_size());
    println!("final states: {:?}", curr_sts.get_max_state());
    return curr_sts.to_stochastic_deterministic_finite_automaton().unwrap();
}


pub fn trial_merge_with_blue_state(
    sts: &mut StochasticTransitionSystem,
    size: usize,
    red_state: usize,
    blue_state: usize,
    red_states: &mut HashSet<usize>,
    blue_states: &mut HashSet<usize>,
) {
    let mut temp_red_states = red_states.clone();
    temp_red_states.insert(blue_state);
    let temp_blue_states = sts.get_blue_states(&temp_red_states);

    for each_red_state in temp_red_states.iter() {
        for each_blue_state in temp_blue_states.iter() {
            let mut new_sts: StochasticTransitionSystem = sts.clone();

            // try to merge the blue into the red state
            new_sts.trial_merge(
                *each_red_state,
                *each_blue_state,
                &mut temp_red_states.clone(),
                &mut temp_blue_states.clone(),
            );

            // get the heuristic
            let size_reduction = size - new_sts.get_size();

            if size_reduction<=1{
                // we should not merge the states, add the blue state to the read states, and put the updated sts back to the heap
                trial_merge_with_blue_state(
                    &mut new_sts,
                    size,
                    *each_red_state,
                    *each_blue_state,
                    &mut red_states.clone(),
                    &mut blue_states.clone(),
                );
            }

            let new_entropic_rel = stochastic_language_from_log.entropic_relevance(
                Box::new(sdfa) as Box<dyn EbiTraitQueriableStochasticLanguage>
            ).unwrap();
            let mut entropic_loss = new_entropic_rel.approximate().unwrap();
            entropic_loss  -= entropic_rel.clone().approximate().unwrap();
            let mut new_heuristic = Fraction::from(size_reduction);
            new_heuristic /= entropic_loss.clone();

            let sts_tuple = TupleData {
                idx_to_sts: curr_idx,
                red_state: *each_red_state,
                blue_state: *each_blue_state,
                red_states: temp_red_states.clone(),
                blue_states: temp_blue_states.clone(),
                heuristic: new_heuristic,
            };
            heap.push(sts_tuple);
        }
    }
}