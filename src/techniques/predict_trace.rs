use crate::{
    ebi_framework::displayable::Displayable,
    ebi_traits::{ebi_trait_stochastic_semantics::EbiTraitStochasticSemantics, ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage,},
    stochastic_semantics::stochastic_semantics::StochasticSemantics,
    techniques::astar_for_prediction,
    semantics::semantics::Semantics,
};
use anyhow::{Result, anyhow};
use ebi_arithmetic::{MaybeExact, Zero};
use ebi_objects::{Activity, ebi_objects::labelled_petri_net::TransitionIndex, ebi_objects::language_of_alignments::Move};
use std::ops::{Add, AddAssign};
use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

const ALPHA: f64 = 1.0;

#[derive(Clone, Debug)]
pub struct StochasticWeightedCost {
    cost: f64,
    probability: f64,
    stochastic_weighted_cost: f64,
}

impl Zero for StochasticWeightedCost {
    fn zero() -> Self {
        StochasticWeightedCost {
            cost: 0.0,
            probability: 1.0,
            stochastic_weighted_cost: 0.0,
        }
    }

    fn is_zero(&self) -> bool {
        self.stochastic_weighted_cost == 0.0
    }
}

impl Add for StochasticWeightedCost {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let probability = &self.probability * &other.probability;
        let cost:f64 = self.cost + other.cost;
        StochasticWeightedCost {
            cost: cost,
            probability: probability.clone(),
            stochastic_weighted_cost:  (1.0 + cost).log10().powf(ALPHA) * (1.0 - probability.log10()).powf(1.0-ALPHA),
        }
    }
}

impl AddAssign for StochasticWeightedCost {
    fn add_assign(&mut self, other: Self) {
        self.cost += other.cost;
        self.probability *= other.probability;
    }
}

impl Ord for StochasticWeightedCost {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let self_stochastic_cost = (1.0 + self.cost).log10().powf(ALPHA) * (1.0 - self.probability.log10()).powf(1.0-ALPHA);
        let other_stochastic_cost = (1.0 + other.cost).log10().powf(ALPHA) * (1.0 - other.probability.log10()).powf(1.0-ALPHA);
        self_stochastic_cost
            .partial_cmp(&other_stochastic_cost)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

impl PartialOrd for StochasticWeightedCost {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for StochasticWeightedCost {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for StochasticWeightedCost {}

pub trait PredictTrace {
    fn predict_trace(
        &self,
        trace: &Vec<Activity>,
        // balance: &Fraction,
    ) -> Result<u32>;

    fn predict_log(
        &self,
        event_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    ) -> Result<u32>;
}

impl PredictTrace for EbiTraitStochasticSemantics {
    fn predict_trace(
        &self,
        trace: &Vec<Activity>,
    ) -> Result<u32> {
        match self {
            EbiTraitStochasticSemantics::Usize(sem) => sem.predict_trace(trace),
            EbiTraitStochasticSemantics::Marking(sem) => sem.predict_trace(trace),
            EbiTraitStochasticSemantics::NodeStates(sem) => sem.predict_trace(trace),
        }
    }

    fn predict_log(
        &self,
        event_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    ) -> Result<u32> {
        match self {
            EbiTraitStochasticSemantics::Usize(sem) => sem.predict_log(event_log),
            EbiTraitStochasticSemantics::Marking(sem) => sem.predict_log(event_log),
            EbiTraitStochasticSemantics::NodeStates(sem) => sem.predict_log(event_log),
        }
    }
}

impl<State: Displayable> dyn StochasticSemantics<StoSemState = State, SemState = State, AliState = State> {
    pub fn predict_trace(
        &self,
        trace: &Vec<Activity>,
        // _balance: &Fraction,
    ) -> Result<u32> {
        // get the prefix and the next activity
        let mut prefix_collection = Vec::new();
        for i in 1..=trace.len() {
            prefix_collection.push((trace[..i].to_vec(), trace.get(i)));
            }
        let mut correct_predictions: i128 = 0;
            
        // iterate the prefix of the trace and align it
        for (prefix, next_activity) in prefix_collection {
            // get the start state
            if let Some(initial_state) = self.get_initial_state() {
                let start = (0, initial_state);

                // successor relation in the model
                let successors = |(trace_index, state): &(usize, State)| {
                    let mut result: Vec<((usize, State), StochasticWeightedCost, usize)> = vec![];

                    // log::debug!("successors of log {} model {}", trace_index, state);
                    if trace_index < &prefix.len() {
                        //we can do a log move
                        // log::debug!("\tlog move {}", trace[*trace_index]);

                        result.push((
                            (trace_index + 1, state.clone()),
                            StochasticWeightedCost {
                                cost: 1.0,
                                probability: 1.0,
                                stochastic_weighted_cost: 0.0
                            },
                            1,
                        ));
                    }

                    //walk through the enabled transitions in the model
                    for transition in self.get_enabled_transitions(&state) {
                        let total_weight = self
                            .get_total_weight_of_enabled_transitions(&state)
                            .unwrap();

                        let mut new_state = state.clone();
                        // log::debug!("\t\tnew state before {}", new_state);
                        let _ = self.execute_transition(&mut new_state, transition);
                        // log::debug!("\t\tnew state after {}", new_state);

                        let transition_weight = self.get_transition_weight(&state, transition);
                        let transition_probability: f64 = (transition_weight / &total_weight).approx().unwrap();

                        if let Some(activity) = self.get_transition_activity(transition) {
                            //non-silent model move
                            result.push((
                                (*trace_index, new_state.clone()),
                                StochasticWeightedCost {
                                    cost: 1.0,
                                    probability: transition_probability,
                                    stochastic_weighted_cost: 0.0
                                },
                                0,
                            ));
                            // log::debug!("\tmodel move t{} {} to {}", transition, activity, new_state);

                            //which may also be a synchronous move
                            if trace_index < &prefix.len() && activity == prefix[*trace_index] {
                                //synchronous move
                                // log::debug!("\tsynchronous move t{} {} to {}", transition, activity, new_state);
                                result.push((
                                    (trace_index + 1, new_state),
                                    StochasticWeightedCost {
                                        cost: 0.0,
                                        probability: transition_probability,
                                        stochastic_weighted_cost: 0.0
                                    },
                                    2,
                                ));
                            }
                        } else {
                            //silent move
                            result.push((
                                (*trace_index, new_state),
                                StochasticWeightedCost {
                                    cost: 0.0,
                                    probability: transition_probability,
                                    stochastic_weighted_cost: 0.0
                                },
                                2,
                            ));
                        }
                    }

                    // log::debug!("successors of {} {}: {:?}", trace_index, state, result);
                    result
                };

                //function that returns a heuristic on how far we are still minimally from a final state
                let heuristic = |_astate: &(usize, State)| StochasticWeightedCost::zero();

                //function that returns whether we are in a final synchronous product state
                let success = |(trace_index, state): &(usize, State)| {
                    trace_index == &prefix.len() && self.is_final_state(&state)
                };
                match astar_for_prediction::astar(&start, successors, heuristic, success) {
                    Some((path, _cost)) => {
                        let predicted_next_activity =  get_next_activity(self, &prefix, path);
                        if predicted_next_activity.is_none() && next_activity.is_none() {
                            println!("Correctly predicted next activity to be None after prefix {:?}", prefix);
                            correct_predictions += 1;
                            continue;
                        }
                        else if predicted_next_activity.is_none() && !next_activity.is_none() {
                            println!("Incorrectly predicted next activity to be None after prefix {:?}", prefix);
                            continue;
                        }
                        else if !predicted_next_activity.is_none() && next_activity.is_none() {
                            println!("Incorrectly predicted next activity after prefix {:?}, actual next activity is None", prefix);
                        }
                        else if predicted_next_activity.unwrap() == *next_activity.unwrap() {
                            correct_predictions += 1;
                            println!("Correctly predicted next activity {:?} after prefix {:?}", predicted_next_activity, prefix);
                        } else {
                            println!("Incorrectly predicted next activity {:?} to be after prefix {:?}, actual next activity is {:?}", predicted_next_activity, prefix, next_activity);
                        } 
                    }
                    None => {
                        println!("no alignment found for prefix {:?}", prefix);
                    },
                }
            } 
        }
        println!("Accuracy so far: {} / {}", correct_predictions, trace.len());
        Ok(correct_predictions as u32)
        // ---------------------------------------------------------------------
        // get the start state
        // if let Some(initial_state) = self.get_initial_state() {
        //     let start = (0, initial_state);

        //     // successor relation in the model
        //     let successors = |(trace_index, state): &(usize, State)| {
        //         let mut result: Vec<((usize, State), StochasticWeightedCost, usize)> = vec![];

        //         // log::debug!("successors of log {} model {}", trace_index, state);
        //         if trace_index < &trace.len() {
        //             //we can do a log move
        //             // log::debug!("\tlog move {}", trace[*trace_index]);

        //             result.push((
        //                 (trace_index + 1, state.clone()),
        //                 StochasticWeightedCost {
        //                     cost: 1.0,
        //                     probability: 1.0,
        //                     stochastic_weighted_cost: 0.0
        //                 },
        //                 1,
        //             ));
        //         }

        //         //walk through the enabled transitions in the model
        //         for transition in self.get_enabled_transitions(&state) {
        //             let total_weight = self
        //                 .get_total_weight_of_enabled_transitions(&state)
        //                 .unwrap();

        //             let mut new_state = state.clone();
        //             // log::debug!("\t\tnew state before {}", new_state);
        //             let _ = self.execute_transition(&mut new_state, transition);
        //             // log::debug!("\t\tnew state after {}", new_state);

        //             let transition_weight = self.get_transition_weight(&state, transition);
        //             let transition_probability: f64 = (transition_weight / &total_weight).approx().unwrap();

        //             if let Some(activity) = self.get_transition_activity(transition) {
        //                 //non-silent model move
        //                 result.push((
        //                     (*trace_index, new_state.clone()),
        //                     StochasticWeightedCost {
        //                         cost: 1.0,
        //                         probability: transition_probability,
        //                         stochastic_weighted_cost: 0.0
        //                     },
        //                     0,
        //                 ));
        //                 // log::debug!("\tmodel move t{} {} to {}", transition, activity, new_state);

        //                 //which may also be a synchronous move
        //                 if trace_index < &trace.len() && activity == trace[*trace_index] {
        //                     //synchronous move
        //                     // log::debug!("\tsynchronous move t{} {} to {}", transition, activity, new_state);
        //                     result.push((
        //                         (trace_index + 1, new_state),
        //                         StochasticWeightedCost {
        //                             cost: 0.0,
        //                             probability: transition_probability,
        //                             stochastic_weighted_cost: 0.0
        //                         },
        //                         2,
        //                     ));
        //                 }
        //             } else {
        //                 //silent move
        //                 result.push((
        //                     (*trace_index, new_state),
        //                     StochasticWeightedCost {
        //                         cost: 0.0,
        //                         probability: transition_probability,
        //                         stochastic_weighted_cost: 0.0
        //                     },
        //                     2,
        //                 ));
        //             }
        //         }

        //         // log::debug!("successors of {} {}: {:?}", trace_index, state, result);
        //         result
        //     };

        //     //function that returns a heuristic on how far we are still minimally from a final state
        //     let heuristic = |_astate: &(usize, State)| StochasticWeightedCost::zero();

        //     //function that returns whether we are in a final synchronous product state
        //     let success = |(trace_index, state): &(usize, State)| {
        //         trace_index == &trace.len() && self.is_final_state(&state)
        //     };


        //     match astar_for_prediction::astar(&start, successors, heuristic, success) {
        //         Some((path, _cost)) => {
        //             let moves = transform_alignment(self, &trace, path)?;
        //             let mut alignments = LanguageOfAlignments::new(self.activity_key().clone());
        //             alignments.push(moves);
        //             Ok(alignments)
        //         }
        //         None => Err(anyhow!("no alignment found")),
        //     }
        // } else {
        //     Err(anyhow!("Cannot align from an empty language."))
        // }
    }
}


impl<State: Displayable> dyn StochasticSemantics<StoSemState = State, SemState = State, AliState = State> {
    pub fn predict_log(
        &self,
        event_log: Box<dyn EbiTraitFiniteStochasticLanguage>,
    ) -> Result<u32> {
        let mut accuracy = 0.0;
        for (trace, probability) in event_log.iter_trace_probability() {
            // get the prefix and the next activity
            let mut prefix_collection = Vec::new();
            for i in 1..= trace.len() {
                prefix_collection.push((trace[..i].to_vec(), trace.get(i)));
                }
            let mut correct_predictions: i128 = 0;
                
            // iterate the prefix of the trace and align it
            for (prefix, next_activity) in prefix_collection {
                // get the start state
                if let Some(initial_state) = self.get_initial_state() {
                    let start = (0, initial_state);

                    // successor relation in the model
                    let successors = |(trace_index, state): &(usize, State)| {
                        let mut result: Vec<((usize, State), StochasticWeightedCost, usize)> = vec![];

                        // log::debug!("successors of log {} model {}", trace_index, state);
                        if trace_index < &prefix.len() {
                            //we can do a log move
                            // log::debug!("\tlog move {}", trace[*trace_index]);

                            result.push((
                                (trace_index + 1, state.clone()),
                                StochasticWeightedCost {
                                    cost: 1.0,
                                    probability: 1.0,
                                    stochastic_weighted_cost: 0.0
                                },
                                1,
                            ));
                        }

                        //walk through the enabled transitions in the model
                        for transition in self.get_enabled_transitions(&state) {
                            let total_weight = self
                                .get_total_weight_of_enabled_transitions(&state)
                                .unwrap();

                            let mut new_state = state.clone();
                            // log::debug!("\t\tnew state before {}", new_state);
                            let _ = self.execute_transition(&mut new_state, transition);
                            // log::debug!("\t\tnew state after {}", new_state);

                            let transition_weight = self.get_transition_weight(&state, transition);
                            let transition_probability: f64 = (transition_weight / &total_weight).approx().unwrap();

                            if let Some(activity) = self.get_transition_activity(transition) {
                                //non-silent model move
                                result.push((
                                    (*trace_index, new_state.clone()),
                                    StochasticWeightedCost {
                                        cost: 1.0,
                                        probability: transition_probability,
                                        stochastic_weighted_cost: 0.0
                                    },
                                    0,
                                ));
                                // log::debug!("\tmodel move t{} {} to {}", transition, activity, new_state);

                                //which may also be a synchronous move
                                if trace_index < &prefix.len() && activity == prefix[*trace_index] {
                                    //synchronous move
                                    // log::debug!("\tsynchronous move t{} {} to {}", transition, activity, new_state);
                                    result.push((
                                        (trace_index + 1, new_state),
                                        StochasticWeightedCost {
                                            cost: 0.0,
                                            probability: transition_probability,
                                            stochastic_weighted_cost: 0.0
                                        },
                                        2,
                                    ));
                                }
                            } else {
                                //silent move
                                result.push((
                                    (*trace_index, new_state),
                                    StochasticWeightedCost {
                                        cost: 0.0,
                                        probability: transition_probability,
                                        stochastic_weighted_cost: 0.0
                                    },
                                    2,
                                ));
                            }
                        }

                        // log::debug!("successors of {} {}: {:?}", trace_index, state, result);
                        result
                    };

                    //function that returns a heuristic on how far we are still minimally from a final state
                    let heuristic = |_astate: &(usize, State)| StochasticWeightedCost::zero();

                    //function that returns whether we are in a final synchronous product state
                    let success = |(trace_index, state): &(usize, State)| {
                        trace_index == &prefix.len() && self.is_final_state(&state)
                    };
                    match astar_for_prediction::astar(&start, successors, heuristic, success) {
                        Some((path, _cost)) => {
                            let predicted_next_activity =  get_next_activity(self, &prefix, path);
                            if predicted_next_activity.is_none() && next_activity.is_none() {
                                println!("Correctly predicted next activity to be None after prefix {:?}", prefix);
                                correct_predictions += 1;
                                continue;
                            }
                            else if predicted_next_activity.is_none() && !next_activity.is_none() {
                                println!("Incorrectly predicted next activity to be None after prefix {:?}", prefix);
                                continue;
                            }
                            else if !predicted_next_activity.is_none() && next_activity.is_none() {
                                println!("Incorrectly predicted next activity after prefix {:?}, actual next activity is None", prefix);
                            }
                            else if predicted_next_activity.unwrap() == *next_activity.unwrap() {
                                correct_predictions += 1;
                                println!("Correctly predicted next activity: {:?} after prefix {:?}", predicted_next_activity, prefix);
                            } else {
                                println!("Incorrectly predicted next activity {:?} to be after prefix {:?}, actual next activity is {:?}", predicted_next_activity, prefix, next_activity);
                            } 
                        }
                        None => {
                            println!("no alignment found for prefix {:?}", prefix);
                        },
                    }
                } 
            }
            accuracy += correct_predictions as f64 / trace.len() as f64 * probability.clone().approx().unwrap() as f64;
        }
        println!("accuracy for prediction: {}", accuracy);
        Ok(1 as u32)

    }
}

pub fn get_next_activity<T, State>(
    semantics: &T,
    prefix: &Vec<Activity>,
    states: Vec<(usize, State)>,
) -> Option<Activity>
where
    T: Semantics<SemState = State> + ?Sized,
    State: Display + Debug + Clone + Hash + Eq,
{
    // log::debug!("transform alignment {:?}", states);
    let mut alignment = vec![];

    let mut it = states.into_iter();

    let (mut previous_trace_index, mut previous_state) = it.next().unwrap();

    let mut counter = 0;
    let mut find_the_next = false;
    let mut next_activity= None;
    for (trace_index, state) in it {
        // log::debug!("transform {} from {} to {}", trace_index, previous_state, state);
        if trace_index != previous_trace_index {
            //we did a move on the log
            let activity = prefix[previous_trace_index];

            if previous_state == state {
                //we did not move on the model => log move
                alignment.push(Move::LogMove(activity));
                counter += 1;
                if counter == prefix.len() {
                    find_the_next = true;
                }

            } else {
                //we moved on the model => synchronous move
                // let transition =
                //     find_transition_with_label(semantics, &previous_state, &state, activity)
                //         .with_context(|| {
                //             format!(
                //                 "Map synchronous move from {} {} to {} {} with label {}",
                //                 previous_trace_index, previous_state, trace_index, state, activity
                //             )
                //         })?;
                // alignment.push(Move::SynchronousMove(activity, transition));
                if find_the_next {
                    next_activity = Some(activity);
                    break;
                }
                counter += 1;
                if counter == prefix.len() {
                    find_the_next = true;
                }
            }
        } else {
            //we did not do a move on the log

            if let Some(transition) =
                is_there_a_silent_transition_enabled(semantics, &previous_state, &state)
            {
                //there is a silent transition enabled, which is the cheapest
                alignment.push(Move::SilentMove(transition));
            } else {
                //otherwise, we take an arbitrary labelled model move
                let transition = find_labelled_transition(semantics, &previous_state, &state);
                // alignment.push(Move::ModelMove(
                //     semantics.get_transition_activity(transition).unwrap(),
                //     transition,
                // ));
                if find_the_next {
                    next_activity = semantics.get_transition_activity(transition.unwrap());
                    break;
                }
            }
        }

        previous_trace_index = trace_index;
        previous_state = state;

        // log::debug!("prefix: {:?}", alignment);
    }
    next_activity
}

pub fn find_transition_with_label<T, FS>(
    semantics: &T,
    from: &FS,
    to: &FS,
    label: Activity,
) -> Result<TransitionIndex>
where
    T: Semantics<SemState = FS> + ?Sized,
    FS: Display + Debug + Clone + Hash + Eq,
{
    // log::debug!("find transition with label {}", label);
    for transition in semantics.get_enabled_transitions(from) {
        // log::debug!("transition {} is enabled", transition);
        if semantics.get_transition_activity(transition) == Some(label) {
            let mut from = from.clone();
            semantics.execute_transition(&mut from, transition)?;
            if &from == to {
                return Ok(transition);
            }
        }
    }
    Err(anyhow!(
        "There is no transition with activity {} that brings the model from {} to {}",
        label,
        from,
        to
    ))
}


pub fn is_there_a_silent_transition_enabled<T, FS>(
    semantics: &T,
    from: &FS,
    to: &FS,
) -> Option<TransitionIndex>
where
    T: Semantics<SemState = FS> + ?Sized,
    FS: Display + Debug + Clone + Hash + Eq,
{
    // log::debug!("is there a silent transition enabled from {} to {}", from, to);
    // log::debug!("enabled transitions {:?}", semantics.get_enabled_transitions(from));
    for transition in semantics.get_enabled_transitions(from) {
        if semantics.is_transition_silent(transition) {
            let mut from = from.clone();
            let _ = semantics.execute_transition(&mut from, transition);
            if &from == to {
                // log::debug!("yes");
                return Some(transition);
            }
        }
    }
    // log::debug!("no");
    None
}


pub fn find_labelled_transition<T, FS>(semantics: &T, from: &FS, to: &FS) -> Result<TransitionIndex>
where
    T: Semantics<SemState = FS> + ?Sized,
    FS: Display + Debug + Clone + Hash + Eq,
{
    for transition in semantics.get_enabled_transitions(from) {
        if !semantics.is_transition_silent(transition) {
            let mut from = from.clone();
            semantics.execute_transition(&mut from, transition)?;
            if &from == to {
                return Ok(transition);
            }
        }
    }
    Err(anyhow!(
        "There is no transition with any activity enabled that brings the model from {} to {}",
        from,
        to
    ))
}