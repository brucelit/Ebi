use crate::{
    ebi_framework::displayable::Displayable,
    ebi_traits::ebi_trait_stochastic_semantics::EbiTraitStochasticSemantics,
    stochastic_semantics::stochastic_semantics::StochasticSemantics,
    techniques::{align::transform_alignment, astar_for_prediction},
};
use anyhow::{Result, anyhow};
use ebi_arithmetic::{MaybeExact, Zero};
use ebi_objects::{Activity, LanguageOfAlignments};
use std::ops::{Add, AddAssign};

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
    ) -> Result<LanguageOfAlignments>;
}

impl PredictTrace for EbiTraitStochasticSemantics {
    fn predict_trace(
        &self,
        trace: &Vec<Activity>,
    ) -> Result<LanguageOfAlignments> {
        match self {
            EbiTraitStochasticSemantics::Usize(sem) => sem.predict_trace(trace),
            EbiTraitStochasticSemantics::Marking(sem) => sem.predict_trace(trace),
            EbiTraitStochasticSemantics::NodeStates(sem) => sem.predict_trace(trace),
        }
    }
}

impl<State: Displayable> dyn StochasticSemantics<StoSemState = State, SemState = State, AliState = State> {
    pub fn predict_trace(
        &self,
        trace: &Vec<Activity>,
        // _balance: &Fraction,
    ) -> Result<LanguageOfAlignments> {
        let mut prefix_collection = Vec::new();
        for i in 1..=trace.len() {
                prefix_collection.push(trace[..i].to_vec());
            }
            
        // iterate the prefix of the trace and align it
        for prefix in prefix_collection {
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
                        let moves = transform_alignment(self, &prefix, path)?;
                        let mut alignments = LanguageOfAlignments::new(self.activity_key().clone());
                        alignments.push(moves);
                        println!("Predicted next activity after prefix {:?}: {:?}", prefix, alignments.alignments);
                        // return the non-silent model move or synchronous move from alignment after the last log move
                        let last_trace = prefix.
                    }
                    None => {
                        println!("no alignment found for prefix {:?}", prefix);
                    },
                }
            } 
        }

        // ---------------------------------------------------------------------
            // get the start state
            if let Some(initial_state) = self.get_initial_state() {
                let start = (0, initial_state);

                // successor relation in the model
                let successors = |(trace_index, state): &(usize, State)| {
                    let mut result: Vec<((usize, State), StochasticWeightedCost, usize)> = vec![];

                    // log::debug!("successors of log {} model {}", trace_index, state);
                    if trace_index < &trace.len() {
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
                            if trace_index < &trace.len() && activity == trace[*trace_index] {
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
                    trace_index == &trace.len() && self.is_final_state(&state)
                };


                match astar_for_prediction::astar(&start, successors, heuristic, success) {
                    Some((path, _cost)) => {
                        let moves = transform_alignment(self, &trace, path)?;
                        let mut alignments = LanguageOfAlignments::new(self.activity_key().clone());
                        alignments.push(moves);
                        Ok(alignments)
                    }
                    None => Err(anyhow!("no alignment found")),
                }
            } else {
                Err(anyhow!("Cannot align from an empty language."))
            }
        }
}
