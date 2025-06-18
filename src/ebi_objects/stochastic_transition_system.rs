use crate::{
    ebi_framework::{
        activity_key::{Activity, ActivityKey}, infoable::Infoable
    }, ebi_objects::stochastic_deterministic_finite_automaton::StochasticDeterministicFiniteAutomaton, ebi_traits::ebi_trait_graphable::EbiTraitGraphable, math::{fraction::Fraction, traits::Signed}
    
};
use anyhow::{Result, anyhow};
use layout::topo::layout::VisualGraph;
use num::Float;
use std::{
    cmp::{max, Ordering},
    collections::{HashMap, HashSet, VecDeque},
    fmt
};

pub const FORMAT_SPECIFICATION: &str = "A stochastic transition system is a JSON structure with the top level being an object.
    This object contains the following key-value pairs:
    \\begin{itemize}
    \\item \\texttt{initialState} being the index of the initial state. This field is optional: if omitted, the stochastic transition system has an empty stochastic language.
    \\item \\texttt{transitions} being a list of transitions. 
    Each transition is an object with \\texttt{from} being the source state index of the transition, 
    \\texttt{to} being the target state index of the transition, 
    \\texttt{label} being the activity of the transition, and
    \\texttt{freq} being the frequency of the transition (may be given as a fraction in a string or a float value. Must be $\\leq 1$). 
    Silent transitions are not supported.
    The file format supports deadlocks and livelocks.
    The frequency that a trace terminates in a state is the sum of frequency of the incoming transitions - the sum probability of the outgoing transitions of the state.
    \\end{itemize}";

#[derive(Debug)]
pub struct StochasticTransitionSystem {
    activity_key: ActivityKey,
    initial_state: usize,
    max_state: usize,
    sources: Vec<usize>, //transition -> source of arc
    targets: Vec<usize>, //transition -> target of arc
    activities: Vec<Activity>, //transition -> activity of arc (every transition is labelled)
    transition_frequencies: Vec<usize>, //transition -> frequency of arc
    state_terminating_frequencies: Vec<usize>, //state -> termination frequency
}

impl StochasticTransitionSystem {

    pub fn new() -> Self {
        Self {
            activity_key: ActivityKey::new(),
            max_state: 0,
            initial_state: 0,
            sources: vec![],
            targets: vec![],
            activities: vec![],
            transition_frequencies: vec![],
            state_terminating_frequencies: vec![0],
        }
    }


/// Merges a blue state into a red state following the merging rules:
/// 1. The transition from red to blue becomes a self-loop on red
/// 2. All outgoing transitions from blue are inherited by red
/// 3. If red already has an outgoing transition with the same label as blue's transition,
///    the target states are merged and frequencies are combined
/// 4. Updates the red_states and blue_states sets accordingly
pub fn trial_merge(
    &mut self,
    red_state: usize,
    blue_state: usize,
    red_states: &mut HashSet<usize>,
    blue_states: &mut HashSet<usize>
) -> Result<()> {
    // println!("Merging red state {} with blue state {}", red_state, blue_state);
    // println!("Before merge - sources: {:?}", self.sources);
    // println!("Before merge - targets: {:?}", self.targets);
    
    // Find the transition from red_state to blue_state
    let mut red_to_blue_transition = None;
    for i in 0..self.sources.len() {
        if self.sources[i] == red_state && self.targets[i] == blue_state {
            red_to_blue_transition = Some(i);
            break;
        }
    }
    
    let transition_idx = red_to_blue_transition
        .ok_or_else(|| anyhow!("No transition found from red state {} to blue state {}", red_state, blue_state))?;
    
    // Step 1: Convert red->blue transition to self-loop on red
    self.targets[transition_idx] = red_state;
    // println!("Step 1: Converted transition {} to self-loop", transition_idx);
    
    // Step 2: Collect all outgoing transitions from blue_state
    let mut blue_outgoing = Vec::new();
    for i in 0..self.sources.len() {
        if self.sources[i] == blue_state {
            blue_outgoing.push((self.activities[i], self.targets[i], self.transition_frequencies[i]));
        }
    }
    
    // println!("Blue state {} outgoing transitions: {:?}", blue_state, blue_outgoing);
    
    // Step 3: Remove all transitions from blue_state first
    self.remove_transitions_from_state(blue_state);
    
    // Step 4: Process each collected outgoing transition
    for (blue_activity, blue_target, blue_frequency) in blue_outgoing {
        // println!("Processing blue transition: {}->{}({:?}, freq={})", blue_state, blue_target, blue_activity, blue_frequency);
        
        // Check if red_state already has an outgoing transition with the same activity
        let mut red_existing = None;
        for j in 0..self.sources.len() {
            if self.sources[j] == red_state && self.activities[j] == blue_activity {
                red_existing = Some(j);
                break;
            }
        }
        
        if let Some(red_idx) = red_existing {
            let red_target = self.targets[red_idx];
            // println!("Found conflicting transition: red {}->{}({:?}) vs blue {}->{}({:?})", 
            //         red_state, red_target, blue_activity, blue_state, blue_target, blue_activity);
            
            // Combine frequencies
            self.transition_frequencies[red_idx] += blue_frequency;
            
            // If targets are different, merge them
            if blue_target != red_target {
                // println!("Merging target states {} and {}", blue_target, red_target);
                self.merge_states(blue_target, red_target, red_states, blue_states)?;
            }
        } else {
            // No conflict, add new transition from red to blue_target
            // println!("Adding new transition {}->{}({:?}) with freq {}", 
            //         red_state, blue_target, blue_activity, blue_frequency);
            
            // Find the correct position to insert (maintain sorted order if needed)
            let insert_pos = self.sources.len(); // For now, just append
            
            self.sources.insert(insert_pos, red_state);
            self.targets.insert(insert_pos, blue_target);
            self.activities.insert(insert_pos, blue_activity);
            self.transition_frequencies.insert(insert_pos, blue_frequency);
        }
    }
    
    // Step 5: Transfer terminating frequency
    self.state_terminating_frequencies[red_state] += self.state_terminating_frequencies[blue_state];
    self.state_terminating_frequencies[blue_state] = 0;
    
    // Step 6: Update state sets
    blue_states.remove(&blue_state);
    
    // println!("After merge - sources: {:?}", self.sources);
    // println!("After merge - targets: {:?}", self.targets);
    
    Ok(())
}

/// Helper function to remove all transitions from a given state
fn remove_transitions_from_state(&mut self, state: usize) {
    let mut i = 0;
    while i < self.sources.len() {
        if self.sources[i] == state {
            self.sources.remove(i);
            self.targets.remove(i);
            self.activities.remove(i);
            self.transition_frequencies.remove(i);
            // Don't increment i since we removed an element
        } else {
            i += 1;
        }
    }
}

/// Helper function to merge two states by redirecting all references from source_state to target_state
fn merge_states(
    &mut self,
    source_state: usize,
    target_state: usize,
    red_states: &mut HashSet<usize>,
    blue_states: &mut HashSet<usize>
) -> Result<()> {
    // println!("  Merging states: {} -> {}", source_state, target_state);
    
    if source_state == target_state {
        return Ok(()); // Nothing to merge
    }
    
    // Redirect all transitions pointing TO source_state to point to target_state
    for i in 0..self.targets.len() {
        if self.targets[i] == source_state {
            // println!("    Redirecting incoming transition to {} -> {}", source_state, target_state);
            self.targets[i] = target_state;
        }
    }
    
    // Collect outgoing transitions from source_state
    let mut source_outgoing = Vec::new();
    for i in 0..self.sources.len() {
        if self.sources[i] == source_state {
            source_outgoing.push((self.activities[i], self.targets[i], self.transition_frequencies[i]));
        }
    }
    
    // Remove all transitions from source_state
    self.remove_transitions_from_state(source_state);
    
    // Process each outgoing transition
    for (activity, dest, frequency) in source_outgoing {
        // println!("    Processing outgoing from {}: {}->{}({:?}, freq={})", source_state, source_state, dest, activity, frequency);
        
        // Check if target_state already has a transition with the same activity
        let mut target_existing = None;
        for j in 0..self.sources.len() {
            if self.sources[j] == target_state && self.activities[j] == activity {
                target_existing = Some(j);
                break;
            }
        }
        
        if let Some(target_idx) = target_existing {
            let existing_dest = self.targets[target_idx];
            // println!("    Found existing transition from {}: {}->{}({:?})", target_state, target_state, existing_dest, activity);
            
            // Merge with existing transition
            self.transition_frequencies[target_idx] += frequency;
            
            // If destinations are different, merge them too
            if existing_dest != dest {
                // println!("    Recursively merging destinations {} and {}", dest, existing_dest);
                self.merge_states(dest, existing_dest, red_states, blue_states)?;
            }
        } else {
            // Add new transition from target_state to dest
            // println!("    Adding new transition {}->{}({:?}) with freq {}", target_state, dest, activity, frequency);
            
            let insert_pos = self.sources.len();
            self.sources.insert(insert_pos, target_state);
            self.targets.insert(insert_pos, dest);
            self.activities.insert(insert_pos, activity);
            self.transition_frequencies.insert(insert_pos, frequency);
        }
    }
    
    // Transfer terminating frequency
    self.state_terminating_frequencies[target_state] += self.state_terminating_frequencies[source_state];
    self.state_terminating_frequencies[source_state] = 0;
    
    // Update state sets
    blue_states.remove(&source_state);
    red_states.remove(&source_state);
    
    // Update initial state if necessary
    if self.initial_state == source_state {
        self.initial_state = target_state;
    }
    
    // println!("  Completed merging {} -> {}", source_state, target_state);
    Ok(())
}

    // Find all blue states that are direct successors of red states
    pub fn get_blue_states(&self, red_states: &HashSet<usize>) -> HashSet<usize> {
        let mut blue_states = HashSet::new();
        
        // Find all states that are direct successors of red states
        for i in 0..self.sources.len() {
            if red_states.contains(&self.sources[i]) {
                let target = self.targets[i];
                if !red_states.contains(&target) {
                    blue_states.insert(target);
                }
            }
        }
        
        blue_states
    }

    pub fn to_stochastic_deterministic_finite_automaton(&self) -> Result<StochasticDeterministicFiniteAutomaton> {
        // Create a new SDFA with the same activity key
        let mut sdfa = StochasticDeterministicFiniteAutomaton::new();
        sdfa.set_activity_key(&self.activity_key);
        sdfa.set_initial_state(Some(self.initial_state));
        
        // Calculate total outgoing frequency for each state
        // This includes both outgoing transitions and termination frequencies
        let mut total_frequencies = vec![0usize; self.max_state + 1];
        let mut state_map = HashMap::new();

        // Map states from sts to sdfa
        let mut source_idx = 0;
        for i in 0..self.sources.len() {
            let source = self.sources[i];
            if !state_map.contains_key(&source){
                state_map.insert(source, source_idx);
                source_idx += 1;
                total_frequencies[state_map[&source]] = self.state_terminating_frequencies[source];
                total_frequencies[state_map[&source]] += self.transition_frequencies[i];
            } 
            else{
                total_frequencies[state_map[&source]] += self.transition_frequencies[i];
            }
        }

        for i in 0..self.targets.len() {
            let target: usize = self.targets[i];
            if !state_map.contains_key(&target){
                state_map.insert(target, source_idx);
                source_idx += 1;
            }
        }

        // Add transitions with computed probabilities to the SDFA
        for i in 0..self.sources.len() {
            let source = self.sources[i];
            let target = self.targets[i];
            let activity = self.activities[i];
            let frequency = self.transition_frequencies[i];
            // Skip transitions from states with zero total frequency (to avoid division by zero)
            if total_frequencies[state_map[&source]] == 0 {
                continue;
            }
            // Convert frequency to probability
            let probability= Fraction::from((frequency, total_frequencies[state_map[&source]]));

            // Add transition to SDFA
            sdfa.add_transition(state_map[&source], activity, state_map[&target], probability)?;
        }

        // Verify that our conversion resulted in valid probabilities
        for state in self.sources.iter() {
            if total_frequencies[state_map[&state]] > 0 {
                let termination_prob = sdfa.get_termination_probability(state_map[&state]);
                if termination_prob.is_negative() {
                    return Err(anyhow!("Invalid probability distribution for state {}: termination probability is negative", state_map[&state]));
                }
            }
        }

        Ok(sdfa)
    }

    pub fn get_activity_key(&mut self) -> &mut ActivityKey {
        &mut self.activity_key
    }

    pub fn get_sources(&self) -> &Vec<usize> {
        &self.sources
    }

    pub fn get_targets(&self) -> &Vec<usize> {
        &self.targets
    }

    pub fn get_activities(&self) -> &Vec<Activity> {
        &self.activities
    }

    pub fn get_transition_frequencies(&self) -> &Vec<usize> {
        &self.transition_frequencies
    }

    pub fn get_state_terminating_frequencies(&self) -> &Vec<usize> {
        &self.state_terminating_frequencies
    }

    pub fn get_size(&self) -> usize {
        // The size of the stochastic transition system is the number of transitions and states
            let mut unique_set = HashSet::new();
    
            // Add all elements from both lists
            for item in self.sources.iter() {
                unique_set.insert(item);
            }
            for item in self.targets.iter() {
                unique_set.insert(item);
            }
            
            unique_set.len() + self.sources.len()
    }

    pub fn get_max_state(&self) -> usize {
        self.max_state
    }

    pub fn set_initial_state(&mut self, state: usize) {
        self.initial_state = state;
    }

    pub fn add_transition(&mut self, source: usize, activity: Activity, target: usize, frequency: usize) -> Result<()> {
        // self.ensure_states(max(source, target), frequency);
        let (found, from) = self.binary_search(source, self.activity_key.get_id_from_activity(activity));
        if found {
            //edge already present
            Err(anyhow!("tried to insert an edge that would violate the determinism of the 
        stochastic transition system"))
        } else {
            self.sources.insert(from, source);
            self.targets.insert(from, target);
            self.activities.insert(from, activity);
            self.state_terminating_frequencies[source] -= frequency;

            self.transition_frequencies.insert(from, frequency);
            Ok(())
        }
    }

    /**
     * Adds the frequency to the transition. Returns the target state, which may be new.
     */
    pub fn take_or_add_transition(&mut self, source_state: usize, activity: Activity, frequency: usize) -> usize {
        // self.state_terminating_frequencies[source_state] += &frequency;

        let (found, transition) = self.binary_search(source_state, self.activity_key.get_id_from_activity(activity));
        if found {
            self.transition_frequencies[transition] += &frequency;
            self.state_terminating_frequencies[self.targets[transition]] += &frequency;
            return self.targets[transition];
        } else {
            let target = self.add_state();
            self.sources.insert(transition, source_state);
            self.targets.insert(transition, target);
            self.activities.insert(transition, activity);
            self.transition_frequencies.insert(transition, frequency);
            self.state_terminating_frequencies.push(frequency);
            return target;
        }
    }

    pub fn get_initial_state(&self) -> usize {
        self.initial_state
    }

    pub fn get_and_add_initial_state(&mut self) -> usize {
        self.state_terminating_frequencies[self.initial_state] += 1;
        self.initial_state
    }

    pub fn update_terminating_frequency(&mut self, state: usize, frequency: usize) {
        self.state_terminating_frequencies[state] -= frequency;
    }


    pub fn get_termination_frequency(&self, state: usize) -> usize {
        self.state_terminating_frequencies[state]
    }

    pub fn add_state(&mut self) -> usize {
		self.max_state += 1;
        return self.max_state;
	}

    fn compare(source1: usize, activity1: usize, source2: usize, activity2: Activity) -> Ordering {
		if source1 < source2 {
			return Ordering::Greater;
		} else if source1 > source2 {
			return Ordering::Less;
		} else if activity2 > activity1 {
			return Ordering::Greater;
		} else if activity2 < activity1 {
			return Ordering::Less;
		} else {
			return Ordering::Equal;
		}
	}

    fn binary_search(&self, source: usize, activity: usize) -> (bool, usize) {
        if self.sources.is_empty() {
            return (false, 0);
        }


        let mut size = self.sources.len();
        let mut left = 0;
        let mut right = size;
        while left < right {
            let mid = left + size / 2;

            let cmp = Self::compare(source, activity, self.sources[mid], self.activities[mid]);

            left = if cmp == Ordering::Less { mid + 1 } else { left };
            right = if cmp == Ordering::Greater { mid } else { right };
            if cmp == Ordering::Equal {
                assert!(mid < self.sources.len());
                return (true, mid);
            }

            size = right - left;
        }

        assert!(left <= self.sources.len());
        (false, left)
	}

    pub fn set_activity_key(&mut self, activity_key: &ActivityKey) {
        self.activity_key = activity_key.clone();
    }

    pub fn get_trace_number(&self) -> usize {
        let mut trace_number = 0;
        // The trace number is the number of transitions in the system
        for i in 0..self.sources.len(){
            if self.sources[i] == self.initial_state {
                trace_number += self.transition_frequencies[i];
            }
        }
        trace_number
    }


   /// Filters out edges with frequencies below the given threshold percentage
    /// Prioritizes least frequent edges, with ties broken by choosing deeper edges (longer prefix)
      /// Filters out edges with frequencies below the given threshold percentage
    /// Prioritizes least frequent edges, with ties broken by choosing deeper edges (longer prefix)
    pub fn filter_by_percentage(&mut self, filter_percentage: f64) -> Result<()> {
        if filter_percentage <= 0.0 || filter_percentage >= 100.0 {
            return Err(anyhow!("Filter percentage must be between 0 and 100"));
        }
        
        // Calculate total frequency and target amount to filter
        let total_frequency: usize = self.get_trace_number();
        let target_filter_amount = (total_frequency as f64 * filter_percentage / 100.0).ceil() as usize;
        
        // Calculate depth (distance from initial state) for each state
        let state_depths = self.calculate_state_depths();
        
        // Create list of edges with their priorities
        let mut edge_priorities = self.create_edge_priority_list(&state_depths);
        
        // Sort by frequency (ascending), then by depth (descending for ties)
        edge_priorities.sort_by(|a, b| {
            let freq_cmp = a.frequency.cmp(&b.frequency);
            if freq_cmp == std::cmp::Ordering::Equal {
                // If frequencies are equal, prefer deeper edges (higher depth)
                b.depth.cmp(&a.depth)
            } else {
                freq_cmp
            }
        });
        
        // Filter edges iteratively until target percentage is reached
        let mut filtered_amount = 0;
        let mut edges_to_remove = Vec::new();
        
        for edge_priority in edge_priorities {
            if filtered_amount >= target_filter_amount {
                break;
            }
            
            edges_to_remove.push(edge_priority.transition_index);
            filtered_amount += edge_priority.frequency;
        }
        
        // Remove the selected edges and update frequencies
        self.remove_edges_and_update_frequencies(edges_to_remove)?;
        
        // Clean up unreachable states
        self.remove_unreachable_states()?;
        
        Ok(())
    }
    
    /// Calculates the depth (shortest distance from initial state) for each state
    fn calculate_state_depths(&self) -> HashMap<usize, usize> {
        let mut depths = HashMap::new();
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        
        // Start BFS from initial state
        queue.push_back((self.initial_state, 0));
        depths.insert(self.initial_state, 0);
        visited.insert(self.initial_state);
        
        while let Some((current_state, current_depth)) = queue.pop_front() {
            // Find all outgoing transitions from current state
            for (transition_idx, &source) in self.sources.iter().enumerate() {
                if source == current_state {
                    let target = self.targets[transition_idx];
                    if !visited.contains(&target) {
                        visited.insert(target);
                        let target_depth = current_depth + 1;
                        depths.insert(target, target_depth);
                        queue.push_back((target, target_depth));
                    }
                }
            }
        }
        
        depths
    }
    
    /// Creates a list of edge priorities with frequency and depth information
    fn create_edge_priority_list(&self, state_depths: &HashMap<usize, usize>) -> Vec<EdgePriority> {
        let mut edge_priorities = Vec::new();
        
        for (transition_idx, &frequency) in self.transition_frequencies.iter().enumerate() {
            let source_state = self.sources[transition_idx];
            let target_state = self.targets[transition_idx];
            
            // The depth of an edge is the depth of its target state
            // (representing how deep this edge goes into the system)
            let depth = state_depths.get(&target_state).unwrap_or(&0);
            
            edge_priorities.push(EdgePriority {
                transition_index: transition_idx,
                frequency,
                depth: *depth,
                source_state,
                target_state,
            });
        }
        
        edge_priorities
    }
    
    /// Removes selected edges and updates frequencies along their prefixes
    fn remove_edges_and_update_frequencies(&mut self, edges_to_remove: Vec<usize>) -> Result<()> {
        // First, update frequencies for the edges we're removing
        for &edge_idx in &edges_to_remove {
            if edge_idx >= self.sources.len() {
                continue; // Skip if already removed
            }
            
            let frequency_to_remove = self.transition_frequencies[edge_idx];
            let target_state = self.targets[edge_idx];
            
            // Update terminating frequency of target state
            if target_state < self.state_terminating_frequencies.len() {
                self.state_terminating_frequencies[target_state] = 
                    self.state_terminating_frequencies[target_state].saturating_sub(frequency_to_remove);
            }
            
            // Find and update all edges in the prefix path leading to this edge
            self.update_prefix_frequencies(self.sources[edge_idx], frequency_to_remove);
        }
        
        // Sort indices in descending order to avoid index shifting issues when removing
        let mut sorted_indices = edges_to_remove;
        sorted_indices.sort_by(|a, b| b.cmp(a));
        
        // Remove the selected edges
        for &edge_idx in &sorted_indices {
            if edge_idx < self.sources.len() {
                self.sources.remove(edge_idx);
                self.targets.remove(edge_idx);
                self.activities.remove(edge_idx);
                self.transition_frequencies.remove(edge_idx);
            }
        }
        
        // Now remove any additional edges that have become zero frequency due to prefix updates
        self.remove_zero_frequency_transitions();
        
        Ok(())
    }
    
    /// Removes all transitions with zero frequency
    fn remove_zero_frequency_transitions(&mut self) {
        let mut indices_to_remove = Vec::new();
        
        for (i, &freq) in self.transition_frequencies.iter().enumerate() {
            if freq == 0 {
                indices_to_remove.push(i);
            }
        }
        
        // Remove in reverse order to avoid index shifting
        indices_to_remove.reverse();
        for idx in indices_to_remove {
            self.sources.remove(idx);
            self.targets.remove(idx);
            self.activities.remove(idx);
            self.transition_frequencies.remove(idx);
        }
    }
    
    /// Updates frequencies along the prefix path leading to a removed edge
    fn update_prefix_frequencies(&mut self, target_state: usize, frequency_to_remove: usize) {
        let mut current_state = target_state;
        let mut visited = HashSet::new();
        
        // Trace back from the target state to the initial state
        while current_state != self.initial_state && !visited.contains(&current_state) {
            visited.insert(current_state);
            
            // Find incoming edges to current_state
            let mut found_incoming = false;
            for (i, &target) in self.targets.iter().enumerate() {
                if target == current_state && self.transition_frequencies[i] >= frequency_to_remove {
                    // Reduce frequency of this incoming edge
                    self.transition_frequencies[i] -= frequency_to_remove;
                    current_state = self.sources[i];
                    found_incoming = true;
                    break;
                }
            }
            
            if !found_incoming {
                break;
            }
        }
    }
    
    /// Removes states that are no longer reachable from the initial state
    fn remove_unreachable_states(&mut self) -> Result<()> {
        let reachable_states = self.find_reachable_states();
        let mut state_mapping = HashMap::new();
        let mut new_state_counter = 0;
        
        // Create mapping from old states to new states (only for reachable states)
        for &state in &reachable_states {
            state_mapping.insert(state, new_state_counter);
            new_state_counter += 1;
        }
        
        // Filter out transitions that reference unreachable states
        let mut transitions_to_keep = Vec::new();
        for i in 0..self.sources.len() {
            let source = self.sources[i];
            let target = self.targets[i];
            
            if state_mapping.contains_key(&source) && state_mapping.contains_key(&target) {
                transitions_to_keep.push(i);
            }
        }
        
        // Keep only transitions between reachable states
        let mut new_sources = Vec::new();
        let mut new_targets = Vec::new();
        let mut new_activities = Vec::new();
        let mut new_frequencies = Vec::new();
        
        for &i in &transitions_to_keep {
            new_sources.push(state_mapping[&self.sources[i]]);
            new_targets.push(state_mapping[&self.targets[i]]);
            new_activities.push(self.activities[i]);
            new_frequencies.push(self.transition_frequencies[i]);
        }
        
        self.sources = new_sources;
        self.targets = new_targets;
        self.activities = new_activities;
        self.transition_frequencies = new_frequencies;
        
        // Update initial state
        if let Some(&new_initial) = state_mapping.get(&self.initial_state) {
            self.initial_state = new_initial;
        } else {
            return Err(anyhow!("Initial state became unreachable after filtering"));
        }
        
        // Compact terminating frequencies array
        let mut new_terminating_frequencies = vec![0; new_state_counter];
        for (&old_state, &new_state) in &state_mapping {
            if old_state < self.state_terminating_frequencies.len() {
                new_terminating_frequencies[new_state] = self.state_terminating_frequencies[old_state];
            }
        }
        self.state_terminating_frequencies = new_terminating_frequencies;
        
        // Update max_state
        self.max_state = if new_state_counter > 0 { new_state_counter - 1 } else { 0 };
        
        Ok(())
    }
    
    /// Finds all states reachable from the initial state
    fn find_reachable_states(&self) -> Vec<usize> {
        let mut reachable = HashSet::new();
        let mut stack = vec![self.initial_state];
        
        while let Some(state) = stack.pop() {
            if reachable.insert(state) {
                // Find all states reachable from this state
                for (i, &source) in self.sources.iter().enumerate() {
                    if source == state {
                        let target = self.targets[i];
                        if !reachable.contains(&target) {
                            stack.push(target);
                        }
                    }
                }
            }
        }
        
        let mut reachable_vec: Vec<usize> = reachable.into_iter().collect();
        reachable_vec.sort();
        reachable_vec
    }
    
    /// Gets statistics about the current system state
    pub fn get_filtering_stats(&self) -> FilteringStats {
        let total_frequency: usize = self.transition_frequencies.iter().sum();
        let num_states = self.find_reachable_states().len();
        let num_transitions = self.sources.len();
        
        FilteringStats {
            total_frequency,
            num_states,
            num_transitions,
            avg_transition_frequency: if num_transitions > 0 {
                total_frequency as f64 / num_transitions as f64
            } else {
                0.0
            },
        }
    }
    
    /// Debug method to show edge priorities for testing
    pub fn debug_edge_priorities(&self) -> Vec<EdgePriority> {
        let state_depths = self.calculate_state_depths();
        let mut edge_priorities = self.create_edge_priority_list(&state_depths);
        
        edge_priorities.sort_by(|a, b| {
            let freq_cmp = a.frequency.cmp(&b.frequency);
            if freq_cmp == std::cmp::Ordering::Equal {
                b.depth.cmp(&a.depth)
            } else {
                freq_cmp
            }
        });
        
        edge_priorities
    }
    
    /// Minimizes the stochastic transition system using Hopcroft's algorithm
    /// Returns the number of states that were merged
     /// Minimizes the stochastic transition system using Hopcroft's algorithm
    /// Returns the number of states that were merged
    pub fn minimize_hopcroft(&mut self) -> Result<usize> {
        // println!("Starting Hopcroft minimization...");
        
        let original_state_count = self.max_state + 1;
        
        // Step 1: Create initial partition based on terminating frequencies
        let mut partition = self.create_initial_partition();
        
        // Step 2: Initialize worklist with all partition blocks and all activities
        let mut worklist = self.initialize_worklist(&partition);
        
        // Step 3: Refine partition until no more refinements possible
        while let Some((block, activity)) = worklist.pop_front() {
            let new_blocks = self.refine_partition(&mut partition, &block, activity);
            
            // Add new refinements to worklist
            for new_block in new_blocks {
                for &act in self.get_all_activities().iter() {
                    worklist.push_back((new_block.clone(), act));
                }
            }
        }
        
        // println!("Final partition has {} blocks", partition.len());
        // for (i, block) in partition.iter().enumerate() {
        //     println!("  Block {}: {:?}", i, block);
        // }
        
        // Step 4: Merge states within each partition block
        let merged_count = self.merge_partition_blocks(partition)?;
        
        // println!("Hopcroft minimization completed. Merged {} states", merged_count);
        Ok(original_state_count - (self.max_state + 1))
    }
    
    /// Creates initial partition based on state behavior
    /// States are grouped by their "signature" - whether they are accept states,
    /// and for non-accept states, by their outgoing activity set
    fn create_initial_partition(&self) -> Vec<HashSet<usize>> {
        let mut signature_to_states: HashMap<String, HashSet<usize>> = HashMap::new();
        
        // Get all states that exist in the system
        let mut all_states = HashSet::new();
        for &source in &self.sources {
            all_states.insert(source);
        }
        for &target in &self.targets {
            all_states.insert(target);
        }
        
        for &state in &all_states {
            let signature = self.compute_state_signature(state);
            signature_to_states
                .entry(signature.clone())
                .or_insert_with(HashSet::new)
                .insert(state);
        }
        
        let partition: Vec<HashSet<usize>> = signature_to_states.values().cloned().collect();
        
        // println!("Initial partition based on state signatures:");
        for (i, block) in partition.iter().enumerate() {
            if let Some(&representative) = block.iter().next() {
                let signature = self.compute_state_signature(representative);
                // println!("  Block {} (signature='{}'): {:?}", i, signature, block);
            }
        }
        
        partition
    }
    
    /// Computes a signature for a state based on its outgoing transitions
    /// Accept states (appear in targets but not in sources) get the signature "ACCEPT"
    /// Other states get a signature based on their outgoing activities
    fn compute_state_signature(&self, state: usize) -> String {
        // Check if this is an accept state (in targets but not in sources)
        let appears_in_targets = self.targets.contains(&state);
        let appears_in_sources = self.sources.contains(&state);
        
        if appears_in_targets && !appears_in_sources {
            // This is an accept state - all accept states should be equivalent initially
            // println!("State {} is an accept state (in targets but not in sources)", state);
            return "ACCEPT".to_string();
        }
        
        // For non-accept states, collect all outgoing activities
        let mut outgoing_activities = HashSet::new();
        for i in 0..self.sources.len() {
            if self.sources[i] == state {
                outgoing_activities.insert(self.activities[i]);
            }
        }
        
        if outgoing_activities.is_empty() {
            // State exists but has no outgoing transitions (might be initial state with no outgoing)
            "NO_OUTGOING".to_string()
        } else {
            // Create signature from sorted activities
            let mut activities: Vec<Activity> = outgoing_activities.into_iter().collect();
            activities.sort();
            
            let activity_labels: Vec<String> = activities.iter()
                .map(|a| self.activity_key.get_activity_label(a).to_string())
                .collect();
            
            format!("ACTIVITIES[{}]", activity_labels.join(","))
        }
    }
    
    /// Initializes the worklist with all (block, activity) pairs
    fn initialize_worklist(&self, partition: &[HashSet<usize>]) -> VecDeque<(HashSet<usize>, Activity)> {
        let mut worklist = VecDeque::new();
        let activities = self.get_all_activities();
        
        for block in partition {
            for &activity in &activities {
                worklist.push_back((block.clone(), activity));
            }
        }
        
        // println!("Initialized worklist with {} items", worklist.len());
        worklist
    }
    
    /// Gets all distinct activities in the transition system
    fn get_all_activities(&self) -> HashSet<Activity> {
        self.activities.iter().cloned().collect()
    }
    
    /// Refines the partition based on a splitter block and activity
    /// Returns new blocks created during refinement
    fn refine_partition(
        &self,
        partition: &mut Vec<HashSet<usize>>,
        splitter: &HashSet<usize>,
        activity: Activity
    ) -> Vec<HashSet<usize>> {
        let mut new_blocks = Vec::new();
        let mut blocks_to_remove = Vec::new();
        
        // println!("  Refining with splitter {:?} and activity {:?}", 
        //         splitter, self.activity_key.get_activity_label(&activity));
        
        // For each block in the current partition
        for (block_idx, block) in partition.iter().enumerate() {
            // Split this block based on transitions to the splitter
            let (block_to_splitter, block_not_to_splitter) = 
                self.split_block_by_transitions(block, splitter, activity);
            
            // If both parts are non-empty, we have a refinement
            if !block_to_splitter.is_empty() && !block_not_to_splitter.is_empty() {
                // println!("    Splitting block {:?}", block);
                // println!("      States with transitions to splitter: {:?}", block_to_splitter);
                // println!("      States without transitions to splitter: {:?}", block_not_to_splitter);
                
                blocks_to_remove.push(block_idx);
                new_blocks.push(block_to_splitter);
                new_blocks.push(block_not_to_splitter);
            } else {
                // println!("    Block {:?} not split (to_splitter={}, not_to_splitter={})", 
                //         block, block_to_splitter.len(), block_not_to_splitter.len());
            }
        }
        
        // Remove split blocks (in reverse order to maintain indices)
        for &idx in blocks_to_remove.iter().rev() {
            partition.remove(idx);
        }
        
        // Add new blocks to partition
        for new_block in &new_blocks {
            partition.push(new_block.clone());
        }
        
        // if !new_blocks.is_empty() {
        //     println!("    Created {} new blocks", new_blocks.len());
        // }
        
        new_blocks
    }
    
    /// Splits a block into states that have transitions to the splitter vs those that don't
    fn split_block_by_transitions(
        &self,
        block: &HashSet<usize>,
        splitter: &HashSet<usize>,
        activity: Activity
    ) -> (HashSet<usize>, HashSet<usize>) {
        let mut to_splitter = HashSet::new();
        let mut not_to_splitter = HashSet::new();
        
        for &state in block {
            let has_transition_to_splitter = self.has_transition_to_block(state, activity, splitter);
            
            if has_transition_to_splitter {
                to_splitter.insert(state);
            } else {
                not_to_splitter.insert(state);
            }
        }
        
        (to_splitter, not_to_splitter)
    }
    
    /// Checks if a state has a transition with given activity to any state in the target block
    fn has_transition_to_block(
        &self,
        state: usize,
        activity: Activity,
        target_block: &HashSet<usize>
    ) -> bool {
        for i in 0..self.sources.len() {
            if self.sources[i] == state && 
               self.activities[i] == activity && 
               target_block.contains(&self.targets[i]) {
                return true;
            }
        }
        false
    }
    
    /// Merges states within each partition block while preserving frequencies
    fn merge_partition_blocks(&mut self, partition: Vec<HashSet<usize>>) -> Result<usize> {
        let mut total_merged = 0;
        
        for (block_id, block) in partition.into_iter().enumerate() {
            if block.len() <= 1 {
                continue;
            }
            
            let mut states: Vec<usize> = block.into_iter().collect();
            states.sort();
            
            // Choose the smallest state number as representative
            let representative = states[0];
            
            // println!("Merging partition block {}: {:?} -> {}", block_id, states, representative);
            
            // Merge all other states into the representative
            for &state_to_merge in &states[1..] {
                self.merge_state_preserving_frequencies(state_to_merge, representative)?;
                total_merged += 1;
            }
        }
        
        // Clean up the transition system
        self.cleanup_after_minimization()?;
        
        Ok(total_merged)
    }
    
    /// Merges one state into another while carefully preserving frequencies
    fn merge_state_preserving_frequencies(
        &mut self,
        source_state: usize,
        target_state: usize
    ) -> Result<()> {
        if source_state == target_state {
            return Ok(());
        }
                
        // Redirect all incoming transitions from other states
        for i in 0..self.targets.len() {
            if self.targets[i] == source_state {
                self.targets[i] = target_state;
            }
        }
        
        // Collect outgoing transitions from source_state
        let mut source_transitions = Vec::new();
        for i in 0..self.sources.len() {
            if self.sources[i] == source_state {
                source_transitions.push((
                    self.activities[i],
                    self.targets[i],
                    self.transition_frequencies[i]
                ));
            }
        }
        
        // Remove all transitions from source_state
        self.remove_transitions_from_state(source_state);
        
        // Process each outgoing transition from source_state
        for (activity, dest, frequency) in source_transitions {
            // Check if target_state already has a transition with same activity and destination
            let mut found_match = false;
            
            for j in 0..self.sources.len() {
                if self.sources[j] == target_state && 
                   self.activities[j] == activity && 
                   self.targets[j] == dest {
                    // Combine frequencies
                    self.transition_frequencies[j] += frequency;
                    found_match = true;
                    // println!("    Combined transition {}->{}({:?}): freq {} + {} = {}", 
                    //         target_state, dest, 
                    //         self.activity_key.get_activity_label(&activity),
                    //         self.transition_frequencies[j] - frequency, frequency, 
                    //         self.transition_frequencies[j]);
                    break;
                }
            }
            
            if !found_match {
                // Add new transition from target_state
                let insert_pos = self.sources.len();
                self.sources.insert(insert_pos, target_state);
                self.targets.insert(insert_pos, dest);
                self.activities.insert(insert_pos, activity);
                self.transition_frequencies.insert(insert_pos, frequency);

                // println!("    Added transition {}->{}({:?}): freq {}", 
                //         target_state, dest, 
                //         self.activity_key.get_activity_label(&activity), frequency);
            }
        }
        
        // Combine terminating frequencies
        let source_term_freq = self.state_terminating_frequencies[source_state];
        self.state_terminating_frequencies[target_state] += source_term_freq;
        self.state_terminating_frequencies[source_state] = 0;
        
        if source_term_freq > 0 {
            // println!("    Combined terminating frequencies: {} + {} = {}", 
            //         self.state_terminating_frequencies[target_state] - source_term_freq,
            //         source_term_freq, 
            //         self.state_terminating_frequencies[target_state]);
        }
        
        // Update initial state if necessary
        if self.initial_state == source_state {
            self.initial_state = target_state;
            println!("    Updated initial state: {} -> {}", source_state, target_state);
        }
        
        Ok(())
    }
    
    /// Clean up after minimization: remove duplicates and update max_state
    fn cleanup_after_minimization(&mut self) -> Result<()> {
        // println!("Cleaning up after minimization...");
        
        // Remove duplicate transitions and combine their frequencies
        let mut i = 0;
        while i < self.sources.len() {
            let source = self.sources[i];
            let target = self.targets[i];
            let activity = self.activities[i];
            
            // Look for duplicates
            let mut j = i + 1;
            while j < self.sources.len() {
                if self.sources[j] == source && 
                   self.targets[j] == target && 
                   self.activities[j] == activity {
                    // Found duplicate - combine frequencies
                    self.transition_frequencies[i] += self.transition_frequencies[j];
                    
                    // Remove the duplicate
                    self.sources.remove(j);
                    self.targets.remove(j);
                    self.activities.remove(j);
                    self.transition_frequencies.remove(j);
                    
                    // println!("Removed duplicate transition {}->{}({:?})", 
                    //         source, target, self.activity_key.get_activity_label(&activity));
                } else {
                    j += 1;
                }
            }
            
            i += 1;
        }
        
        // Update max_state to reflect the current maximum state number
        self.max_state = if self.sources.is_empty() && self.targets.is_empty() {
            0
        } else {
            let max_source = self.sources.iter().max().unwrap_or(&0);
            let max_target = self.targets.iter().max().unwrap_or(&0);
            std::cmp::max(*max_source, *max_target)
        };
        
        // println!("Cleanup complete. New max_state: {}", self.max_state);
        Ok(())
    }
    

}

#[derive(Debug, Clone)]
pub struct EdgePriority {
    pub transition_index: usize,
    pub frequency: usize,
    pub depth: usize,
    pub source_state: usize,
    pub target_state: usize,
}

impl std::fmt::Display for EdgePriority {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Edge {}→{}: freq={}, depth={}, idx={}",
            self.source_state, self.target_state, self.frequency, self.depth, self.transition_index
        )
    }
}

#[derive(Debug, Clone)]
pub struct FilteringStats {
    pub total_frequency: usize,
    pub num_states: usize,
    pub num_transitions: usize,
    pub avg_transition_frequency: f64,
}

impl std::fmt::Display for FilteringStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "FilteringStats {{ total_freq: {}, states: {}, transitions: {}, avg_trans_freq: {:.2} }}",
            self.total_frequency, self.num_states, self.num_transitions, self.avg_transition_frequency
        )
    }
}


impl Infoable for StochasticTransitionSystem {
    fn info(&self, f: &mut impl std::io::Write) -> Result<()> {
        writeln!(f, "Number of states\t{}", self.max_state)?;
        writeln!(f, "Number of transitions\t{}", self.sources.len())?;
        writeln!(f, "Number of activities\t{}", self.activity_key.get_number_of_activities())?;

        Ok(write!(f, "")?)
    }
}

impl fmt::Display for StochasticTransitionSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{{")?;
        writeln!(f, "\"initialState\": {},", self.get_initial_state())?;
        writeln!(f, "\"transitions\": [")?;
        for pos in 0..self.sources.len() {
            write!(f, "{{\"from\":{},\"to\":{},\"label\":\"{}\",\"freq\":\"{}\"}}", 
                self.sources[pos], 
                self.targets[pos], 
                self.activity_key.get_activity_label(&self.activities[pos]),
                self.transition_frequencies[pos])?;
            if pos + 1 == self.sources.len() {
                writeln!(f, "")?;
            } else {
                writeln!(f, ",")?;
            }
        }
        writeln!(f, "]}}")?;
        Ok(())
    }
}

impl EbiTraitGraphable for StochasticTransitionSystem {
    fn to_dot(&self) -> Result<layout::topo::layout::VisualGraph> {
        log::info!("to_dot for StochasticTransitionSystem");
        let mut graph = VisualGraph::new(layout::core::base::Orientation::LeftToRight);

        // Create a place for each state
        let mut places = vec![];
        for state in 0..=self.max_state {
            places.push(<dyn EbiTraitGraphable>::create_place(
                &mut graph,
                &format!("{}", self.state_terminating_frequencies[state]),
            ));
        }

        // Add transitions as edges between places
        for pos in 0..self.sources.len() {
            let source = places[self.sources[pos]];
            let target = places[self.targets[pos]];
            let frequency = &self.transition_frequencies[pos];
            let activity = self.activity_key.get_activity_label(&self.activities[pos]);

            <dyn EbiTraitGraphable>::create_edge(
                &mut graph,
                &source,
                &target,
                &format!("{}, {}", activity, frequency),
            );
        }

        Ok(graph)
    }
}
   

impl<'a> IntoIterator for &'a StochasticTransitionSystem {
    type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

    type IntoIter = StochasticTransitionSystemIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            it_sources:  self.sources.iter(),
            it_targets:  self.targets.iter(),
            it_activities: self.activities.iter(),
            it_transition_frequencies: self.transition_frequencies.iter()
        }
    }
}

pub struct StochasticTransitionSystemIterator<'a> {
    it_sources: std::slice::Iter<'a, usize>,
    it_targets: std::slice::Iter<'a, usize>,
    it_activities: std::slice::Iter<'a, Activity>,
    it_transition_frequencies: std::slice::Iter<'a, usize>
}

impl<'a> Iterator for StochasticTransitionSystemIterator<'a> {
    type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(source) = self.it_sources.next() {
            let target = self.it_targets.next().unwrap();
            let activity = self.it_activities.next().unwrap();
            let frequency = self.it_transition_frequencies.next().unwrap();
            let result = Some((source, target, activity, frequency));
            result
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a mut StochasticTransitionSystem {
    type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

    type IntoIter = StochasticTransitionSystemMutIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            it_sources:  self.sources.iter(),
            it_targets:  self.targets.iter(),
            it_activities: self.activities.iter(),
            it_transition_frequencies: self.transition_frequencies.iter()
        }
    }
}

pub struct StochasticTransitionSystemMutIterator<'a> {
    it_sources: std::slice::Iter<'a, usize>,
    it_targets: std::slice::Iter<'a, usize>,
    it_activities: std::slice::Iter<'a, Activity>,
    it_transition_frequencies: std::slice::Iter<'a, usize>
}

impl<'a> Iterator for StochasticTransitionSystemMutIterator<'a> {
    type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(source) = self.it_sources.next() {
            let target = self.it_targets.next().unwrap();
            let activity = self.it_activities.next().unwrap();
            let frequency = self.it_transition_frequencies.next().unwrap();
            let result = Some((source, target, activity, frequency));
            result
        } else {
            None
        }
    }
}

impl Clone for StochasticTransitionSystem {
    fn clone(&self) -> Self {
        Self {
            activity_key: self.activity_key.clone(),
            initial_state: self.initial_state,
            max_state: self.max_state,
            sources: self.sources.clone(),
            targets: self.targets.clone(),
            activities: self.activities.clone(),
            transition_frequencies: self.transition_frequencies.clone(),
            state_terminating_frequencies: self.state_terminating_frequencies.clone(),
        }
    }
}