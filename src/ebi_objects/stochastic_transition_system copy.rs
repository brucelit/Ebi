// use crate::{
//     ebi_framework::{
//         activity_key::{Activity, ActivityKey}, infoable::Infoable
//     }, ebi_objects::stochastic_deterministic_finite_automaton::StochasticDeterministicFiniteAutomaton, ebi_traits::ebi_trait_graphable::EbiTraitGraphable, math::{fraction::Fraction, traits::Signed}
    
// };
// use anyhow::{Result, anyhow};
// use layout::topo::layout::VisualGraph;
// use std::{
//     cmp::{max, Ordering},
//     collections::{HashMap, HashSet},
//     fmt
// };

// pub const FORMAT_SPECIFICATION: &str = "A stochastic transition system is a JSON structure with the top level being an object.
//     This object contains the following key-value pairs:
//     \\begin{itemize}
//     \\item \\texttt{initialState} being the index of the initial state. This field is optional: if omitted, the SDFA has an empty stochastic language.
//     \\item \\texttt{transitions} being a list of transitions. 
//     Each transition is an object with \\texttt{from} being the source state index of the transition, 
//     \\texttt{to} being the target state index of the transition, 
//     \\texttt{label} being the activity of the transition, and
//     \\texttt{freq} being the frequency of the transition (may be given as a fraction in a string or a float value. Must be $\\leq 1$). 
//     Silent transitions are not supported.
//     The file format supports deadlocks and livelocks.
//     The frequency that a trace terminates in a state is the sum of frequency of the incoming transitions - the sum probability of the outgoing transitions of the state.
//     \\end{itemize}
//     For instance:
//     \\lstinputlisting[language=json, style=boxed]{../testfiles/aa-ab-ba.sdfa}";

// #[derive(Debug)]
// pub struct StochasticTransitionSystem {
//     activity_key: ActivityKey,
//     initial_state: usize,
//     max_state: usize,
//     sources: Vec<usize>, //transition -> source of arc
//     targets: Vec<usize>, //transition -> target of arc
//     activities: Vec<Activity>, //transition -> activity of arc (every transition is labelled)
//     transition_frequencies: Vec<usize>, //transition -> frequency of arc
//     state_terminating_frequencies: Vec<usize>, //state -> termination frequency
// }

// impl StochasticTransitionSystem {

//     pub fn new() -> Self {
//         Self {
//             activity_key: ActivityKey::new(),
//             max_state: 0,
//             initial_state: 0,
//             sources: vec![],
//             targets: vec![],
//             activities: vec![],
//             transition_frequencies: vec![],
//             state_terminating_frequencies: vec![0],
//         }
//     }

//     /// Breaks a connection between a red state and a blue state.
//     /// Instead of merging them, we:
//     /// 1. Remove the transition from red to blue
//     /// 2. Add a self-loop on the red state with the same activity and frequency
//     /// 3. Add the transition frequency to the red state's terminating frequency
//     ///
//     /// # Parameters
//     /// * `red_state` - The red state in the connection
//     /// * `blue_state` - The blue state to disconnect
//     /// * `red_states` - Set of current red states (not modified)
//     /// * `blue_states` - Set of current blue states (may be modified if blue state becomes disconnected)
//     ///
//     /// # Returns
//     /// * `Result<()>` - Success or an error message
//     pub fn trial_merge<'a>(
//         &mut self,
//         red_state: usize,
//         blue_state: usize,
//         red_states: &'a mut HashSet<usize>,
//         blue_states: &'a mut HashSet<usize>
//     )  {

//         // Find transitions from red_state to blue_state
//         let mut transitions_to_remove = Vec::new();
//         let mut transition_info = Vec::new(); // (activity, frequency)
        
//         for i in 0..self.sources.len() {
//             if red_states.contains(&self.sources[i])  && self.targets[i] == blue_state {
//                 transitions_to_remove.push(i);
//                 transition_info.push((self.activities[i], self.transition_frequencies[i]));
//             }
//         }
        
//         if !transitions_to_remove.is_empty(){
//                     // Step 3: Remove the original transitions
//             // Sort in descending order to avoid index shifting during removal
//             transitions_to_remove.sort_by(|a, b| b.cmp(a)); 
            
//             for idx in transitions_to_remove {
//                 self.sources.remove(idx);
//                 self.targets.remove(idx);
//                 self.activities.remove(idx);
//                 self.transition_frequencies.remove(idx);
//             }
            
//             // Step 1: Add self-loops on the red state
//             for (activity, frequency) in &transition_info {
//                 // Check if a self-loop with this activity already exists
//                 let mut existing_self_loop = false;
//                 let mut existing_idx = 0;
                
//                 for i in 0..self.sources.len() {
//                     if self.sources[i] == red_state && self.targets[i] == red_state && self.activities[i] == *activity {
//                         existing_self_loop = true;
//                         existing_idx = i;
//                         break;
//                     }
//                 }
//                 if existing_self_loop {
//                     // If a self-loop already exists, increase its frequency
//                     self.transition_frequencies[existing_idx] += frequency;
//                 } else {
//                     // Otherwise, create a new self-loop
//                     // println!("Adding self-loop on red state {} with activity {} and frequency {}", red_state, activity, frequency);
//                     let (found, idx) = self.binary_search(red_state, self.activity_key.get_id_from_activity(*activity));
//                     if !found {
//                         self.sources.insert(idx, red_state);
//                         self.targets.insert(idx, red_state);
//                         self.activities.insert(idx, *activity);
//                         self.transition_frequencies.insert(idx, *frequency);
//                     }
//                 }
                
//             }
//         }
        
//         // Step 4: Update blue_states
//         // Check if blue_state is still connected to any red state
//         let mut blue_state_still_connected = false;
        
//         for i in 0..self.sources.len() {
//             // Check if any red state has a transition to blue_state
//             if red_states.contains(&self.sources[i]) && self.targets[i] == blue_state {
//                 blue_state_still_connected = true;
//                 break;
//             }
//         }
        
//         // If blue_state is no longer connected to any red state, it's no longer blue
//         if !blue_state_still_connected {
//             blue_states.remove(&blue_state);
            
//             // Find any new blue states (direct successors of red states)
//             let new_blue_states = self.get_blue_states(red_states);
//             for state in new_blue_states {
//                 if !blue_states.contains(&state) {
//                     blue_states.insert(state);
//                 }
//             }
//         }

//         // Get all states reachable from the red state
//         let mut reachable_states_from_red = HashSet::new();
//         // let mut edge_num_from_red = 0;
//         // let mut edge_num_from_blue = 0;

//         let mut queue = vec![red_state];
//         while let Some(state) = queue.pop() {
//             if reachable_states_from_red.insert(state) {
//                 for i in 0..self.sources.len() {
//                     if self.sources[i] == state && !reachable_states_from_red.contains(&self.targets[i]) {
//                         queue.push(self.targets[i]);
//                     }
//                     // if self.sources[i] == state {
//                     //     edge_num_from_red+=1;
//                     // }
//                 }
//             }
//         }

//         // Get all states reachable from the blue state
//         let mut reachable_states_from_blue = HashSet::new();
//         let mut queue = vec![blue_state];
//         while let Some(state) = queue.pop() {
//             if reachable_states_from_blue.insert(state) {
//                 for i in 0..self.sources.len() {
//                     if self.sources[i] == state && !reachable_states_from_blue.contains(&self.targets[i]) {
//                         queue.push(self.targets[i]);
//                     }
//                     // if self.sources[i] == state {
//                     //     edge_num_from_blue+=1;
//                     // }
//                 }
//             }
//         }        
//         let mut state2remove = HashSet::new();
//         state2remove = self.merge_states(red_state, 
//         blue_state, 
//         state2remove);
//         for state in state2remove.clone() {
//             for i in 0..self.sources.len() {
//                 if self.sources[i] == state {
//                     self.sources.remove(i);
//                     self.targets.remove(i);
//                     self.activities.remove(i);
//                     self.transition_frequencies.remove(i);
//                     break;
//                 }
//             }
//         }
//         self.max_state -= state2remove.len();
//         // update the blue states
//     }


//     pub fn finish_merge<'a>(
//         &mut self,
//         red_state: usize,
//         blue_state: usize,
//         red_states: &'a mut HashSet<usize>,
//         blue_states: &'a mut HashSet<usize>
//     ){
//         // Find transitions from red_state to blue_state
//         let mut transitions_to_remove = Vec::new();
//         let mut transition_info = Vec::new(); // (activity, frequency)

//         for i in 0..self.sources.len() {
//             if red_states.contains(&self.sources[i])  && self.targets[i] == blue_state {
//                 transitions_to_remove.push(i);
//                 transition_info.push((self.activities[i], self.transition_frequencies[i]));
//             }
//         }
        
//         if !transitions_to_remove.is_empty(){
//                     // Step 3: Remove the original transitions
//             // Sort in descending order to avoid index shifting during removal
//             transitions_to_remove.sort_by(|a, b| b.cmp(a)); 
            
//             for idx in transitions_to_remove {
//                 self.sources.remove(idx);
//                 self.targets.remove(idx);
//                 self.activities.remove(idx);
//                 self.transition_frequencies.remove(idx);
//             }
            
//             // Step 1: Add self-loops on the red state
//             for (activity, frequency) in &transition_info {
//                 // Check if a self-loop with this activity already exists
//                 let mut existing_self_loop = false;
//                 let mut existing_idx = 0;
                
//                 for i in 0..self.sources.len() {
//                     if self.sources[i] == red_state && self.targets[i] == red_state && self.activities[i] == *activity {
//                         existing_self_loop = true;
//                         existing_idx = i;
//                         break;
//                     }
//                 }
//                 if existing_self_loop {
//                     // If a self-loop already exists, increase its frequency
//                     self.transition_frequencies[existing_idx] += frequency;
//                 } else {
//                     // Otherwise, create a new self-loop
//                     // println!("Adding self-loop on red state {} with activity {} and frequency {}", red_state, activity, frequency);
//                     let (found, idx) = self.binary_search(red_state, self.activity_key.get_id_from_activity(*activity));
//                     if !found {
//                         self.sources.insert(idx, red_state);
//                         self.targets.insert(idx, red_state);
//                         self.activities.insert(idx, *activity);
//                         self.transition_frequencies.insert(idx, *frequency);
//                     }
//                 }
                
//             }
//         }
        
//         let mut state2remove = HashSet::new();
//         state2remove = self.merge_states(red_state, 
//         blue_state, 
//         state2remove);
//         for state in state2remove.clone() {
//             for i in 0..self.sources.len() {
//                 if self.sources[i] == state {
//                     self.sources.remove(i);
//                     self.targets.remove(i);
//                     self.activities.remove(i);
//                     self.transition_frequencies.remove(i);
//                     break;
//                 }
//             }
//         }
//         self.state_terminating_frequencies[red_state] += self.state_terminating_frequencies[blue_state];
//         self.max_state -= state2remove.len();
//     }
    

//     pub fn merge_states(&mut self,
//         source_state1: usize,
//         source_state2: usize,
//         mut state2remove:  HashSet<usize>) ->  HashSet<usize> {
//             // Merge the frequency of blue into the red
//             self.state_terminating_frequencies[source_state1] += self.state_terminating_frequencies[source_state2];

//             // iterate over outgoing transitions of the blue state
//             for i in 0..self.sources.len() {
//                 if self.sources[i] == source_state2 {
//                     let mut find_common = false;

//                     // iterate over the outgoing transitions of the red state
//                     for j in 0..self.sources.len() {
//                         if self.sources[j] == source_state1 && self.activities[j] == self.activities[i] {
//                             state2remove.insert(source_state2);
//                             find_common= true;
//                             // merge the frequencies
//                             self.transition_frequencies[j] += self.transition_frequencies[i];

//                             // continue the merging with both target states
//                             let target1 = self.targets[j];
//                             let target2 = self.targets[i];
//                             state2remove = self.merge_states(
//                                 target1, 
//                                 target2, state2remove);                           
//                         }
//                     }
//                     if find_common == false{
//                         self.sources[i] = source_state1;
//                         state2remove.insert(source_state2);
//                     }
//                 }
//             }
//             state2remove
//         }


//     /// Find all blue states that are direct successors of red states
//     pub fn get_blue_states(&self, red_states: &HashSet<usize>) -> HashSet<usize> {
//         let mut blue_states = HashSet::new();
        
//         // Find all states that are direct successors of red states
//         for i in 0..self.sources.len() {
//             if red_states.contains(&self.sources[i]) {
//                 let target = self.targets[i];
//                 if !red_states.contains(&target) {
//                     blue_states.insert(target);
//                 }
//             }
//         }
        
//         blue_states
//     }
    
//     /// Initialize state coloring for grammatical inference
//     /// Returns (red_states, blue_states)
//     pub fn initialize_state_colors(&self) -> HashSet<usize> {
//         let mut red_states = HashSet::new();
        
//         // Initially, only the initial state is red
//         red_states.insert(self.initial_state);
        
//         // Blue states are direct successors of red states
//         let blue_states = self.get_blue_states(&red_states);
    
//         blue_states
//     }
    
//     /// Print the current state color sets (for debugging)
//     pub fn print_state_colors(&self, red_states: &HashSet<usize>, blue_states: &HashSet<usize>) {
        
//         // Calculate white states (all states not in red or blue)
//         let mut all_states = HashSet::new();
//         for i in 0..self.sources.len() {
//             all_states.insert(self.sources[i]);
//             all_states.insert(self.targets[i]);
//         }
        
//         let mut white_states = all_states.clone();
//         for &state in red_states {
//             white_states.remove(&state);
//         }
//         for &state in blue_states {
//             white_states.remove(&state);
//         }
//         println!("White states: {:?}", white_states);
//     }
    

//     pub fn filter_by_frequency_and_depth(&mut self, transitions_to_remove: usize) -> Result<()> {
//         // If we have fewer transitions than requested to remove, return an error
//         let total_transitions = self.sources.len();
//         if total_transitions <= transitions_to_remove {
//             return Err(anyhow!("Cannot remove {} transitions, only {} transitions exist", 
//                               transitions_to_remove, total_transitions));
//         }
        
//         println!("Starting to filter {} transitions out of {}", transitions_to_remove, total_transitions);
        
//         // Track how many transitions we've removed
//         let mut removed_count = 0;
        
//         // Maintain a set of removed states to avoid processing them
//         let mut removed_states = HashSet::new();
        
//         // Continue filtering until we've removed the requested number of transitions
//         while removed_count < transitions_to_remove {
//             // Find the candidate transition to remove (least frequent, max depth)
//             if let Some((transition_idx, frequency)) = self.find_candidate_transition(&removed_states) {
//                 // First find the path from initial state to the source of this transition
//                 let path = self.find_path_to_state(self.sources[transition_idx], &removed_states);
                
//                 // Remove the transition and update frequencies
//                 self.remove_transition_and_update(transition_idx, frequency, &path, &mut removed_states);
                
//                 // Update our counter
//                 removed_count += 1;
                
//                 // Log progress occasionally
//                 if removed_count % 100 == 0 || removed_count == transitions_to_remove {
//                     println!("Removed {} transitions out of {}", removed_count, transitions_to_remove);
//                 }
//             } else {
//                 return Err(anyhow!("Could only remove {} transitions of the requested {}", 
//                                   removed_count, transitions_to_remove));
//             }
//         }
        
//         // Now clean up the data structure by removing the transitions and states
//         self.clean_after_filtering(&removed_states);
        
//         println!("Successfully removed {} transitions", removed_count);
//         Ok(())
//     }
    
//     /// Find the path from the initial state to the given state
//     /// Returns a vector of transition indices that lead to this state
//     fn find_path_to_state(&self, target_state: usize, removed_states: &HashSet<usize>) -> Vec<usize> {
//         if target_state == self.initial_state {
//             return Vec::new(); // Empty path for initial state
//         }
        
//         // Use BFS to find the path
//         let mut queue = Vec::new();
//         let mut visited = HashSet::new();
//         let mut parent_transition = HashMap::new(); // state -> (parent_state, transition_idx)
        
//         queue.push(self.initial_state);
//         visited.insert(self.initial_state);
        
//         while let Some(current) = queue.pop() {
//             // If we found the target state, construct the path
//             if current == target_state {
//                 let mut path = Vec::new();
//                 let mut state = current;
                
//                 while state != self.initial_state {
//                     let (parent, trans_idx) = parent_transition[&state];
//                     path.push(trans_idx);
//                     state = parent;
//                 }
                
//                 // Reverse to get the path from initial state to target
//                 path.reverse();
//                 return path;
//             }
            
//             // Check all outgoing transitions
//             for i in 0..self.sources.len() {
//                 if self.sources[i] == current {
//                     let next_state = self.targets[i];
                    
//                     if !visited.contains(&next_state) && !removed_states.contains(&next_state) {
//                         visited.insert(next_state);
//                         queue.push(next_state);
//                         parent_transition.insert(next_state, (current, i));
//                     }
//                 }
//             }
//         }
        
//         // If we get here, there's no path to the target state
//         Vec::new()
//     }
    
//     /// Find a candidate transition to remove: the least frequent transition at the maximum depth
//     fn find_candidate_transition(&self, removed_states: &HashSet<usize>) -> Option<(usize, usize)> {
//         // Step 1: Find the maximum depth of any remaining transition
//         let depths = self.calculate_state_depths(removed_states);
//         let mut max_depth = 0;
        
//         for i in 0..self.sources.len() {
//             let source = self.sources[i];
            
//             // Skip transitions where source or target is already removed
//             if removed_states.contains(&source) || removed_states.contains(&self.targets[i]) {
//                 continue;
//             }
            
//             if let Some(&source_depth) = depths.get(&source) {
//                 if source_depth > max_depth {
//                     max_depth = source_depth;
//                 }
//             }
//         }
        
//         // Step 2: Among transitions at the maximum depth, find the one with minimum frequency
//         let mut min_frequency = usize::MAX;
//         let mut candidate = None;
        
//         for i in 0..self.sources.len() {
//             let source = self.sources[i];
//             let target = self.targets[i];
            
//             // Skip transitions where source or target is already removed
//             if removed_states.contains(&source) || removed_states.contains(&target) {
//                 continue;
//             }
            
//             if let Some(&source_depth) = depths.get(&source) {
//                 if source_depth == max_depth {
//                     let frequency = self.transition_frequencies[i];
//                     if frequency < min_frequency {
//                         min_frequency = frequency;
//                         candidate = Some((i, frequency));
//                     }
//                 }
//             }
//         }
        
//         candidate
//     }
    
//     /// Calculate the depth of each state from the initial state, skipping removed states
//     fn calculate_state_depths(&self, removed_states: &HashSet<usize>) -> HashMap<usize, usize> {
//         let mut depths = HashMap::new();
//         let mut queue = Vec::new();
        
//         // Start with the initial state at depth 0
//         depths.insert(self.initial_state, 0);
//         queue.push(self.initial_state);
        
//         while let Some(state) = queue.pop() {
//             let current_depth = depths[&state];
            
//             // Look for outgoing transitions
//             for i in 0..self.sources.len() {
//                 if self.sources[i] == state {
//                     let target = self.targets[i];
                    
//                     // Skip removed states
//                     if removed_states.contains(&target) {
//                         continue;
//                     }
                    
//                     // If target's depth is not yet set, or we found a shorter path
//                     if !depths.contains_key(&target) || depths[&target] > current_depth + 1 {
//                         depths.insert(target, current_depth + 1);
//                         queue.push(target);
//                     }
//                 }
//             }
//         }
        
//         depths
//     }
    
//     /// Remove a transition and update frequencies along the path
//     fn remove_transition_and_update(
//         &mut self,
//         transition_idx: usize,
//         frequency: usize,
//         path: &[usize],
//         removed_states: &mut HashSet<usize>
//     ) {
//         // Mark the target state as removed
//         let target_state = self.targets[transition_idx];
//         removed_states.insert(target_state);
        
//         // Update frequencies along the path
//         for &path_transition in path {
//             // Deduct the frequency from this transition
//             if self.transition_frequencies[path_transition] > frequency {
//                 self.transition_frequencies[path_transition] -= frequency;
//             } else {
//                 // If the frequency would go to zero or negative, mark the transition's target for removal too
//                 // (this shouldn't happen in a properly formed tree, but let's be safe)
//                 self.transition_frequencies[path_transition] = 0;
//                 removed_states.insert(self.targets[path_transition]);
//             }
//         }
        
//         // Set the frequency of the removed transition to zero (we'll physically remove it later)
//         self.transition_frequencies[transition_idx] = 0;
        
//         // Also update the state terminating frequencies for states along the path
//         // In a tree structure, each state's frequency equals its incoming transition frequency
//         for i in 0..self.sources.len() {
//             // For each source on our path, we need to update its terminating frequency
//             if path.contains(&i) || i == transition_idx {
//                 let source = self.sources[i];
//                 if source < self.state_terminating_frequencies.len() {
//                     // Deduct the frequency from the terminating frequency
//                     if self.state_terminating_frequencies[source] >= frequency {
//                         self.state_terminating_frequencies[source] -= frequency;
//                     } else {
//                         self.state_terminating_frequencies[source] = 0;
//                     }
//                 }
//             }
//         }
//     }
    
//     /// Clean up the data structure after filtering
//     fn clean_after_filtering(&mut self, removed_states: &HashSet<usize>) {
//         // Create new vectors for the filtered data
//         let mut new_sources = Vec::new();
//         let mut new_targets = Vec::new();
//         let mut new_activities = Vec::new();
//         let mut new_transition_frequencies = Vec::new();
        
//         // Keep only transitions where neither source nor target is removed
//         for i in 0..self.sources.len() {
//             let source = self.sources[i];
//             let target = self.targets[i];
            
//             if !removed_states.contains(&source) && !removed_states.contains(&target) {
//                 new_sources.push(source);
//                 new_targets.push(target);
//                 new_activities.push(self.activities[i]);
//                 new_transition_frequencies.push(self.transition_frequencies[i]);
//             }
//         }
        
//         // Update the vectors
//         self.sources = new_sources;
//         self.targets = new_targets;
//         self.activities = new_activities;
//         self.transition_frequencies = new_transition_frequencies;
        
//         // Identify all connected states (states with incoming or outgoing transitions)
//         let mut connected_states = HashSet::new();
        
//         // Always include the initial state
//         connected_states.insert(self.initial_state);
        
//         // Add states with transitions
//         for i in 0..self.sources.len() {
//             connected_states.insert(self.sources[i]);
//             connected_states.insert(self.targets[i]);
//         }
        
//         // Create a combined set of all states to remove
//         let mut all_states_to_remove = removed_states.clone();
        
//         // Add isolated states (states that aren't connected and aren't the initial state)
//         for state in 0..=self.max_state {
//             if !connected_states.contains(&state) && state != self.initial_state {
//                 all_states_to_remove.insert(state);
//             }
//         }
        
//         // Set terminating frequencies for removed and isolated states to zero
//         for &state in &all_states_to_remove {
//             if state < self.state_terminating_frequencies.len() {
//                 self.state_terminating_frequencies[state] = 0;
//             }
//         }
        
//         // Update max_state based on remaining states
//         if !all_states_to_remove.is_empty() {
//             let mut max_remaining_state = 0;
//             for state in 0..=self.max_state {
//                 if !all_states_to_remove.contains(&state) {
//                     max_remaining_state = max(max_remaining_state, state);
//                 }
//             }
//             self.max_state = max_remaining_state;
//         }
        
//         // Log the cleanup results
//         let isolated_count = all_states_to_remove.len() - removed_states.len();
//         println!("Cleaned up after filtering: removed {} states from transitions, identified {} isolated states",
//                  removed_states.len(), isolated_count);
//     }


//     pub fn to_stochastic_deterministic_finite_automaton(&self) -> Result<StochasticDeterministicFiniteAutomaton> {
//         // Create a new SDFA with the same activity key
//         let mut sdfa = StochasticDeterministicFiniteAutomaton::new();
//         sdfa.set_activity_key(&self.activity_key);
//         sdfa.set_initial_state(Some(self.initial_state));
        
//         // Calculate total outgoing frequency for each state
//         // This includes both outgoing transitions and termination frequencies
//         let mut total_frequencies = vec![0usize; self.max_state + 1];
//         let mut state_map = HashMap::new();

//         // Map states from sts to sdfa
//         let mut source_idx = 0;
//         for i in 0..self.sources.len() {
//             let source = self.sources[i];
//             if !state_map.contains_key(&source){
//                 state_map.insert(source, source_idx);
//                 source_idx += 1;
//                 total_frequencies[state_map[&source]] = self.state_terminating_frequencies[source];
//                 total_frequencies[state_map[&source]] += self.transition_frequencies[i];
//             } 
//             else{
//                 total_frequencies[state_map[&source]] += self.transition_frequencies[i];
//             }
//         }

//         for i in 0..self.targets.len() {
//             let target: usize = self.targets[i];
//             if !state_map.contains_key(&target){
//                 state_map.insert(target, source_idx);
//                 source_idx += 1;
//             }
//         }

//         // Add transitions with computed probabilities to the SDFA
//         for i in 0..self.sources.len() {
//             let source = self.sources[i];
//             let target = self.targets[i];
//             let activity = self.activities[i];
//             let frequency = self.transition_frequencies[i];
//             // Skip transitions from states with zero total frequency (to avoid division by zero)
//             if total_frequencies[state_map[&source]] == 0 {
//                 continue;
//             }
//             // Convert frequency to probability
//             let probability= Fraction::from((frequency, total_frequencies[state_map[&source]]));

//             // Add transition to SDFA
//             sdfa.add_transition(state_map[&source], activity, state_map[&target], probability)?;
//         }

//         // Verify that our conversion resulted in valid probabilities
//         for state in self.sources.iter() {
//             if total_frequencies[state_map[&state]] > 0 {
//                 let termination_prob = sdfa.get_termination_probability(state_map[&state]);
//                 if termination_prob.is_negative() {
//                     return Err(anyhow!("Invalid probability distribution for state {}: termination probability is negative", state_map[&state]));
//                 }
//             }
//         }

//         Ok(sdfa)
//     }

//     pub fn get_activity_key(&mut self) -> &mut ActivityKey {
//         &mut self.activity_key
//     }

//     pub fn get_sources(&self) -> &Vec<usize> {
//         &self.sources
//     }

//     pub fn get_targets(&self) -> &Vec<usize> {
//         &self.targets
//     }

//     pub fn get_activities(&self) -> &Vec<Activity> {
//         &self.activities
//     }

//     pub fn get_transition_frequencies(&self) -> &Vec<usize> {
//         &self.transition_frequencies
//     }

//     pub fn get_state_terminating_frequencies(&self) -> &Vec<usize> {
//         &self.state_terminating_frequencies
//     }

//     pub fn get_size(&self) -> usize {
//         // The size of the stochastic transition system is the number of transitions and states
//         self.max_state + self.transition_frequencies.len()
//     }

//     pub fn get_max_state(&self) -> usize {
//         self.max_state
//     }

//     pub fn set_initial_state(&mut self, state: usize) {
//         self.initial_state = state;
//     }

//     // fn ensure_states(&mut self, new_max_state: usize, frequency: usize) {
//     //     if new_max_state > self.max_state {
//     //         self.state_terminating_frequencies.extend(vec![frequency; new_max_state - self.max_state]);
//     //         self.max_state = new_max_state;
//     //         assert!(self.state_terminating_frequencies.len() == self.max_state + 1)
//     //     }
//     // }

//     pub fn add_transition(&mut self, source: usize, activity: Activity, target: usize, frequency: usize) -> Result<()> {
//         // self.ensure_states(max(source, target), frequency);
//         let (found, from) = self.binary_search(source, self.activity_key.get_id_from_activity(activity));
//         if found {
//             //edge already present
//             Err(anyhow!("tried to insert an edge that would violate the determinism of the 
//         stochastic transition system"))
//         } else {
//             self.sources.insert(from, source);
//             self.targets.insert(from, target);
//             self.activities.insert(from, activity);
//             println!("source: {}", source);
//             println!("target: {}", target);
//             println!("frequency: {}", frequency);
//             self.state_terminating_frequencies[source] -= frequency;
//             println!("state terminating frequencies: {}", self.state_terminating_frequencies[source]);
//             self.transition_frequencies.insert(from, frequency);
//             Ok(())
//         }
//     }

//     /**
//      * Adds the frequency to the transition. Returns the target state, which may be new.
//      */
//     pub fn take_or_add_transition(&mut self, source_state: usize, activity: Activity, frequency: usize) -> usize {
//         // self.state_terminating_frequencies[source_state] += &frequency;

//         let (found, transition) = self.binary_search(source_state, self.activity_key.get_id_from_activity(activity));
//         if found {
//             self.transition_frequencies[transition] += &frequency;
//             self.state_terminating_frequencies[self.targets[transition]] += &frequency;
//             return self.targets[transition];
//         } else {
//             let target = self.add_state();
//             self.sources.insert(transition, source_state);
//             self.targets.insert(transition, target);
//             self.activities.insert(transition, activity);
//             self.transition_frequencies.insert(transition, frequency);
//             self.state_terminating_frequencies.push(frequency);
//             return target;
//         }
//     }

//     pub fn get_initial_state(&self) -> usize {
//         self.initial_state
//     }

//     pub fn get_and_add_initial_state(&mut self) -> usize {
//         self.state_terminating_frequencies[self.initial_state] += 1;
//         self.initial_state
//     }

//     pub fn update_terminating_frequency(&mut self, state: usize, frequency: usize) {
//         self.state_terminating_frequencies[state] -= frequency;
//     }


//     pub fn get_termination_frequency(&self, state: usize) -> usize {
//         self.state_terminating_frequencies[state]
//     }

//     pub fn add_state(&mut self) -> usize {
// 		self.max_state += 1;
//         return self.max_state;
// 	}

//     fn compare(source1: usize, activity1: usize, source2: usize, activity2: Activity) -> Ordering {
// 		if source1 < source2 {
// 			return Ordering::Greater;
// 		} else if source1 > source2 {
// 			return Ordering::Less;
// 		} else if activity2 > activity1 {
// 			return Ordering::Greater;
// 		} else if activity2 < activity1 {
// 			return Ordering::Less;
// 		} else {
// 			return Ordering::Equal;
// 		}
// 	}

//     fn binary_search(&self, source: usize, activity: usize) -> (bool, usize) {
//         if self.sources.is_empty() {
//             return (false, 0);
//         }


//         let mut size = self.sources.len();
//         let mut left = 0;
//         let mut right = size;
//         while left < right {
//             let mid = left + size / 2;

//             let cmp = Self::compare(source, activity, self.sources[mid], self.activities[mid]);

//             left = if cmp == Ordering::Less { mid + 1 } else { left };
//             right = if cmp == Ordering::Greater { mid } else { right };
//             if cmp == Ordering::Equal {
//                 assert!(mid < self.sources.len());
//                 return (true, mid);
//             }

//             size = right - left;
//         }

//         assert!(left <= self.sources.len());
//         (false, left)
// 	}

//     // fn read_number(json: &Value, field: &str) -> Result<usize> {
//     //     match &json[field] {
//     //         Value::Null => return Err(anyhow!("field not found")),
//     //         Value::Bool(_) => return Err(anyhow!("field is a boolean, where number expected")),
//     //         Value::Number(n) => {
//     //             if !n.is_u64() {
//     //                 return Err(anyhow!("number is not an integer"))
//     //             }
//     //             return Ok(usize::try_from(n.as_u64().unwrap())?);
//     //         },
//     //         Value::String(_) => return Err(anyhow!("field is a literal, where number expected")),
//     //         Value::Array(_) => return Err(anyhow!("field is a list, where number expected")),
//     //         Value::Object(_) => return Err(anyhow!("field is an object, where number expected")),
//     //     }
//     // }

//     // fn read_list<'a>(json: &'a Value, field: &str) -> Result<&'a Vec<Value>> {
//     //     match &json[field] {
//     //         Value::Null => return Err(anyhow!("field not found")),
//     //         Value::Bool(_) => return Err(anyhow!("field is a boolean, where list expected")),
//     //         Value::Number(_) => return Err(anyhow!("field is a number, where list expected")),
//     //         Value::String(_) => return Err(anyhow!("field is a literal, where list expected")),
//     //         Value::Array(arr) => return Ok(&arr),
//     //         Value::Object(_) => return Err(anyhow!("field is an object, where list expected")),
//     //     }
//     // }

//     // fn read_string<'a>(json: &'a Value, field: &str) -> Result<String> {
//     //     match &json[field] {
//     //         Value::Null => return Err(anyhow!("field not found")),
//     //         Value::Bool(_) => return Err(anyhow!("field is a boolean, where literal expected")),
//     //         Value::Number(n) => return Ok(n.to_string()),
//     //         Value::String(s) => return Ok(s.to_string()),
//     //         Value::Array(_) => return Err(anyhow!("field is a list, where literal expected")),
//     //         Value::Object(_) => return Err(anyhow!("field is an object, where literal expected")),
//     //     }
//     // }

//     pub fn set_activity_key(&mut self, activity_key: &ActivityKey) {
//         self.activity_key = activity_key.clone();
//     }
// }


// impl Infoable for StochasticTransitionSystem {
//     fn info(&self, f: &mut impl std::io::Write) -> Result<()> {
//         writeln!(f, "Number of states\t{}", self.max_state)?;
//         writeln!(f, "Number of transitions\t{}", self.sources.len())?;
//         writeln!(f, "Number of activities\t{}", self.activity_key.get_number_of_activities())?;

//         Ok(write!(f, "")?)
//     }
// }

// impl fmt::Display for StochasticTransitionSystem {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         writeln!(f, "{{")?;
//         writeln!(f, "\"initialState\": {},", self.get_initial_state())?;
//         writeln!(f, "\"transitions\": [")?;
//         for pos in 0..self.sources.len() {
//             write!(f, "{{\"from\":{},\"to\":{},\"label\":\"{}\",\"freq\":\"{}\"}}", 
//                 self.sources[pos], 
//                 self.targets[pos], 
//                 self.activity_key.get_activity_label(&self.activities[pos]),
//                 self.transition_frequencies[pos])?;
//             if pos + 1 == self.sources.len() {
//                 writeln!(f, "")?;
//             } else {
//                 writeln!(f, ",")?;
//             }
//         }
//         writeln!(f, "]}}")?;
//         Ok(())
//     }
// }

// impl EbiTraitGraphable for StochasticTransitionSystem {
//     fn to_dot(&self) -> Result<layout::topo::layout::VisualGraph> {
//         log::info!("to_dot for StochasticTransitionSystem");
//         let mut graph = VisualGraph::new(layout::core::base::Orientation::LeftToRight);

//         // Create a place for each state
//         let mut places = vec![];
//         for state in 0..=self.max_state {
//             places.push(<dyn EbiTraitGraphable>::create_place(
//                 &mut graph,
//                 &format!("{}", self.state_terminating_frequencies[state]),
//             ));
//         }

//         // Add transitions as edges between places
//         for pos in 0..self.sources.len() {
//             let source = places[self.sources[pos]];
//             let target = places[self.targets[pos]];
//             let frequency = &self.transition_frequencies[pos];
//             let activity = self.activity_key.get_activity_label(&self.activities[pos]);

//             <dyn EbiTraitGraphable>::create_edge(
//                 &mut graph,
//                 &source,
//                 &target,
//                 &format!("{}, {}", activity, frequency),
//             );
//         }

//         Ok(graph)
//     }
// }
   
//     // fn to_dot(&self) -> layout::topo::layout::VisualGraph {
//     //     let mut graph = VisualGraph::new(layout::core::base::Orientation::LeftToRight);
//     //     let mut set = HashSet::new();
    
//     //     // Insert all elements from both vectors
//     //     for item in self.targets.clone() {
//     //         set.insert(item.clone());
//     //     }
//     //     for item in self.sources.clone() {
//     //         set.insert(item.clone());
//     //     }
//     //     // Convert back to vector
//     //     let state_set:Vec<usize> = set.into_iter().collect();

//     //     let mut places: Vec<layout::adt::dag::NodeHandle> = vec![];
//     //     for state in state_set.clone() {
//     //         println!("state: {}", state);
//     //         places.push(<dyn Dottable>::create_place(&mut graph, &format!("{}", self.state_terminating_frequencies[state])));
//     //     }

//     //     println!("sources: {:?}", self.sources);
//     //     println!("targets: {:?}", self.targets);
//     //     for pos in 0..self.sources.len() {
//     //         println!("pos: {}", pos);
//     //         println!("source: {}", self.sources[pos]);
//     //         println!("target: {}", self.targets[pos]);
//     //         let source = places[self.sources[pos]];
//     //         let target = places[self.targets[pos]];
//     //         let frequency = &self.transition_frequencies[pos];
//     //         let activity = self.activity_key.get_activity_label(&self.activities[pos]);
            
//     //         <dyn Dottable>::create_edge(&mut graph, &source, &target, &format!("{}, {}", activity, frequency.to_string()));
//     //     }

//     //     return graph;
//     // }


// impl<'a> IntoIterator for &'a StochasticTransitionSystem {
//     type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

//     type IntoIter = StochasticTransitionSystemIterator<'a>;

//     fn into_iter(self) -> Self::IntoIter {
//         Self::IntoIter {
//             it_sources:  self.sources.iter(),
//             it_targets:  self.targets.iter(),
//             it_activities: self.activities.iter(),
//             it_transition_frequencies: self.transition_frequencies.iter()
//         }
//     }
// }

// pub struct StochasticTransitionSystemIterator<'a> {
//     it_sources: std::slice::Iter<'a, usize>,
//     it_targets: std::slice::Iter<'a, usize>,
//     it_activities: std::slice::Iter<'a, Activity>,
//     it_transition_frequencies: std::slice::Iter<'a, usize>
// }

// impl<'a> Iterator for StochasticTransitionSystemIterator<'a> {
//     type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

//     fn next(&mut self) -> Option<Self::Item> {
//         if let Some(source) = self.it_sources.next() {
//             let target = self.it_targets.next().unwrap();
//             let activity = self.it_activities.next().unwrap();
//             let frequency = self.it_transition_frequencies.next().unwrap();
//             let result = Some((source, target, activity, frequency));
//             result
//         } else {
//             None
//         }
//     }
// }

// impl<'a> IntoIterator for &'a mut StochasticTransitionSystem {
//     type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

//     type IntoIter = StochasticTransitionSystemMutIterator<'a>;

//     fn into_iter(self) -> Self::IntoIter {
//         Self::IntoIter {
//             it_sources:  self.sources.iter(),
//             it_targets:  self.targets.iter(),
//             it_activities: self.activities.iter(),
//             it_transition_frequencies: self.transition_frequencies.iter()
//         }
//     }
// }

// pub struct StochasticTransitionSystemMutIterator<'a> {
//     it_sources: std::slice::Iter<'a, usize>,
//     it_targets: std::slice::Iter<'a, usize>,
//     it_activities: std::slice::Iter<'a, Activity>,
//     it_transition_frequencies: std::slice::Iter<'a, usize>
// }

// impl<'a> Iterator for StochasticTransitionSystemMutIterator<'a> {
//     type Item = (&'a usize, &'a usize, &'a Activity, &'a usize);

//     fn next(&mut self) -> Option<Self::Item> {
//         if let Some(source) = self.it_sources.next() {
//             let target = self.it_targets.next().unwrap();
//             let activity = self.it_activities.next().unwrap();
//             let frequency = self.it_transition_frequencies.next().unwrap();
//             let result = Some((source, target, activity, frequency));
//             result
//         } else {
//             None
//         }
//     }
// }

// impl Clone for StochasticTransitionSystem {
//     fn clone(&self) -> Self {
//         Self {
//             activity_key: self.activity_key.clone(),
//             initial_state: self.initial_state,
//             max_state: self.max_state,
//             sources: self.sources.clone(),
//             targets: self.targets.clone(),
//             activities: self.activities.clone(),
//             transition_frequencies: self.transition_frequencies.clone(),
//             state_terminating_frequencies: self.state_terminating_frequencies.clone(),
//         }
//     }
// }