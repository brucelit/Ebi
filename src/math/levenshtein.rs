use crate::{activity_key::Activity, math::fraction::Fraction};

pub fn normalised(trace1: &Vec<Activity>, trace2: &Vec<Activity>) -> Fraction {
    let dist = strsim::generic_levenshtein(trace1, trace2);

    Fraction::from((dist, trace1.len().max(trace2.len())))
}

pub fn distance(trace1: &Vec<usize>, trace2: &Vec<usize>) -> usize {
    strsim::generic_levenshtein(trace1, trace2)
}