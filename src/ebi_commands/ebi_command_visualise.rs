use layout::backends::svg::SVGWriter;

use crate::ebi_framework::{dottable::Dottable, ebi_command::EbiCommand, ebi_input::{EbiInput, EbiInputType}, ebi_object::{EbiObject, EbiObjectType}, ebi_output::{EbiOutput, EbiOutputType}};

pub const EBI_VISUALISE: EbiCommand = EbiCommand::Group { 
    name_short: "vis", 
    name_long: Some("visualise"),
    explanation_short: "Visualse an object.", 
    explanation_long: None,
    children: &[
        &EBI_VISUALISE_SVG,
        &EBI_VISUALISE_TEXT
    ]
};

pub const EBI_VISUALISE_TEXT: EbiCommand = EbiCommand::Command { 
    name_short: "txt", 
    name_long: Some("text"),
    explanation_short: "Visualise an object as text.",
    explanation_long: None, 
    latex_link: None, 
    cli_command: None, 
    exact_arithmetic: false, 
    input_types: &[ &[&EbiInputType::AnyObject] ], 
    input_names: &[ "FILE" ], 
    input_helps: &[ "Any file that can be visualised textually." ], 
    execute: |mut inputs, _| {
        let result = match inputs.remove(0) {
                EbiInput::Object(EbiObject::StochasticLabelledPetriNet(slpn), _) => slpn.to_string(),
                EbiInput::Object(EbiObject::LabelledPetriNet(lpn), _) => lpn.to_string(),
                EbiInput::Object(EbiObject::FiniteStochasticLanguage(lang), _) => lang.to_string(),
                EbiInput::Object(EbiObject::StochasticDeterministicFiniteAutomaton(sdfa), _) => sdfa.to_string(),
                EbiInput::Object(EbiObject::EventLog(log), _) => log.to_string(),
                EbiInput::Object(EbiObject::FiniteLanguage(language), _) => language.to_string(),
                EbiInput::Object(EbiObject::DirectlyFollowsModel(d), _) => d.to_string(),
                EbiInput::Object(EbiObject::Alignments(a), _) => a.to_string(),
                EbiInput::Object(EbiObject::DeterministicFiniteAutomaton(s), _) => s.to_string(),
                EbiInput::Object(EbiObject::Executions(s), _) => s.to_string(),
                EbiInput::FileHandler(_) => unreachable!(),
                EbiInput::Trait(_, _) => unreachable!(),
                EbiInput::String(_) => unreachable!(),
                EbiInput::Usize(_) => unreachable!(),
                EbiInput::Fraction(_) => unreachable!(),
        };
        Ok(EbiOutput::String(result))
    }, 
    output_type: &EbiOutputType::String
};

pub const EBI_VISUALISE_SVG: EbiCommand = EbiCommand::Command { 
    name_short: "svg", 
    name_long: None, 
    explanation_short: "Visualise an object as scalable vector graphics.",
    explanation_long: None, 
    latex_link: None, 
    cli_command: None, 
    exact_arithmetic: true, 
    input_types: &[ 
        &[
            &EbiInputType::Object(EbiObjectType::StochasticLabelledPetriNet),
            &EbiInputType::Object(EbiObjectType::StochasticDeterministicFiniteAutomaton),
            &EbiInputType::Object(EbiObjectType::LabelledPetriNet),
            &EbiInputType::Object(EbiObjectType::DirectlyFollowsModel),
            &EbiInputType::Object(EbiObjectType::DeterministicFiniteAutomaton),
        ] ], 
    input_names: &[ "FILE" ], 
    input_helps: &[ "Any file that can be visualised as a graph." ], 
    execute: |mut inputs, _| {
        let mut result = match inputs.remove(0) {
            EbiInput::Object(EbiObject::LabelledPetriNet(lpn), _) => lpn.to_dot(),
            EbiInput::Object(EbiObject::StochasticLabelledPetriNet(slpn), _) => slpn.to_dot(),
            EbiInput::Object(EbiObject::StochasticDeterministicFiniteAutomaton(sdfa), _) => sdfa.to_dot(),
            EbiInput::Object(EbiObject::DirectlyFollowsModel(dfm), _) => dfm.to_dot(),
            EbiInput::Object(EbiObject::EventLog(_), _) => unreachable!(),
            EbiInput::Object(EbiObject::FiniteLanguage(_), _) => unreachable!(),
            EbiInput::Object(EbiObject::FiniteStochasticLanguage(_), _) => unreachable!(),
            EbiInput::Object(EbiObject::Alignments(_), _) => unreachable!(),
            EbiInput::Object(EbiObject::DeterministicFiniteAutomaton(dfa), _) => dfa.to_dot(),
            EbiInput::Object(EbiObject::Executions(_), _) => unreachable!(),
            EbiInput::FileHandler(_) => unreachable!(),
            EbiInput::Trait(_, _) => unreachable!(),
            EbiInput::String(_) => unreachable!(),
            EbiInput::Usize(_) => unreachable!(),
            EbiInput::Fraction(_) => unreachable!(),
        };

        let mut svg = SVGWriter::new();
        result.do_it(false, false, false, &mut svg);

        return Ok(EbiOutput::SVG(svg.finalize()));
    
    }, 
    output_type: &EbiOutputType::SVG
};