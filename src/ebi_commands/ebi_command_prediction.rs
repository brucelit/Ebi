use anyhow::{Context, anyhow};
use clap::{Arg, ArgAction, value_parser};
use ebi_objects::{EbiObject, EbiObjectType, HasActivityKey};

use crate::{
    ebi_framework::{
        ebi_command::EbiCommand,
        ebi_input::EbiInputType,
        ebi_output::{EbiOutput, EbiOutputType},
        ebi_trait::EbiTrait,
    },
    ebi_traits::{
        ebi_trait_stochastic_semantics::EbiTraitStochasticSemantics,
    },
    math::constant_fraction::ConstFraction,
    techniques::predict_trace::PredictTrace,
};

pub const EBI_PREDICTION: EbiCommand = EbiCommand::Group {
    name_short: "pred",
    name_long: Some("prediction"),
    explanation_short: "Predict the next sequence of events in a running trace.",
    explanation_long: None,
    children: &[
        &EBI_PREDICTION_NEXT,
        // &EBI_PREDICTION_SUFFIX,
    ],
};

pub const EBI_PREDICTION_NEXT: EbiCommand = EbiCommand::Command {
    name_short: "pretra",
    name_long: Some("predict-trace"),
    library_name: "ebi_commands::ebi_command_prediction::EBI_PREDICT_TRACE",
    explanation_short: "Compute the most likely next activity of a trace given the stochastic model.",
    explanation_long: None,
    latex_link: None,
    cli_command: Some(|command| {
        command.arg(
            Arg::new("trace")
                .action(ArgAction::Set)
                .value_name("TRACE")
                .help("The trace.")
                .required(true)
                .value_parser(value_parser!(String))
                .num_args(0..),
        )
    }),
    exact_arithmetic: false,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::StochasticSemantics)],
        &[&EbiInputType::Fraction(
            Some(ConstFraction::zero()),
            Some(ConstFraction::one()),
            None,
        )],
    ],
    input_names: &["FILE", "VALUE"],
    input_helps: &[
        "The model.",
        "Balance between 0 (=only consider deviations) to 1 (=only consider weight in the model)",
    ],
    execute: |mut inputs, cli_matches| {
        let mut semantics = inputs.remove(0).to_type::<EbiTraitStochasticSemantics>()?;
        // let balance = inputs.remove(0).to_type::<Fraction>()?;
        if let Some(x) = cli_matches.unwrap().get_many::<String>("trace") {
            let t: Vec<&String> = x.collect();
            let trace = t
                .into_iter()
                .map(|activity| activity.as_str())
                .collect::<Vec<_>>();
            let trace = semantics.activity_key_mut().process_trace_ref(&trace);

            log::trace!("predict the trace {:?} given the model", trace);

            let result = semantics
                .predict_trace(&trace)
                .with_context(|| format!("cannot explain the trace {:?}", trace))?;
            return Ok(EbiOutput::Object(EbiObject::LanguageOfAlignments(result)));
        } else {
            return Err(anyhow!("no trace given"));
        }
    },
    output_type: &EbiOutputType::ObjectType(EbiObjectType::LanguageOfAlignments),
};





