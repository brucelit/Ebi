use anyhow::{Context, anyhow};
use clap::{Arg, ArgAction, value_parser};
use ebi_objects::HasActivityKey;
use ebi_arithmetic::{Fraction, Zero,fraction::{fraction_enum::FractionEnum}};
use crate::{
    ebi_framework::{
        ebi_command::EbiCommand,
        ebi_input::EbiInputType,
        ebi_output::{EbiOutput, EbiOutputType},
        ebi_trait::EbiTrait,
    },
    ebi_traits::{
        ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage,
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
        &EBI_PREDICTION_TRACE_NEXT,
        &EBI_PREDICTION_LOG_NEXT,
        // &EBI_PREDICTION_SUFFIX,
    ],
};

pub const EBI_PREDICTION_TRACE_NEXT: EbiCommand = EbiCommand::Command {
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

            let num = semantics
                .predict_trace(&trace)
                .with_context(|| format!("cannot explain the trace {:?}", trace))?;
            Ok(EbiOutput::Fraction(FractionEnum::from((num, trace.len()))))
        } else {
            Err(anyhow!("no trace given"))
        }
    },
    output_type: &EbiOutputType::Fraction,
};

pub const EBI_PREDICTION_LOG_NEXT: EbiCommand = EbiCommand::Command {
    name_short: "prelog",
    name_long: Some("predict-log"),
    library_name: "ebi_commands::ebi_command_prediction::EBI_PREDICT_LOG",
    explanation_short: "Compute the most likely next activity of a log given the stochastic model.",
    explanation_long: None,
    latex_link: None,
    cli_command: None,
    exact_arithmetic: false,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Trait(EbiTrait::StochasticSemantics)],
    ],
    input_names: &["FILE_1", "FILE_2"],
    input_helps: &[
        "The model.",
        "The event log",
    ],
    execute: |mut inputs, _| {
        let log = inputs
        .remove(0)
        .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;

        let semantics = inputs.remove(0).to_type::<EbiTraitStochasticSemantics>()?;

        let _ = semantics
            .predict_log(log)
            .context("cannot make prediction for the log")?;
        Ok(EbiOutput::Fraction(Fraction::zero()))
    },
    output_type: &EbiOutputType::Fraction,
};








