use anyhow::Context;
use clap::{Arg, ArgAction, Command, value_parser};
use std::io::Write;

use crate::{
    ebi_framework::{
        ebi_command::EbiCommand,
        ebi_input::EbiInputType,
        ebi_output::{EbiOutput, EbiOutputType},
        ebi_trait::EbiTrait,
    },
    ebi_info,
    ebi_traits::{
        ebi_trait_event_log::EbiTraitEventLog,
        ebi_trait_finite_stochastic_language::EbiTraitFiniteStochasticLanguage,
    },
    math::fraction::{Fraction, FractionNotParsedYet},
    techniques::{permutation_test::PermutationTest, bootstrap_test::{
        BootstrapTest, StatisticalTestsLogCategoricalAttribute,
    }},
};

use super::ebi_command_association::{self, number_of_samples};

macro_rules! p_value {
    () => {
        Fraction::from((1usize, 20usize))
    };
}
pub(crate) use p_value;

pub const EBI_TEST: EbiCommand = EbiCommand::Group {
    name_short: "tst",
    name_long: Some("test"),
    explanation_short: "Test a hypothesis.",
    explanation_long: None,
    children: &[&EBI_BOOTSTRAP_TEST, &EBI_TEST_LOG_ATTRIBUTE],
};

pub const EBI_TEST_LOG_ATTRIBUTE: EbiCommand = EbiCommand::Command {
    name_short: "lcat",
    name_long: Some("log-categorical-attribute"),
    explanation_short: "Test the hypothesis that the sub-logs defined by the categorical attribute are derived from identical processes.",
    explanation_long: Some(concat!(
        "Test the hypothesis that the sub-logs defined by the categorical attribute are derived from identical processes.; ",
        number_of_samples!(),
        " samples are taken."
    )),
    latex_link: Some("\\cite{DBLP:journals/tkde/LeemansMPH23}"),
    cli_command: Some(|command| {
        cli_p_value(ebi_command_association::cli_number_of_samples(command))
    }),
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::EventLog)],
        &[&EbiInputType::String],
    ],
    input_names: &["FILE", "ATTRIBUTE"],
    input_helps: &[
        "The event log for which the test is to be performed.",
        concat!(
            concat!(
                "The trace attribute for which the test is to be performed. The trace attributes of a log can be found using `Ebi ",
                ebi_info!()
            ),
            "`."
        ),
    ],
    execute: |mut inputs, cli_matches| {
        let event_log = inputs.remove(0).to_type::<dyn EbiTraitEventLog>()?;
        let attribute = inputs.remove(0).to_type::<String>()?;
        let number_of_samples = cli_matches.unwrap().get_one::<usize>("samples").unwrap();
        let p_value = cli_matches
            .unwrap()
            .get_one::<FractionNotParsedYet>("pvalue")
            .unwrap()
            .try_into()
            .context("Parsing p value")?;

        let (value, sustained) = event_log
            .log_categorical_attribute(*number_of_samples, &attribute, &p_value)
            .with_context(|| format!("attribute {}", attribute))?;

        let mut f = vec![];
        writeln!(f, "p-value \t {}", value)?;
        writeln!(
            f,
            "Null-hypothesis: the sub-logs defined by trace attribute `{}` all follow identical processes.",
            attribute
        )?;

        if sustained {
            writeln!(
                f,
                "The data does not provide enough evidence to make a claim on this hypothesis.\nDo not reject the null-hypothesis."
            )?;
        } else {
            writeln!(
                f,
                "The data provides enough evidence to conclude that at least one value of the categorical attribute `{}` associates with a difference in process.\nReject the null-hypothesis.",
                attribute
            )?;
        };

        Ok(EbiOutput::String(String::from_utf8(f).unwrap()))
    },
    output_type: &EbiOutputType::String,
};

pub const EBI_BOOTSTRAP_TEST: EbiCommand = EbiCommand::Command {
    name_short: "btst",
    name_long: Some("bootstrap-test"),
    explanation_short: "Test the hypothesis that the logs are derived from identical processes.",
    explanation_long: Some(concat!(
        "Test the hypothesis that the logs are derived from identical processes; ",
        number_of_samples!(),
        " samples are taken."
    )),
    latex_link: Some("\\cite{DBLP:journals/tkde/LeemansMPH23}"),
    cli_command: Some(|command| {
        cli_p_value(ebi_command_association::cli_number_of_samples(command))
    }),
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
    ],
    input_names: &["LANG_1", "LANG_2"],
    input_helps: &[
        "The first event log for which the test is to be performed.",
        "The first event log for which the test is to be performed.",
    ],
    execute: |mut inputs, cli_matches| {
        let mut log1 = inputs
            .remove(0)
            .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let mut log2: Box<dyn EbiTraitFiniteStochasticLanguage + 'static> =
            inputs
                .remove(0)
                .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let number_of_samples = cli_matches.unwrap().get_one::<usize>("samples").unwrap();
        let p_value = cli_matches
            .unwrap()
            .get_one::<FractionNotParsedYet>("pvalue")
            .unwrap()
            .try_into()
            .context("Parsing p value")?;

        let (value, sustained) = log1
            .bootstrap_test(log2.as_mut(), *number_of_samples, &p_value)
            .with_context(|| format!("performing test"))?;

        let mut f = vec![];
        writeln!(f, "p-value \t {}", value)?;
        writeln!(f, "Null-hypothesis: the logs follow identical processes.")?;

        if sustained {
            writeln!(
                f,
                "The data does not provide enough evidence to make a claim on this hypothesis.\nDo not reject the null-hypothesis."
            )?;
        } else {
            writeln!(
                f,
                "The data provides enough evidence to conclude that the logs are derived from different processes.\nReject the null-hypothesis.",
            )?;
        };

        Ok(EbiOutput::String(String::from_utf8(f).unwrap()))
    },
    output_type: &EbiOutputType::String,
};

pub const EBI_PERMUTATION_TEST: EbiCommand = EbiCommand::Command {
    name_short: "perm",
    name_long: Some("permutation-test"),
    explanation_short: "Test the hypothesis that the logs are derived from identical processes.",
    explanation_long: Some(concat!(
        "Test the hypothesis that the logs are derived from identical processes; ",
        number_of_samples!(),
        " samples are taken."
    )),
    latex_link: None,
    cli_command: Some(|command| {
        cli_p_value(ebi_command_association::cli_number_of_samples(command))
    }),
    exact_arithmetic: true,
    input_types: &[
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
        &[&EbiInputType::Trait(EbiTrait::FiniteStochasticLanguage)],
    ],
    input_names: &["LANG_1", "LANG_2"],
    input_helps: &[
        "The first event log for which the test is to be performed.",
        "The first event log for which the test is to be performed.",
    ],
    execute: |mut inputs, cli_matches| {
        let mut log1 = inputs
            .remove(0)
            .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let mut log2: Box<dyn EbiTraitFiniteStochasticLanguage + 'static> =
            inputs
                .remove(0)
                .to_type::<dyn EbiTraitFiniteStochasticLanguage>()?;
        let number_of_samples = cli_matches.unwrap().get_one::<usize>("samples").unwrap();
        let p_value = cli_matches
            .unwrap()
            .get_one::<FractionNotParsedYet>("pvalue")
            .unwrap()
            .try_into()
            .context("Parsing p value")?;

        let (value, sustained) = log1
            .permutation_test(log2.as_mut(), *number_of_samples, &p_value)
            .with_context(|| format!("performing test"))?;

        let mut f = vec![];
        writeln!(f, "p-value \t {}", value)?;
        writeln!(f, "Null-hypothesis: the logs follow identical processes.")?;

        if sustained {
            writeln!(
                f,
                "The data does not provide enough evidence to make a claim on this hypothesis.\nDo not reject the null-hypothesis."
            )?;
        } else {
            writeln!(
                f,
                "The data provides enough evidence to conclude that the logs are derived from different processes.\nReject the null-hypothesis.",
            )?;
        };

        Ok(EbiOutput::String(String::from_utf8(f).unwrap()))
    },
    output_type: &EbiOutputType::String,
};

pub fn cli_p_value(command: Command) -> Command {
    command.arg(
        Arg::new("pvalue")
            .action(ArgAction::Set)
            .value_name("NUMBER")
            .short('p')
            .long("p-value")
            .help("Use threshold p-value.")
            .default_value(p_value!().to_string())
            .value_parser(value_parser!(FractionNotParsedYet))
            .required(false),
    )
}
