use std::{collections::BTreeSet, fmt::{self, Display}, fs::File, io::Write, path::PathBuf};
use anyhow::{Context, Result};

use crate::{ebi_objects::{alignments::{Alignments, EBI_ALIGNMENTS}, compressed_event_log::{CompressedEventLog, EBI_COMPRESSED_EVENT_LOG}, deterministic_finite_automaton::{DeterministicFiniteAutomaton, EBI_DETERMINISTIC_FINITE_AUTOMATON}, directly_follows_model::{DirectlyFollowsModel, EBI_DIRCTLY_FOLLOWS_MODEL}, executions::{Executions, EBI_EXECUTIONS}, finite_language::{FiniteLanguage, EBI_FINITE_LANGUAGE}, finite_stochastic_language::{FiniteStochasticLanguage, EBI_FINITE_STOCHASTIC_LANGUAGE}, labelled_petri_net::{LabelledPetriNet, EBI_LABELLED_PETRI_NET}, stochastic_deterministic_finite_automaton::{StochasticDeterministicFiniteAutomaton, EBI_STOCHASTIC_DETERMINISTIC_FINITE_AUTOMATON}, stochastic_labelled_petri_net::{StochasticLabelledPetriNet, EBI_STOCHASTIC_LABELLED_PETRI_NET}}, math::{fraction::Fraction, log_div::LogDiv, root::ContainsRoot, root_log_div::RootLogDiv}};

use super::{ebi_command::{EbiCommand, EBI_COMMANDS}, ebi_file_handler::{EbiFileHandler, EBI_FILE_HANDLERS}, ebi_object::{EbiObject, EbiObjectType}, exportable::Exportable, prom_link::{JavaObjectHandler, JAVA_OBJECT_HANDLERS_CONTAINSROOT, JAVA_OBJECT_HANDLERS_FRACTION, JAVA_OBJECT_HANDLERS_LOGDIV, JAVA_OBJECT_HANDLERS_ROOTLOGDIV, JAVA_OBJECT_HANDLERS_STRING, JAVA_OBJECT_HANDLERS_SVG, JAVA_OBJECT_HANDLERS_USIZE}};

pub enum EbiOutput {
    Object(EbiObject),
    String(String),
    SVG(String),
    Usize(usize),
    Fraction(Fraction),
    LogDiv(LogDiv),
    ContainsRoot(ContainsRoot),
    RootLogDiv(RootLogDiv)
}

impl EbiOutput {
    pub fn get_type(&self) -> EbiOutputType {
        match self {
            EbiOutput::Object(o) => EbiOutputType::ObjectType(o.get_type()),
            EbiOutput::String(_) => EbiOutputType::String,
            EbiOutput::SVG(_) => EbiOutputType::SVG,
            EbiOutput::Usize(_) => EbiOutputType::Usize,
            EbiOutput::Fraction(_) => EbiOutputType::Fraction,
            EbiOutput::LogDiv(_) => EbiOutputType::LogDiv,
            EbiOutput::ContainsRoot(_) => EbiOutputType::ContainsRoot,
            EbiOutput::RootLogDiv(_) => EbiOutputType::RootLogDiv
        }
    }
}

#[derive(PartialEq,Eq)]
pub enum EbiOutputType {
    ObjectType(EbiObjectType),
    String,
    SVG,
    Usize,
    Fraction,
    LogDiv,
    ContainsRoot,
    RootLogDiv
}

impl EbiOutputType {

    /**
     * Get all commands that output this type.
     */
    pub fn get_applicable_commands(&self) -> BTreeSet<Vec<&'static EbiCommand>> {
        let mut result = EBI_COMMANDS.get_command_paths();
        result.retain(|path| {
            if let EbiCommand::Command { output_type: output, .. } = path[path.len() - 1] {
                if output == &self {
                    return true;
                }
            }
            false
        });
        result
    }

    /**
     * Returns all exporters that can handle this output type
     */
    pub fn get_exporters(&self) -> Vec<EbiExporter> {
        match self {
            EbiOutputType::ObjectType(etype) => {
                let mut result = vec![];
                for file_handler in EBI_FILE_HANDLERS {
                    for exporter in file_handler.object_exporters {
                        if &exporter.get_type() == etype {
                            result.push(EbiExporter::Object(exporter, file_handler))
                        }
                    }
                }
                result
            } ,
            EbiOutputType::String => vec![EbiExporter::String],
            EbiOutputType::SVG => vec![EbiExporter::SVG],
            EbiOutputType::Usize => vec![EbiExporter::Usize],
            EbiOutputType::Fraction => vec![EbiExporter::Fraction],
            EbiOutputType::LogDiv => vec![EbiExporter::LogDiv],
            EbiOutputType::ContainsRoot => vec![EbiExporter::ContainsRoot],
            EbiOutputType::RootLogDiv => vec![EbiExporter::RootLogDiv],
        }
    }

    pub fn get_default_exporter(&self) -> EbiExporter {
        match self {
            EbiOutputType::ObjectType(EbiObjectType::Alignments) => EbiExporter::Object(&EbiObjectExporter::Alignments(Alignments::export_from_object), &EBI_ALIGNMENTS),
            EbiOutputType::ObjectType(EbiObjectType::DeterministicFiniteAutomaton) => EbiExporter::Object(&EbiObjectExporter::DeterministicFiniteAutomaton(DeterministicFiniteAutomaton::export_from_object), &EBI_DETERMINISTIC_FINITE_AUTOMATON),
            EbiOutputType::ObjectType(EbiObjectType::DirectlyFollowsModel) => EbiExporter::Object(&EbiObjectExporter::DirectlyFollowsModel(DirectlyFollowsModel::export_from_object), &EBI_DIRCTLY_FOLLOWS_MODEL),
            EbiOutputType::ObjectType(EbiObjectType::EventLog) => EbiExporter::Object(&EbiObjectExporter::EventLog(CompressedEventLog::export_from_object), &EBI_COMPRESSED_EVENT_LOG),
            EbiOutputType::ObjectType(EbiObjectType::Executions) => EbiExporter::Object(&EbiObjectExporter::Executions(Executions::export_from_object), &EBI_EXECUTIONS),
            EbiOutputType::ObjectType(EbiObjectType::FiniteLanguage) => EbiExporter::Object(&EbiObjectExporter::FiniteLanguage(FiniteLanguage::export_from_object), &EBI_FINITE_LANGUAGE),
            EbiOutputType::ObjectType(EbiObjectType::FiniteStochasticLanguage) => EbiExporter::Object(&EbiObjectExporter::FiniteStochasticLanguage(FiniteStochasticLanguage::export_from_object), &EBI_FINITE_STOCHASTIC_LANGUAGE),
            EbiOutputType::ObjectType(EbiObjectType::LabelledPetriNet) => EbiExporter::Object(&EbiObjectExporter::LabelledPetriNet(LabelledPetriNet::export_from_object), &EBI_LABELLED_PETRI_NET),
            EbiOutputType::ObjectType(EbiObjectType::StochasticDeterministicFiniteAutomaton) => EbiExporter::Object(&&EbiObjectExporter::StochasticDeterministicFiniteAutomaton(StochasticDeterministicFiniteAutomaton::export_from_object), &EBI_STOCHASTIC_DETERMINISTIC_FINITE_AUTOMATON),
            EbiOutputType::ObjectType(EbiObjectType::StochasticLabelledPetriNet) => EbiExporter::Object(&&EbiObjectExporter::StochasticLabelledPetriNet(StochasticLabelledPetriNet::export_from_object), &EBI_STOCHASTIC_LABELLED_PETRI_NET),
            EbiOutputType::String => EbiExporter::String,
            EbiOutputType::SVG => EbiExporter::SVG,
            EbiOutputType::Usize => EbiExporter::Usize,
            EbiOutputType::Fraction => EbiExporter::Fraction,
            EbiOutputType::LogDiv => EbiExporter::LogDiv,
            EbiOutputType::ContainsRoot => EbiExporter::ContainsRoot,
            EbiOutputType::RootLogDiv => EbiExporter::RootLogDiv,
        }
    }

    pub fn exporters_as_strings_with_articles(&self, last_connector: &str) -> String {
        let mut list = self.get_exporters().into_iter().map(|exp| exp.get_article().to_string() + " " + &exp.to_string()).collect::<Vec<_>>();
        if list.len() == 1 {
            return list.remove(0)
        }
        let (last, list) = list.split_last().unwrap();
        format!("{} {} {}", list.join(", "), last_connector, last)
    }
}

impl Display for EbiOutputType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EbiOutputType::ObjectType(t) => t.fmt(f),
            EbiOutputType::String => Display::fmt(&"text", f),
            EbiOutputType::SVG => Display::fmt(&"svg", f),
            EbiOutputType::Usize => Display::fmt(&"integer", f),
            EbiOutputType::Fraction => Display::fmt(&"fraction", f),
            EbiOutputType::LogDiv => Display::fmt(&"logarithm", f),
            EbiOutputType::ContainsRoot => Display::fmt(&"root", f),
            EbiOutputType::RootLogDiv => Display::fmt(&"rootlog", f),
        }
    }
}

#[derive(Debug,Clone,Eq,PartialEq,Hash)]
pub enum EbiExporter {
    Object(&'static EbiObjectExporter, &'static EbiFileHandler),
    String,
    SVG,
    Usize,
    Fraction,
    LogDiv,
    ContainsRoot,
    RootLogDiv
}

impl EbiExporter {
    pub fn export_from_object(&self, output: EbiOutput, f: &mut dyn std::io::Write) -> Result<()> {
        match (self, output) {
            (EbiExporter::Object(exporter, _), object) => exporter.export(object, f),
            (EbiExporter::String, EbiOutput::String(object)) => object.export(f),
            (EbiExporter::String, _) => unreachable!(),
            (EbiExporter::SVG, EbiOutput::SVG(object)) => object.export(f),
            (EbiExporter::SVG, _) => unreachable!(),
            (EbiExporter::Usize, EbiOutput::Usize(object)) => object.export(f),
            (EbiExporter::Usize, _) => unreachable!(),
            (EbiExporter::Fraction, EbiOutput::Fraction(object)) => object.export(f),
            (EbiExporter::Fraction, _) => unreachable!(),
            (EbiExporter::LogDiv, EbiOutput::LogDiv(object)) => object.export(f),
            (EbiExporter::LogDiv, _) => unreachable!(),
            (EbiExporter::ContainsRoot, EbiOutput::ContainsRoot(object)) => object.export(f),
            (EbiExporter::ContainsRoot, _) => unreachable!(),
            (EbiExporter::RootLogDiv, EbiOutput::RootLogDiv(object)) => object.export(f),
            (EbiExporter::RootLogDiv, _) => unreachable!(),
        }
    }

    pub fn get_article(&self) -> &str {
        match self {
            EbiExporter::Object(_, file_handler) => file_handler.article,
            EbiExporter::String => "",
            EbiExporter::SVG => "an",
            EbiExporter::Usize => "an",
            EbiExporter::Fraction => "a",
            EbiExporter::LogDiv => "a",
            EbiExporter::ContainsRoot => "a",
            EbiExporter::RootLogDiv => "a",
        }
    }

    pub fn get_name(&self) -> &str  {
        match self {
            EbiExporter::Object(_, file_handler) => file_handler.name,
            EbiExporter::String => "string",
            EbiExporter::SVG => "SVG",
            EbiExporter::Usize => "integer",
            EbiExporter::Fraction => "fraction",
            EbiExporter::LogDiv => "logdiv",
            EbiExporter::ContainsRoot => "containsroot",
            EbiExporter::RootLogDiv => "rootlogdiv",
        }
    }

    pub fn get_java_object_handlers(&self) -> &'static [JavaObjectHandler] {
        match self {
            EbiExporter::Object(_, file_handler) => file_handler.java_object_handlers,
            EbiExporter::String => JAVA_OBJECT_HANDLERS_STRING,
            EbiExporter::SVG => JAVA_OBJECT_HANDLERS_SVG,
            EbiExporter::Usize => JAVA_OBJECT_HANDLERS_USIZE,
            EbiExporter::Fraction => JAVA_OBJECT_HANDLERS_FRACTION,
            EbiExporter::LogDiv => JAVA_OBJECT_HANDLERS_LOGDIV,
            EbiExporter::ContainsRoot => JAVA_OBJECT_HANDLERS_CONTAINSROOT,
            EbiExporter::RootLogDiv => JAVA_OBJECT_HANDLERS_ROOTLOGDIV,
        }
    }

    pub fn get_extension(&self) -> &str {
        match self {
            EbiExporter::Object(_, file_handler) => file_handler.file_extension,
            EbiExporter::String => "txt",
            EbiExporter::SVG => "svg",
            EbiExporter::Usize => "int",
            EbiExporter::Fraction => "frac",
            EbiExporter::LogDiv => "logdiv",
            EbiExporter::ContainsRoot => "croot",
            EbiExporter::RootLogDiv => "rldiv",
        }
    }
}

impl Display for EbiExporter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EbiExporter::Object(_, file_handler) => Display::fmt(file_handler, f),
            EbiExporter::String => Display::fmt(&"text", f),
            EbiExporter::SVG => Display::fmt(&"svg", f),
            EbiExporter::Usize => Display::fmt(&"integer", f),
            EbiExporter::Fraction => Display::fmt(&"fraction", f),
            EbiExporter::LogDiv => Display::fmt(&"logarithm", f),
            EbiExporter::ContainsRoot => Display::fmt(&"root", f),
            EbiExporter::RootLogDiv => Display::fmt(&"rootlog", f)
        }
    }
}

#[derive(Debug,Eq,PartialEq,Hash)]
pub enum EbiObjectExporter {
    EventLog(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    DirectlyFollowsModel(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    FiniteLanguage(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    FiniteStochasticLanguage(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    LabelledPetriNet(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    StochasticDeterministicFiniteAutomaton(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    StochasticLabelledPetriNet(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    Alignments(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    DeterministicFiniteAutomaton(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
    Executions(fn(object: EbiOutput, &mut dyn std::io::Write) -> Result<()>),
}

impl EbiObjectExporter {
    pub fn get_type(&self) -> EbiObjectType {
        match self {
            EbiObjectExporter::EventLog(_) => EbiObjectType::EventLog,
            EbiObjectExporter::DirectlyFollowsModel(_) => EbiObjectType::DirectlyFollowsModel,
            EbiObjectExporter::FiniteLanguage(_) => EbiObjectType::FiniteLanguage,
            EbiObjectExporter::FiniteStochasticLanguage(_) => EbiObjectType::FiniteStochasticLanguage,
            EbiObjectExporter::LabelledPetriNet(_) => EbiObjectType::LabelledPetriNet,
            EbiObjectExporter::StochasticDeterministicFiniteAutomaton(_) => EbiObjectType::StochasticDeterministicFiniteAutomaton,
            EbiObjectExporter::StochasticLabelledPetriNet(_) => EbiObjectType::StochasticLabelledPetriNet,
            EbiObjectExporter::Alignments(_) => EbiObjectType::Alignments,
            EbiObjectExporter::DeterministicFiniteAutomaton(_) => EbiObjectType::DeterministicFiniteAutomaton,
            EbiObjectExporter::Executions(_) => EbiObjectType::Executions,
        }
    }

    pub fn export(&self, object: EbiOutput, f: &mut dyn std::io::Write) -> Result<()> {
        match self {
            EbiObjectExporter::EventLog(exporter) => (exporter)(object, f),
            EbiObjectExporter::DirectlyFollowsModel(exporter) => (exporter)(object, f),
            EbiObjectExporter::FiniteLanguage(exporter) => (exporter)(object, f),
            EbiObjectExporter::FiniteStochasticLanguage(exporter) => (exporter)(object, f),
            EbiObjectExporter::LabelledPetriNet(exporter) => (exporter)(object, f),
            EbiObjectExporter::StochasticDeterministicFiniteAutomaton(exporter) => (exporter)(object, f),
            EbiObjectExporter::StochasticLabelledPetriNet(exporter) => (exporter)(object, f),
            EbiObjectExporter::Alignments(exporter) => (exporter)(object, f),
            EbiObjectExporter::DeterministicFiniteAutomaton(exporter) => (exporter)(object, f),
            EbiObjectExporter::Executions(exporter) => (exporter)(object, f),
        }
    }
}

pub fn export_object(to_file: &PathBuf, object: EbiOutput, exporter: EbiExporter) -> Result<()> {
    let file = File::create(to_file).with_context(|| format!("Writing result to file {:?}.", to_file))?;
    let mut writer = std::io::BufWriter::new(&file);
    exporter.export_from_object(object, &mut writer).with_context(|| format!("Writing result to file {:?}.", to_file))?;
    return writer.flush().with_context(|| format!("writing result to file {:?}", to_file));
}

pub fn export_to_string(object: EbiOutput, exporter: EbiExporter) -> Result<String> {
    let mut f = vec![];
    exporter.export_from_object(object, &mut f)?;
    Ok(String::from_utf8(f)?)
}

impl Display for EbiObjectExporter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.get_type().to_string())
    }
}