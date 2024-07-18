use clap::Parser;
use libc::{SIGPIPE, SIG_DFL};

use klassify::breakpoint;
use klassify::build;
use klassify::classify;
use klassify::extract;
use klassify::extract_bam;
use klassify::info;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
struct Args {
    #[clap(subcommand)]
    subcommand: SubCommand,
}

#[derive(Parser, Debug)]
enum SubCommand {
    #[clap(about = "Build reference kmer table")]
    Build(build::BuildArgs),
    #[clap(about = "Classify reads")]
    Classify(classify::ClassifyArgs),
    #[clap(about = "Print details about the kmer table")]
    Info(info::InfoArgs),
    #[clap(about = "Detect breakpoints")]
    Breakpoint(breakpoint::BreakpointArgs),
    #[clap(about = "Extract reads")]
    Extract(extract::ExtractArgs),
    #[clap(about = "Extract reads from BAM")]
    ExtractBam(extract_bam::ExtractBamArgs),
}

fn main() {
    env_logger::init_from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    );

    // RUST annoyingly panicked when piping results to `head` (which isn't technically an error)
    // so we need to disable this.
    unsafe {
        libc::signal(SIGPIPE, SIG_DFL);
    }

    let args = Args::parse();
    match args.subcommand {
        SubCommand::Build(build) => {
            build::build(&build.fasta_files);
        }
        SubCommand::Classify(classify) => {
            classify::classify(&classify.bincode_file, &classify.reads_file);
        }
        SubCommand::Info(info) => {
            info::info(&info.bincode_file);
        }
        SubCommand::Breakpoint(breakpoint) => {
            breakpoint::breakpoint(&breakpoint.bincode_file, &breakpoint.fasta_files);
        }
        SubCommand::Extract(extract) => {
            extract::extract(&extract.reads_tsv, &extract.fasta_files);
        }
        SubCommand::ExtractBam(extract_bam) => {
            extract_bam::extract_bam(&extract_bam.regions_file, &extract_bam.bam_file);
        }
    }
}
