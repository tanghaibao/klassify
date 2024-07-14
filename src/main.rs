use clap::Parser;
use libc::{SIGPIPE, SIG_DFL};

use klassify::build;
use klassify::classify;
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
    }
}
