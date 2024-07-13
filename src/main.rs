extern crate libc;

use clap::Parser;

use kclassify::build;

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
}

fn main() {
    env_logger::init_from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    );

    // RUST annoyingly panicked when piping results to `head` (which isn't technically an error)
    // so we need to disable this.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    let args = Args::parse();
    match args.subcommand {
        SubCommand::Build(build) => {
            build::build(&build.input);
        }
    }
}
