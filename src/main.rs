use clap::Parser;
use libc::{SIGPIPE, SIG_DFL};

use klassify::breakpoint;
use klassify::build;
use klassify::classify;
use klassify::extract;
use klassify::extract_bam;
use klassify::info;
use klassify::regions;
use klassify::sort_bam;

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
    #[clap(about = "Prepare BAM files and generate depths for each bin")]
    Regions(regions::RegionsArgs),
    #[clap(about = "Sort BAM file by divergence")]
    SortBam(sort_bam::SortBamArgs),
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
            build::build(&build.fasta_files, &build.output_file, build.kmer_size);
        }
        SubCommand::Classify(classify) => {
            classify::classify(
                &classify.bincode_file,
                &classify.reads_file,
                &classify.output_dir,
                classify.prefix_length,
            );
        }
        SubCommand::Info(info) => {
            info::info(&info.bincode_file);
        }
        SubCommand::Breakpoint(breakpoint) => {
            breakpoint::breakpoint(&breakpoint.bincode_file, &breakpoint.fasta_files);
        }
        SubCommand::Extract(extract) => {
            extract::extract(
                &extract.reads_tsv,
                &extract.fasta_files,
                &extract.output_file,
            );
        }
        SubCommand::ExtractBam(extract_bam) => {
            extract_bam::extract_bam(
                &extract_bam.regions_file,
                &extract_bam.bam_file,
                extract_bam.flank_size,
            );
        }
        SubCommand::Regions(regions) => {
            regions::regions(&regions.bam_files, !regions.no_chr_only);
        }
        SubCommand::SortBam(sort_bam) => {
            sort_bam::sort_bam(
                &sort_bam.input,
                &sort_bam.output,
                sort_bam.max_de,
                sort_bam.min_mapq,
            );
        }
    }
}
