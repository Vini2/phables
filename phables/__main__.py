"""
Entrypoint for phables

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from .util import (
    snake_base,
    get_version,
    default_to_output,
    copy_config,
    run_snakemake,
    OrderedCommands,
    print_citation,
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="phables.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="phables.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), '-v', '--version', is_flag=True)
def cli():
    """
        Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.
        Please refer the full documentation available on Read the Docs at https://phables.readthedocs.io/
    """
    pass


help_msg_extra = """
\b
\b
DOWNLOAD AND SETUP DATABASES
phables install
\b
\b
PREPROCESSING DATA
phables preprocess <options> 
\b
\b
RUN PHABLES
phables run <options> 
For more information on Phables please visit:
https://phables.readthedocs.io/
\b
\b
CLUSTER EXECUTION:
phables run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           phables run --input [folder]
Specify threads:    phables run ... --threads [threads]
Disable conda:      phables run ... --no-use-conda 
Change defaults:    phables run ... --snake-default="-k --nolock"
Add Snakemake args: phables run ... --dry-run --keep-going --touch
Specify targets:    phables run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


# Run command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--input",
    "_input",
    help="Path to Hecatomb output",
    type=click.Path(),
    required=True
)
@click.option(
    "--minlength",
    default=2000,
    required=False,
    help="minimum length of circular unitigs to consider",
    type=int,
)
@click.option(
    "--mincov",
    default=10,
    required=False,
    help="minimum coverage of paths to output",
    type=int,
)
@click.option(
    "--compcount",
    default=200,
    required=False,
    help="maximum unitig count to consider a component",
    type=int,
)
@click.option(
    "--maxpaths",
    default=10,
    required=False,
    help="maximum number of paths to resolve for a component",
    type=int,
)
@click.option(
    "--mgfrac",
    default=0.2,
    required=False,
    help="length threshold to consider single copy marker genes",
    type=float,
)
@click.option(
    "--alignscore",
    default=90,
    required=False,
    help="minimum alignment score for phrog annotations",
    type=float,
)
@click.option(
    "--seqidentity",
    default=0.3,
    required=False,
    help="minimum sequence identity for phrog annotations",
    type=float,
)
@common_options
def run(
    _input, 
    minlength, 
    mincov, 
    compcount, 
    maxpaths, 
    mgfrac, 
    alignscore, 
    seqidentity, 
    output, 
    log, 
    **kwargs):
    
    """Run Phables"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input, 
        "minlength": minlength,
        "mincov": mincov,
        "compcount": compcount,
        "maxpaths": maxpaths,
        "mgfrac": mgfrac,
        "alignscore": alignscore,
        "seqidentity": seqidentity,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "phables.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


# Install command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
def install(output, **kwargs):
    """Install databases"""

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join('workflow', 'install.smk')),
        **kwargs
    )


# Test command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
def test(output, **kwargs):
    """Test Phables"""
    test_dir = snake_base(os.path.join('..', 'tests', 'data'))

    # Config to add or update in configfile
    merge_config = {
        "dir": test_dir
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join('workflow', 'test_phables.smk')),
        merge_config=merge_config,
        **kwargs
    )


# Preprocess command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=click.Path(exists=True), required=True)
@click.option("--reads", help="Reads directory", type=click.Path(exists=True), required=True)
@common_options
def preprocess(_input, reads, output, log, **kwargs):
    """Preprocess data"""
    # Config to add or update in configfile
    merge_config = {"input": _input, "reads": reads, "output": output, "log": log}

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join('workflow', 'preprocess.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(install)
cli.add_command(preprocess)
cli.add_command(test)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
