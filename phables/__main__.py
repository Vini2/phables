"""
Entrypoint for phables

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("phables.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("phables.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


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
            "--profile", help="Snakemake profile", default=None, show_default=False
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
        click.option(
            "--system_config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
            type=click.Path(),
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def run_options(func):
    """Command line args for run subcommand etc"""
    options = [
        click.option(
            "--input",
            help="Path to assembly graph file in .GFA format",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "--reads",
            help="Path to directory containing paired-end reads",
            type=click.Path(exists=True),
            required=True,
        ),
        click.option(
            "--minlength",
            default=2000,
            required=False,
            help="minimum length of circular unitigs to consider",
            type=int,
            show_default=True,
        ),
        click.option(
            "--mincov",
            default=10,
            required=False,
            help="minimum coverage of paths to output",
            type=int,
            show_default=True,
        ),
        click.option(
            "--compcount",
            default=200,
            required=False,
            help="maximum unitig count to consider a component",
            type=int,
            show_default=True,
        ),
        click.option(
            "--maxpaths",
            default=10,
            required=False,
            help="maximum number of paths to resolve for a component",
            type=int,
            show_default=True,
        ),
        click.option(
            "--mgfrac",
            default=0.2,
            required=False,
            help="length threshold to consider single copy marker genes",
            type=float,
            show_default=True,
        ),
        click.option(
            "--evalue",
            default=1e-10,
            required=False,
            help="maximum e-value for phrog annotations",
            type=float,
            show_default=True,
        ),
        click.option(
            "--seqidentity",
            default=0.3,
            required=False,
            help="minimum sequence identity for phrog annotations",
            type=float,
            show_default=True,
        ),
        click.option(
            "--covtol",
            default=100,
            required=False,
            help="coverage tolerance for extending subpaths",
            type=int,
            show_default=True,
        ),
        click.option(
            "--alpha",
            default=1.2,
            required=False,
            help="coverage multiplier for flow interval modelling",
            type=float,
            show_default=True,
        ),
        click.option(
            "--longreads",
            help="provide long reads as input (else defaults to short reads)",
            is_flag=True,
            default=False,
            show_default=True,
            required=False,
        ),
        click.option(
            "--prefix",
            help="prefix for genome identifier",
            type=str,
            required=False,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """
    Phables: from fragmented assemblies to high-quality bacteriophage genomes.
    Please refer the full documentation available on Read the Docs at https://phables.readthedocs.io/
    """
    pass


help_msg_extra = """
\b
\b
If you use Phables in your work, please cite Phables as,
\b
Vijini Mallawaarachchi, Michael J Roach, Przemyslaw Decewicz, 
Bhavya Papudeshi, Sarah K Giles, Susanna R Grigson, George Bouras, 
Ryan D Hesse, Laura K Inglis, Abbey L K Hutton, Elizabeth A Dinsdale, 
Robert A Edwards, Phables: from fragmented assemblies to high-quality 
bacteriophage genomes, Bioinformatics, Volume 39, Issue 10, 
October 2023, btad586, https://doi.org/10.1093/bioinformatics/btad586
\b
\b
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
Required:           phables run --input [assembly graph file]
Specify threads:    phables run ... --threads [threads]
Disable conda:      phables run ... --no-use-conda 
Change defaults:    phables run ... --snake-default="-k --nolock"
Add Snakemake args: phables run ... --dry-run --keep-going --touch
Specify targets:    phables run ... print_stages
Available targets:
    all             Run everything (default)
    preprocess      Run preprocessing only
    phables         Run phables (and preprocessing if needed)
    postprocess     Run postprocessing (with preprocessing and phables if needed)
    print_stages    List available stages
"""


# Run command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
@run_options
def run(**kwargs):
    """Run Phables"""

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "phables.smk")),
        merge_config=kwargs,
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
        snakefile_path=snake_base(os.path.join("workflow", "install.smk")),
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
def test(**kwargs):
    """Test Phables"""
    test_dir = snake_base("test_data")

    # Config to add or update in configfile
    merge_config = {"dir": test_dir}

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "test_phables.smk")),
        merge_config=merge_config,
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
cli.add_command(test)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
