"""
Entrypoint for MAGBuild

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click
from shutil import copyfile
from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    try:
        from MAGBuild._version import version
    except Exception:
        version = "0.1.0"
    return version


def print_citation():
    with open(snake_base("../MAGBuild.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value

# This is to ensure that the output directory gets generated if it hasnt already been generated
def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option('--input', '_input', help='Directory of reads', type=click.Path(), required=False, default='test/illumina-subset', show_default=True),
        click.option('--output', 'output', help='Output directory', type=click.Path(),
                     default='sphae.out', show_default=True),
        click.option("--configfile", default="config.yaml", show_default=False, callback=default_to_output,
                     help="Custom config file [default: config.yaml]"),
        click.option('--threads', help='Number of threads to use', default=1, show_default=True),
        click.option('--profile', help='Snakemake profile', default=None, show_default=False),
        click.option('--db_dir', 'db_dir', help='Custom database directory', type=click.Path(), required=False),
        click.option('--temp-dir', help='Temp directory', required=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'],
                     help="Customise Snakemake runtime args", show_default=True),
        click.option("--log", default="sphae.log", callback=default_to_output, hidden=True,),
        click.option("--system-config", default=snake_base(os.path.join("config", "config.yaml")),hidden=True,),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Assembling pure culture phages from both Illumina and Nanopore sequencing technology
    \b
    For more options, run:
    MAGBuild command --help"""
    pass


help_msg_run = """
\b
RUN EXAMPLES 
I dont know, when I do, I will update this :)
"""


@click.option('--output', 'output', help='Output directory', type=click.Path(), default='sphae.out', show_default=True)
@click.option("--configfile", default="config.yaml", show_default=False, callback=default_to_output,help="Custom config file [default: (outputDir)/config.yaml]",)
@click.option('--threads', help='Number of threads to use', default=1, show_default=True)
@click.option('--profile', help='Snakemake profile', default=None, show_default=False)
@click.option('--temp-dir', 'temp_dir', help='Temp directory', required=False)
@click.option('--snake-default', multiple=True,default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'], help="Customise Snakemake runtime args", show_default=True)
@click.option("--log", default="sphae.log", callback=default_to_output, hidden=True,)
@click.option("--system-config", default=snake_base(os.path.join("config", "config.yaml")),hidden=True,)
@click.argument("snake_args", nargs=-1)


@click.command(epilog=help_msg_run, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
@click.option('--sequencing', 'sequencing', help="sequencing method", default='paired', show_default=True, type=click.Choice(['paired', 'longread']))
def run(_input, output, sequencing, no_medaka, temp_dir, configfile, **kwargs):
    """Run MAGBuild"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    merge_config = {
        "args": {
            "input": _input, 
            "output": output, 
            "sequencing": sequencing,
            "configfile": configfile,
            "temp_dir": temp_dir,
        }
    }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'Snakefile')),
        configfile=configfile,
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@click.option('--configfile', default='config.yaml', help='Copy template config to file', show_default=True)
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(config)
cli.add_command(citation)

def main():
    cli()

if __name__ == '__main__':
    main()
