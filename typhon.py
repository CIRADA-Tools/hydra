import re
import os
import sys
sys.path.insert(0,"modules")
from libs import Typhon
from libs import SpecialHelpOrder
from libs import HelpCommandPriority
import click

# define cache directory
def get_this_source_file_directory():
    def path(fullpath):
        return re.sub(r"(.*/).*$",r"\1",fullpath)
    return path(os.path.realpath(__file__))
escritoire = os.path.abspath(get_this_source_file_directory()+"/.escritoire")
cache = os.path.abspath(f"{escritoire}/data/typhon")

cmd=HelpCommandPriority()
@click.group(cls=SpecialHelpOrder)
def cli():
    """\b
       Generates source-finder catalogues by optimizing input parameters.

       Typhon: a powerful controlling force, farther of Cerberus.
    """
    pass

@cli.command(help_priority=cmd.prio())
@click.argument('fits_file',nargs=1)
@click.option('--diagnostics',is_flag=True,default=False,help="Include diagnostics run.")
@click.option('--optimize-rms-box',is_flag=True,default=False,help="Optimize RMS/Step Parameters.")
def optimize(
    fits_file,
    diagnostics,
    optimize_rms_box
):
    """\b
       Optimize multiple source-finder parameters.
    """
    Typhon(
        fits_image_file     = fits_file,
        processing_cache    = cache,
        is_diagnostics      = diagnostics,
        is_optimize_rms_box = optimize_rms_box
    ).optimize()

@cli.command(help_priority=cmd.prio())
@click.argument('fits_file',nargs=1)
@click.option('--use',default=None,type=str,help="Use FITS_FILE_PREFIX.typhon.tar.gz optimization.")
@click.option('--diagnostics',is_flag=True,default=False,help="Include diagnostics run.")
@click.option('--residual',is_flag=True,default=False,help="Output residual and model FITS files.")
def process(
    fits_file,
    use,
    diagnostics,
    residual
):
    """\b
       Processes FITS_FILE using optimized input parameters.

       It will generate a FITS_FILE_PREFIX.typhon.tar.gz if it does not exist,
       or use an alternative one, via, the --use flag.
    """
    Typhon(
        fits_image_file  = fits_file,
        processing_cache = cache,
        use              = use,
        is_diagnostics   = diagnostics,
        is_residual      = residual
    ).process()


if __name__ == "__main__":
    cli()
