import re
import os
import sys
sys.path.insert(0,"../")
from libs.config_tools import cfg_fetch_beam_pars
from libs.fits_tools  import FITS
from libs.fits_tools  import Homados
from libs.fits_tools  import get_fits_image_data
from libs.fits_tools  import sigma_clipper
from libs.click_tools import SpecialHelpOrder
from libs.click_tools import HelpCommandPriority
from libs.exceptions  import print_warning
from libs.exceptions  import print_ok
from libs.fits_tools import invert_fits_image

import numpy as np
import click

# fits io handlling
from astropy.io import fits
from warnings import simplefilter as warning_filter
from astropy.io.fits.verify import VerifyWarning
warning_filter('ignore',category=VerifyWarning)


############################################################
#
#    * * * C O M A N D - L I N E   I N T E R F A C E * * *
#

# initialize our help command
cmd=HelpCommandPriority()
@click.group(cls=SpecialHelpOrder)
def homados():
    """\b
       Adds Gaussian noise to FITS images.

       Homados: the personification of battle noise.
    """
    pass

@homados.command(help_priority=cmd.prio())
@click.argument('fits_deep_image',nargs=1)
@click.option('--output-dir', default=None, type=str, help="Results output directory.")
@click.option('--n-sigma', default=5, type=float, show_default=True, help="Target noise threshold (n\u03C3).")
@click.option('--BMAJ', default=None, type=float, help="Major beam axis [\u00B0].")
@click.option('--BMIN', default=None, type=float, help="Minor beam axis [\u00B0].")
@click.option('--BPA', default=None, type=float, help="Beam position angle [\u00B0].")
@click.option('--no-noise', is_flag=True, default=False, show_default=True, help="Don't Keep noise file.")
@click.option('--verbose', is_flag=True, default=False, show_default=True, help="Verbose output.")
def shallow(
    fits_deep_image,
    output_dir,
    n_sigma,
    bmaj,
    bmin,
    bpa,
    no_noise,
    verbose
):
    """\b
       Creates shallow image from FITS_FILE deep/reference image.
    """
    if not os.path.isfile(fits_deep_image):
        print_warning(f"ERROR: File '{fits_deep_image}' does not exist!")
        print_ok("Bye!")
        return

    # we require beam pars -- check file, command-line, and configuraiton; otherwise, ask user
    fd = FITS(fits_deep_image)
    is_beam_pars = (bmaj != None and bmin != None and bpa != None) or fd.has_beam_shape()
    if not is_beam_pars: # ok, not in file or on command-line -- check configuration
        survey = fd.get_survey()
        beam_pars = cfg_fetch_beam_pars(survey)
        if not beam_pars is None: # it's in the configuration -- grab it
            bmaj = beam_pars['BMAJ']
            bmin = beam_pars['BMIN']
            bpa  = beam_pars['BPA']
            msg = ""
            msg += f"WARNING: Beam shape missing from header: {fits_deep_image}\n"
            msg += f"> Using {survey} configuration defaults:\n"
            msg += f">   BMAJ: {bmaj}\n"
            msg += f">   BMIN: {bmin}\n"
            msg += f">   BPA:  {bpa}"
            print_warning(msg)
            is_beam_pars = True
    del fd # free some memory

    # ok, let's go!
    if is_beam_pars: # process...
        Homados(fits_deep_image,output_dir,verbose).make_noise(n_sigma,{'BMAJ': bmaj, 'BMIN': bmin, 'BPA': bpa},not no_noise)
    else: # whoops!
        msg = ""
        msg += f"Beam shape missing from header: {fits_deep_image}\n"
        msg += f"> Require Inputs:\n"
        msg += f">   --BMAJ\n"
        msg += f">   --BMIN\n"
        msg += f">   --BPA"
        print(msg)

@homados.command(help_priority=cmd.prio())
@click.argument('fits_file',nargs=1)
@click.option('--output-dir',default=None,type=str,help="Results output directory.")
def invert(
    fits_file,
    output_dir
):
    """\b
       Inverts FITS_FILE image.
    """
    inverted_fits_file = invert_fits_image(fits_file,output_dir)
    print(f"Computing stats...")
    stats = sigma_clipper(inverted_fits_file)
    rms  = stats['rms']*10**6
    mean = stats['mean']*10**6
    median = stats['median']*10**6
    print(f"Median: {median} \u03bcJy")
    print(f"Mean:   {mean} \u03bcJy")
    print(f"RMS:    {rms} \u03bcJy")
    print(f"[Done]")


@homados.command(help_priority=cmd.prio())
@click.argument('fits_file',nargs=1)
def statistics(
    fits_file
):
    """\b
       Computes RMS, mean, etc., of FITS_FILE image.
    """
    if os.path.isfile(fits_file):
        print(f"> Computing stats: {fits_file}")
        try:
            stats = sigma_clipper(fits_file)
            rms  = stats['rms']*10**6
            mean = stats['mean']*10**6
            median = stats['median']*10**6
            initial_min = stats['initial_min']*10**6
            initial_max = stats['initial_max']*10**6
            final_min = stats['final_min']*10**6
            final_max = stats['final_max']*10**6
        except:
            print_warning("Whoops! This can't be a fits file.")
            print_ok("Bye!")
            return
        print(f"Median: {median} \u03bcJy")
        print(f"Mean:   {mean} \u03bcJy")
        print(f"RMS:    {rms} \u03bcJy")
        print(f"Initial Range: [{initial_min}, {initial_max}] \u03bcJy")
        print(f"Final Range:   [{final_min}, {final_max}] \u03bcJy")
    else:
        print_warning(f"ERROR: File '{fits_file}' does not exist!")
        print_ok("Bye!")
        return
    print(f"[Done]")



if __name__ == "__main__":
    homados()
