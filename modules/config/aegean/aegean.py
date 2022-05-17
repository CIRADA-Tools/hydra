###############################################################################################
# 
#    * * *   A E G E A N   C O N T A I N E R  P R O C E S S I N G   S C R I P T    * * *
#

import re
import os
import sys
import glob
import click
from shutil import copy2 as copy
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

# note: will need to remove the following lines from .dcr file to use these libs
#   && apt-get autoremove --purge -y python3-pip \
#from astropy.coordinates import SkyCoord

#from astropy.convolution import Gaussian2DKernel
#from astropy.modeling import Gaussian2D

from os.path import abspath

catalogue_header = {
    'island': None,
    'source': None,
    'background': 'Jy/beam',
    'local_rms': 'Jy/beam',
    'ra_str': 'HMS',
    'dec_str': 'DMS',
    'ra': 'deg',
    'err_ra': 'deg',
    'dec': 'deg',
    'err_dec': 'deg',
    'peak_flux': 'Jy/beam',
    'err_peak_flux': 'Jy/beam',
    'int_flux': 'Jy',
    'err_int_flux': 'Jy',
    'a': 'arcsec',
    'err_a': 'arcsec',
    'b': 'arcsec',
    'err_b': 'arcsec',
    'pa': 'deg',
    'err_pa': 'deg',
    'flags': None,
    'residual_mean': None,
    'residual_std': None,
    'uuid': None,
    'psf_a': 'arcsec',
    'psf_b': 'arcsec',
    'psf_pa': 'arcsec'
}

# aegean defaults config
aegean_cfg = {
    'seedclip':  5.0,
    'floodclip': 4.0,
}

def cleenup_fits_file(fits_file):
    if os.path.isfile(fits_file):
        hdul = fits.open(fits_file)
        cols = hdul[1].columns
        for key,value in catalogue_header.items(): 
            if not value is None:
                cols.change_attrib(key,'unit',value)
        hdul.writeto(fits_file,overwrite=True)

def cleanup_reg_file(reg_file):
    if os.path.isfile(reg_file):
        with open(reg_file,"r+") as f:
            file_contents = list()
            for line in f:
                if re.match(r"ellipse",line):
                    # remove the region labels
                    line = re.sub(r"\s*#.*","",line)
                file_contents.append(line)
            f.seek(0)
            f.truncate(0)
            for line in file_contents:
                f.write(line)

# debug: docker shell commands
# python3 aegean.py data/ processing/ results/ J165209+212528_s3arcmin_VLASS.fits --fits --floodclip 4.671794117647059 --seedclip 4.676470588235294 --residual
# python3 aegean.py data/ processing/ results/ emu_simulated_04.fits --fits --floodclip 3.996  --seedclip 4.0 --residual
# python3 aegean.py data processing results emu_simulated_04.pybdsf.residual.sample.inverted.fits --seedclip 4.0 --floodclip 2.0 --residual
# python3 aegean.py data/ processing/ results/ J085542+112459_s3arcmin_VLASS.shallow.fits --seedclip 4.0 --floodclip 2.0 --residual
def create_residual_image(image_filename,catalogue_filename,input_dir,output_dir):
    def get_wcs(header):
        header = hdul[0].header

        # minor header fix
        rotation_matrix_map = {
            'PC01_01': 'PC1_1',
            'PC02_01': 'PC2_1',
            'PC03_01': 'PC3_1',
            'PC04_01': 'PC4_1',
            'PC01_02': 'PC1_2',
            'PC02_02': 'PC2_2',
            'PC03_02': 'PC3_2',
            'PC04_02': 'PC4_2',
            'PC01_03': 'PC1_3',
            'PC02_03': 'PC2_3',
            'PC03_03': 'PC3_3',
            'PC04_03': 'PC4_3',
            'PC01_04': 'PC1_4',
            'PC02_04': 'PC2_4',
            'PC03_04': 'PC3_4',
            'PC04_04': 'PC4_4',
        }
        for key in rotation_matrix_map.keys():
            if key in header:
                header.insert(key,(rotation_matrix_map[key],header[key]),after=True)
                header.remove(key)

        # wcs image header
        wcs = WCS(header)

        # trim to 2d from nd
        naxis = wcs.naxis
        while naxis > 2:
            wcs = wcs.dropaxis(2)
            naxis -= 1

        return wcs

    # file definitions
    image_file     = f"{input_dir}/{image_filename}"
    catalogue_file = f"{input_dir}/{catalogue_filename}"
    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".aegean.residual.model.fits",f"{output_dir}/{image_filename}")
    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".aegean.residual.fits",f"{output_dir}/{image_filename}")
    print(f"> Creating Residual Image:")
    print(f">>    INPUT_IMAGE: {image_file}")
    print(f">>    INPUT_CATALOGUE: {catalogue_file}")

    # open image and catalogue files
    hdul = fits.open(image_file)
    header = hdul[0].header
    wcs = get_wcs(header)
    shape = hdul[0].data.shape
    image = np.squeeze(hdul[0].data)
    CDELT1 = np.abs(wcs.wcs.cdelt[0])
    CDELT2 = np.abs(wcs.wcs.cdelt[1])
    qt = Table.read(catalogue_file)
    #print(qt) # debug

    ## notes: https://astronomy.stackexchange.com/questions/26721/converting-jy-beam-to-jy
    #BMAJ = header['BMAJ']
    #BMIN = header['BMIN']
    #BPA  = np.deg2rad(header['BPA'])
    #BMAJ = BMAJ*np.sqrt((CDELT1*np.sin(BPA))**2+(CDELT2*np.cos(BPA))**2)/(CDELT1*CDELT2)
    #BMIN = BMIN*np.sqrt((CDELT1*np.cos(BPA))**2+(CDELT2*np.sin(BPA))**2)/(CDELT1*CDELT2)
    #beam_size = np.pi*BMAJ*BMIN/(4.0*np.log(2.0))

    # build our model image
    scale = 8
    n_x = image.shape[1]-1
    n_y = image.shape[0]-1
    model = np.zeros(image.shape)
    def make_gaussian(ra_0,dec_0,smaj,smin,pa,peak):
        a = ((np.cos(pa)/smin)**2+(np.sin(pa)/smaj)**2)/2.0
        b = (np.sin(2.0*pa)/(smin)**2-np.sin(2.0*pa)/(smaj)**2)/2.0
        c = ((np.sin(pa)/smin)**2+(np.cos(pa)/smaj)**2)/2.0
        def gaussian(ra,dec):
            return peak*np.exp(-a*(ra-ra_0)**2-b*(ra-ra_0)*(dec-dec_0)-c*(dec-dec_0)**2)
        return gaussian
    for row in qt:
        # compute (ra,dec) in pixel units
        ra  = (row['ra']*eval(f"u.{catalogue_header['ra']}")).to(u.deg).value
        dec = (row['dec']*eval(f"u.{catalogue_header['dec']}")).to(u.deg).value
        r_pix = wcs.all_world2pix([[ra,dec]],0,ra_dec_order=True)[0]
        ra_0  = r_pix[0]
        dec_0 = r_pix[1]

        # compute semi-major/minor axes in pixel units
        pa = (row['pa']*eval(f"u.{catalogue_header['pa']}")).to(u.rad).value
        a  = u.Quantity(row['a'],catalogue_header['a']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        b  = u.Quantity(row['b'],catalogue_header['b']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        smaj = a*np.sqrt((np.sin(pa)/CDELT1)**2+(np.cos(pa)/CDELT2)**2)
        smin = b*np.sqrt((np.cos(pa)/CDELT1)**2+(np.sin(pa)/CDELT2)**2)

        # get peak flux
        #pf = row['peak_flux']/beam_size
        pf = row['peak_flux'] # NB: assumes units of Jy/beam.
        
        # create the gaussian function
        #gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,image[int(np.round(dec_0)),int(np.round(ra_0))]) # debug
        gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,pf)

        # make residuals
        ra_min  = max(0,int(np.floor(ra_0-scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
        ra_max  = min(n_x,int(np.ceil(ra_0+scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
        dec_min = max(0,int(np.floor(dec_0-scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
        dec_max = min(n_y,int(np.ceil(dec_0+scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
        i,j = np.mgrid[ra_min:ra_max+1,dec_min:dec_max+1]
        model[j,i] += gaussian(i,j)

    # create residual image
    image -= model

    # output model image
    print(f">>    OUTPUT_MODEL: {model_file}")
    model.shape = shape
    fits.PrimaryHDU(model,header=header).writeto(model_file,overwrite=True)

    # output residual image
    print(f">>    OUTPUT_IMAGE: {residual_file}")
    image.shape = shape
    fits.PrimaryHDU(image,header=header).writeto(residual_file,overwrite=True)
    print(f"> [Done]")

#def create_residual_image(image_filename,catalogue_filename,input_dir,output_dir):
#    # file definitions
#    image_file     = f"{input_dir}/{image_filename}"
#    catalogue_file = f"{input_dir}/{catalogue_filename}"
#    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".aegean.residual.model.fits",f"{output_dir}/{image_filename}")
#    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".aegean.residual.fits",f"{output_dir}/{image_filename}")
#    print(f"> Creating Residual Image:")
#    print(f">>    INPUT_IMAGE: {image_file}")
#    print(f">>    INPUT_CATALOGUE: {catalogue_file}")
#    print(f">>    OUTPUT_MODEL: {model_file}")
#    print(f">>    OUTPUT_IMAGE: {residual_file}")
#    aeres = "AeRes -c {} -f {} -r {} -m {}".format(
#        catalogue_file,
#        image_file,
#        residual_file,
#        model_file
#    )
#    print(f">> Executing AeRes:")
#    print(f"> $ {aeres}")
#    sys.stdout.flush()
#    os.system(aeres)
#    print(f"> [Done]")

@click.command()
@click.argument('input_dir',nargs=1)
@click.argument('processing_dir',nargs=1)
@click.argument('output_dir',nargs=1)
@click.argument('fits_image_file',nargs=1)
@click.option(
    '--seedclip',
    default      = aegean_cfg['seedclip'],
    type         = float,
    show_default = True,
    help         = "The clipping value (in \u03C3's) for seeding islands."
)
@click.option(
    '--floodclip',
    default      = aegean_cfg['floodclip'],
    type         = float,
    show_default = True,
    help         = "The clipping value (in \u03C3's) for growing islands."
)
@click.option('--box-size', type=int,default=-1,help="Grid RMS Box Size. [default: ~4*beam_size]")
@click.option('--step-size',type=int,default=-1,help="Grid RMS Step Size. [default: 5*box_size]")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and model FITS files.")
@click.option('--dump',is_flag=True,default=False,help="Dump out all processing files.")
def process(
    input_dir,
    processing_dir,
    output_dir,
    fits_image_file,
    seedclip,
    floodclip,
    box_size,
    step_size,
    fits,
    residual,
    dump
):
    """\b
       Aegean image processing tool.

       inputs:

          \b
          INPUT_DIR: location of image_filename.fits
          PROCESSING_DIR: location of scratch directory
          OUTPUT_DIR: location to place results
          FITS_IMAGE_FILE: image_filename.fits (without path)

       outputs:

          \b
          OUTPUT_DIR/image_filename.aegean.csv
          OUTPUT_DIR/image_filename.aegean.reg
    """

    # house cleaning
    input_dir      = re.sub(r"/+$","",input_dir)
    processing_dir = re.sub(r"/+$","",processing_dir)
    output_dir     = re.sub(r"/+$","",output_dir)

    src = f"{input_dir}/{fits_image_file}"
    print(f"Processing: {src}")
    print(f"> Local Context: ")
    print(f">   INPUT_DIR:  {input_dir}")
    print(f">   PROCESSING_DIR: {processing_dir}")
    print(f">   OUTPUT_DIR: {output_dir}")
    print(f"> Supported Flags: ")
    print(f">    --seedclip: {seedclip}")
    print(f">    --floodclip: {floodclip}")

    # check if image file exists
    if not os.path.isfile(src):
        print(f"ERROR: Image file '{src}' not found!")
        print(f"Bye!")
        exit()

    # link image file from data to processing dir
    dst = f"{processing_dir}/{fits_image_file}"
    #print(f"> Getting image file:")
    #print(f"> $ cp {src} {dst}")
    #copy(src,dst)
    print(f"> Linking image file:")
    print(f"> $ ln -s {dst} {src}")
    # this gaurd condition is usefull for debugging.
    if not (os.path.islink(dst) or os.path.exists(dst)):
        os.symlink(src,dst)

    # execute bane
    bane = f"BANE {dst} --cores 1"
    print(f"> Executing BANE:")
    print(f"> $ {bane}")
    if step_size > 0 and box_size > 0 and step_size > box_size:
        print(f"ERROR: step_size > box_size")
        print(f">  box_size: {box_size}")
        print(f"> step_size: {step_size}")
        print(f"Bye!")
        exit()
    if step_size > 0:
        bane += f" --grid {step_size} {step_size}" # 16 
    if box_size > 0:
        bane += f" --box {box_size} {box_size}" # 96
    sys.stdout.flush()
    exit_code = os.system(bane)
    if exit_code != 0:
        print("> Whoops! BANE error!")
        print("Bye!")
        exit(exit_code)
    print(f"> BANE files: {processing_dir}/")
    files = [re.sub(r'^(.*?/)*','',f) for f in glob.glob(f"{processing_dir}/*")]
    for f in files:
        if f != fits_image_file:
            print(f">   o {f}")

    # execute aegean
    cat_file  = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",f".aegean.{'fits' if fits else 'csv'}",fits_image_file)
    reg_file  = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",".aegean.reg",fits_image_file)
    aegean = "aegean --autoload {} --table {} --cores 1 --seedclip {} --floodclip {}".format(
        dst,
        f"{processing_dir}/{cat_file},{processing_dir}/{reg_file}",
        seedclip,
        floodclip
    )
    print(f"> Executing aegean:")
    print(f"> $ {aegean}")
    sys.stdout.flush()
    exit_code = os.system(aegean)
    if exit_code != 0:
        print("> Whoops! Aegean error!")
        print("Bye!")
        exit(exit_code)
    print(f"> aegean files: {processing_dir}/")
    for f in glob.glob(f"{processing_dir}/*"):
        f = re.sub(r'^(.*?/)*','',f)
        if not f in files:
            print(f">   o {f}")

    # copy result from process to data dir
    print(f"> Transfering results:")
    cat_comp_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",f".aegean_comp.{'fits' if fits else 'csv'}",fits_image_file)
    reg_comp_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",".aegean_comp.reg",fits_image_file)
    src = f"{processing_dir}/{cat_comp_file}"
    if fits:
        cleenup_fits_file(src) # add units
    dst = f"{output_dir}/{cat_file}"
    print(f"> $ cp {src} {dst}")
    #try:
    #    copy(src,dst)
    #except Exception as FileNotFoundError:
    #    print("> Whoops! No output!")
    #    print("> Creating empty catalogue.")
    #    csv_file = re.sub(r"\.fits$",".csv",src) if fits else dst
    #    with open(csv_file,"w") as fd:
    #        fd.write(",".join(catalogue_header.keys())+"\n")
    #    if fits:
    #        Table.read(csv_file).write(dst,format='fits',overwrite=True)
    #        cleenup_fits_file(dst) # add units
    if not os.path.isfile(src):
        print("> Whoops! No output!")
        print("> Creating empty catalogue.")
        csv_file = re.sub(r"\.fits$",".csv",src)
        with open(csv_file,"w") as fd:
            fd.write(",".join(catalogue_header.keys())+"\n")
        if fits:
            Table.read(csv_file).write(src,format='fits',overwrite=True)
            cleenup_fits_file(src) # add units
    copy(src,dst)
    #dst = f"{output_dir}/{cat_file}"
    src = f"{processing_dir}/{reg_comp_file}"
    cleanup_reg_file(src) # remove annotations
    dst = f"{output_dir}/{reg_file}"
    print(f"> $ cp {src} {dst}")
    #try:
    #    copy(src,dst)
    #except Exception as FileNotFoundError:
    #    print("> Whoops! No output!")
    #    print("> Creating empty region file.")
    #    with open(dst,"w") as fd:
    #        fd.write("# Empty -- no sources.\n")
    if not os.path.isfile(src):
        print("> Whoops! No output!")
        print("> Creating empty region file.")
        with open(src,"w") as fd:
            fd.write("# Empty -- no sources.\n")
    copy(src,dst)

    # create residual file...
    if residual:
        create_residual_image(fits_image_file,cat_comp_file,processing_dir,output_dir)

    if dump:
        rms_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$","_rms.fits",fits_image_file)
        bkg_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$","_bkg.fits",fits_image_file)
        src = f"{processing_dir}/{rms_file}"
        dst = "{}/{}".format(output_dir,re.sub(r"_(rms\.fits)$",r".aegean.\1",rms_file))
        print(f"> $ cp {src} {dst}")
        copy(src,dst)
        src = f"{processing_dir}/{bkg_file}"
        dst = "{}/{}".format(output_dir,re.sub(r"_(bkg\.fits)$",r".aegean.\1",bkg_file))
        print(f"> $ cp {src} {dst}")
        copy(src,dst)

    # done
    print(f"[Done]")

if __name__ == "__main__":
    process()
