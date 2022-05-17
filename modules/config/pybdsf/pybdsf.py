###############################################################################################
# 
#    * * *   P Y B D S F   C O N T A I N E R  P R O C E S S I N G   S C R I P T    * * *
#

# TO-DO in progress...
#    Debug cmds.
#       External:
#          $ docker run --rm -t -v /Users/susy/cirada/emu_pipeline/.escritoire/data/hydra:/home/pybdsf/data pybdsf /home/pybdsf/data /home/pybdsf/processing /home/pybdsf/data emu_simulated_04.shallow.sample.fits --thresh-pix 2.2794117647058822 --thresh-isl 2.0 --frequency 943490740.7407 --box-size 90 --step-size 45 --residual
#       Local:
#          # python3 pybdsf.py data processing results emu_simulated_04.shallow.sample.fits --thresh-pix 2.2794117647058822 --thresh-isl 2.0 --frequency 943490740.7407 --box-size 90 --step-size 45 --residual
#

import re
import os
import sys
import glob
import click
from shutil import copy2 as copy

import bdsf as bd
from astropy.io import fits
from astropy.table import Table

import numpy as np
from astropy import units as u
from astropy.wcs import WCS

# required to create empty catalogue
catalogue_header = {
    'Source_id': '',
    'Isl_id': '',
    'RA': 'deg',
    'E_RA': 'deg',
    'DEC': 'deg',
    'E_DEC': 'deg',
    'Total_flux': 'Jy',
    'E_Total_flux': 'Jy',
    'Peak_flux': 'Jy/beam',
    'E_Peak_flux': 'Jy/beam',
    'RA_max': 'deg',
    'E_RA_max': 'deg',
    'DEC_max': 'deg',
    'E_DEC_max': 'deg',
    'Maj': 'deg',
    'E_Maj': 'deg',
    'Min': 'deg',
    'E_Min': 'deg',
    'PA': 'deg',
    'E_PA': 'deg',
    'Maj_img_plane': 'deg',
    'E_Maj_img_plane': 'deg',
    'Min_img_plane': 'deg',
    'E_Min_img_plane': 'deg',
    'PA_img_plane': 'deg',
    'E_PA_img_plane': 'deg',
    'DC_Maj': 'deg',
    'E_DC_Maj': 'deg',
    'DC_Min': 'deg',
    'E_DC_Min': 'deg',
    'DC_PA': 'deg',
    'E_DC_PA': 'deg',
    'DC_Maj_img_plane': 'deg',
    'E_DC_Maj_img_plane': 'deg',
    'DC_Min_img_plane': 'deg',
    'E_DC_Min_img_plane': 'deg',
    'DC_PA_img_plane': 'deg',
    'E_DC_PA_img_plane': 'deg',
    'Isl_Total_flux': 'Jy',
    'E_Isl_Total_flux': 'Jy',
    'Isl_rms': 'Jy/beam',
    'Isl_mean': 'Jy/beam',
    'Resid_Isl_rms': 'Jy/beam',
    'Resid_Isl_mean': 'Jy/beam',
    'S_Code': ''
}


# pybdsf defaults config
# notes:
#    [1] https://www.astron.nl/citt/pybdsf/ug_basics.html#quick-start-example
pybdsf_cfg = {
    'atrous_do': True,
    'flagging_opts': True,
    'flag_maxsize_bm': 100,
    'rms_map': False,
    'mean_map': 'zero',
    'thresh_pix': 5.0,
    'thresh_isl': 3.0,
    'interactive': False,
    'quiet': True,
}


def cleenup_fits_file(fits_file):
    if os.path.isfile(fits_file):
        hdul = fits.open(fits_file)
        cols = hdul[1].columns
        for key,value in catalogue_header.items(): 
            if not value is None:
                cols.change_attrib(key,'unit',value)
        hdul.writeto(fits_file,overwrite=True)


def cleanup_csv_file_header(csv_file):
    if os.path.isfile(csv_file):
        with open(csv_file,"r+") as f:
            file_contents = list()
            is_keep = False
            for line in f:
                if re.match(r"^\# Source_id",line):
                    line =  re.sub(r"^# ","",line)
                    is_keep = True
                if is_keep:
                    file_contents.append(line)
            f.seek(0)
            f.truncate(0)
            for line in file_contents:
                f.write(line)

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

def process_image(fits_image_file,thresh_pix,thresh_isl,frequency=None,rms_box=None,quiet=pybdsf_cfg['quiet']):
    # ensure thresh_pix > thresh_pix; otherwise, PyBDSF hangs.
    if thresh_isl >= thresh_pix:
        print(f"WARNING: thresh_isl={thresh_isl} >= thresh_pix={thresh_pix} not recommended!")

    # ok, let's dance!
    kwargs = {}
    if not frequency is None:
        kwargs['frequency'] = frequency
    if not rms_box is None and len(rms_box)==2 and rms_box[0]>0 and rms_box[1]>0:
        kwargs['rms_box'] = rms_box # (step,size) in pixels
    try:
        img = bd.process_image(
            input=fits_image_file,
            atrous_do       = pybdsf_cfg['atrous_do'],
            flagging_opts   = pybdsf_cfg['flagging_opts'],
            flag_maxsize_bm = pybdsf_cfg['flag_maxsize_bm'],
            rms_map         = pybdsf_cfg['rms_map'],
            mean_map        = pybdsf_cfg['mean_map'],
            thresh_pix      = thresh_pix,
            thresh_isl      = thresh_isl,
            interactive     = pybdsf_cfg['interactive'],
            quiet           = quiet,
            **kwargs
        )
    except Exception as e:
        print(f"ERROR: {e}")
        return None

    return img

# debug: docker shell commands
# python3 pybdsf.py data/ processing/ results/ J165209+212528_s3arcmin_VLASS.fits --fits --thresh-isl 3.878470588235294 --thresh-pix 3.8823529411764706 --frequency 2.987741489322E+09 --residual
# python3 pybdsf.py data/ processing/ results/ emu_simulated_04.fits --fits --thresh-isl 3.7185 --thresh-pix 3.7222222222222223 --frequency 9.434907407407E+08 --residual
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
    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".pybdsf.residual.model.fits",f"{output_dir}/{image_filename}")
    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".pybdsf.residual.fits",f"{output_dir}/{image_filename}")
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
    #print(qt.colnames) # debug

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
        # skip bad fits
        if np.isnan(row['Maj']) or np.isnan(row['Min']) or np.isnan(row['PA']):
            continue

        # compute (ra,dec) in pixel units
        ra  = (row['RA']*eval(f"u.{catalogue_header['RA']}")).to(u.deg).value
        dec = (row['DEC']*eval(f"u.{catalogue_header['DEC']}")).to(u.deg).value
        r_pix = wcs.all_world2pix([[ra,dec]],0,ra_dec_order=True)[0]
        ra_0  = r_pix[0]
        dec_0 = r_pix[1]

        # compute semi-major/minor axes in pixel units
        pa = u.Quantity(row['PA'],catalogue_header['PA']).to(u.rad).value
        a  = u.Quantity(row['Maj'],catalogue_header['Maj']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        b  = u.Quantity(row['Min'],catalogue_header['Min']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        smaj = a*np.sqrt((np.sin(pa)/CDELT1)**2+(np.cos(pa)/CDELT2)**2)
        smin = b*np.sqrt((np.cos(pa)/CDELT1)**2+(np.sin(pa)/CDELT2)**2)

        # get peak flux
        pf = row['Peak_flux'] # NB: assumes units of Jy/beam.
        
        # create the gaussian function
        #gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,image[int(np.round(dec_0)),int(np.round(ra_0))]) # debug
        gaussian = make_gaussian(ra_0,dec_0,smaj,smin,pa,pf)

        # make residuals
        #print(f"ra_0={ra_0}, scale={scale}, smaj={smaj} (Maj: {row['Maj']}), smin={smin} (Min: {row['Min']}), pa={pa} (PA: {row['PA']})")
        ra_min  = max(0,int(np.floor(ra_0-scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
        ra_max  = min(n_x,int(np.ceil(ra_0+scale*np.sqrt((smaj*np.sin(pa))**2+(smin*np.cos(pa))**2))))
        dec_min = max(0,int(np.floor(dec_0-scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
        dec_max = min(n_y,int(np.ceil(dec_0+scale*np.sqrt((smaj*np.cos(pa))**2+(smin*np.sin(pa))**2))))
        i,j = np.mgrid[ra_min:ra_max+1,dec_min:dec_max+1]
        model[j,i] += gaussian(i,j)

    # create residual image
    image -= model

    # output model image
    model.shape = shape
    print(f">>    OUTPUT_MODEL: {model_file}")
    fits.PrimaryHDU(model,header=header).writeto(model_file,overwrite=True)

    # output residual image
    image.shape = shape
    print(f">>    OUTPUT_IMAGE: {residual_file}")
    fits.PrimaryHDU(image,header=header).writeto(residual_file,overwrite=True)
    print(f"> [Done]")

@click.command()
@click.argument('input_dir',nargs=1)
@click.argument('processing_dir',nargs=1)
@click.argument('output_dir',nargs=1)
@click.argument('fits_image_file',nargs=1)
@click.option(
    # Notes: Threshold for the island boundary in number of sigma above the mean. Determines extent of
    # island used for fitting
    '--thresh-pix',
    default      = pybdsf_cfg['thresh_pix'],
    type         = float,
    show_default = True, 
    help         = "Island threshold boundary (in \u03C3's)."
)
@click.option(
    # Notes: Source detection threshold: threshold for the island peak in number of sigma above the mean.
    # If false detection rate thresholding is used, this value is ignored and thresh_pix is calculated 
    # inside the program.
    '--thresh-isl',
    default      = pybdsf_cfg['thresh_isl'],
    type         = float,
    show_default = True,
    help         = "Source detection threshold (in \u03C3's)."
)
@click.option(
    '--frequency',
    default      = None,
    type         = float,
    show_default = True,
    help         = "Input frequency (optional: if in fits header)."
)
@click.option('--box-size', type=int,default=None,help="Grid RMS Box Size (requires: --step-size).")
@click.option('--step-size',type=int,default=None,help="Grid RMS Step Size (requires: --box-size).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and model FITS files.")
@click.option('--dump',is_flag=True,default=False,show_default=True,help="Dump out all processing files.")
def process(
    input_dir,
    processing_dir,
    output_dir,
    fits_image_file,
    thresh_pix,
    thresh_isl,
    box_size,
    step_size,
    frequency,
    fits,
    residual,
    dump
):
    """\b
       PyBDSF image processing tool.

       inputs:

          \b
          INPUT_DIR: location of image_filename.fits
          PROCESSING_DIR: location of scratch directory
          OUTPUT_DIR: location to place results
          FITS_IMAGE_FILE: image_filename.fits (without path)

       outputs:

          \b
          OUTPUT_DIR/image_filename.pybdsf.fits
          OUTPUT_DIR/image_filename.pybdsf.reg
    """
    if thresh_isl >= thresh_pix:
        print(f"WARNING: thresh_isl={thresh_isl} > thresh_pix={thresh_pix} not recommended!")

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
    print(f">    --thresh-pix: {thresh_pix}")
    print(f">    --thresh-isl: {thresh_isl}")

    # check if image file exists
    if not os.path.isfile(src):
        print(f"ERROR: Image file '{src}' not found!")
        print(f"Bye!")
        exit()

    # link image file from data to processing dir
    # TO-DO: Handle the Taylor term problem...
    dst = f"{processing_dir}/{fits_image_file}"
    #print(f"> Getting image file:")
    #print(f"> $ cp {src} {dst}")
    #copy(src,dst)
    print(f"> Linking image file:")
    print(f"> $ ln -s {dst} {src}")
    # this gaurd condition is usefull for debugging.
    if not (os.path.islink(dst) or os.path.exists(dst)):
        os.symlink(src,dst)

    # execute pybdsf
    print(f"> Running PyBDSF")
    rms_box = None
    if not box_size is None and not step_size is None:
        if step_size > box_size:
            print(f"ERROR: step_size > box_size")
            print(f">  box_size: {box_size}")
            print(f"> step_size: {step_size}")
            print(f"Bye!")
            exit()
        #rms_box = (step_size,box_size)
        rms_box = (box_size,step_size)
        print(f">> Setting: rms_box = {rms_box}")
    elif (not box_size is None and step_size) is None or (box_size is None and not step_size is None):
        print(f"ERROR: {'box_size' if box_size is None else 'step_size'} is undefined!")
        print(f">  box_size: {box_size}")
        print(f"> step_size: {step_size}")
        print(f"Bye!")
        exit()
    img = process_image(
        fits_image_file = dst,
        thresh_pix = thresh_pix,
        thresh_isl = thresh_isl,
        frequency = frequency,
        rms_box = rms_box
    )

    # create catalogue
    cat_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",f".pybdsf.{'fits' if fits else 'csv'}",fits_image_file)
    print(f"> Writing Catalogue: {cat_file}")
    cat_src = f"{processing_dir}/{cat_file}"
    cat_dst = f"{output_dir}/{cat_file}"
    if isinstance(img,bd.image.Image) and img.write_catalog(outfile=cat_src,format=('fits' if fits else 'csv'),catalog_type='srl',clobber=True):
        if not fits:
            cleanup_csv_file_header(cat_src)
        else:
            cleenup_fits_file(cat_src)
    else:
        print("> Whoops! No output!")
        print("> Creating empty catalogue.")
        csv_file = re.sub(r"\.fits$",".csv",cat_src) if fits else cat_src
        print(f" ==> {csv_file}")
        with open(csv_file,"w") as fd:
            fd.write(",".join(catalogue_header.keys())+"\n")
        if fits:
            Table.read(csv_file).write(cat_src,format='fits',overwrite=True)
            cleenup_fits_file(cat_src)

    # create region file
    reg_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$",".pybdsf.reg",fits_image_file)
    print(f"> Writing Region File: {reg_file}")
    reg_src = f"{processing_dir}/{reg_file}"
    reg_dst = f"{output_dir}/{reg_file}"
    if isinstance(img,bd.image.Image) and img.write_catalog(outfile=reg_src,format='ds9',catalog_type='srl',clobber=True):
        cleanup_reg_file(reg_src)
    else:
        print("> Whoops! No output!")
        print("> Creating empty region file.")
        with open(reg_src,"w") as fd:
            fd.write("# Empty -- no sources.\n")

    print(f"> pybdsf files: {processing_dir}/")
    for f in glob.glob(f"{processing_dir}/*"):
        f = re.sub(r'^(.*?/)*','',f)
        if f != fits_image_file:
            print(f"> o {f}")

    # copy result from process to data dir
    print(f"> Transfering results:")
    print(f"> $ cp {cat_src} {cat_dst}")
    copy(cat_src,cat_dst)
    print(f"> $ cp {reg_src} {reg_dst}")
    copy(reg_src,reg_dst)

    if residual:
        create_residual_image(fits_image_file,cat_file,processing_dir,output_dir)

    if dump:
        log_file = re.sub("\.([Ff][Ii][Tt]([Ss]|))$",r".\1.pybdsf.log",fits_image_file)
        log_src = f"{processing_dir}/{log_file}"
        log_dst = "{}/{}".format(output_dir,re.sub(r"\.[Ff][Ii][Tt]([Ss]|)\.(pybdsf.log)$",r".\2",log_file))
        print(f"> $ cp {log_src} {log_dst}")
        copy(log_src,log_dst)


if __name__ == "__main__":
    process()




