###############################################################################################
# 
#    * * *   S E L A V Y   C O N T A I N E R  P R O C E S S I N G   S C R I P T    * * *
#

import re
import os
import sys
import glob
import click
from shutil import copy2 as copy

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS


# required to create empty catalogue
catalogue_header = {
    'island_id': '',
    'component_id': '',
    'component_name': '',
    'ra_hms_cont': 'HMS',
    'dec_dms_cont': 'DMS ',
    'ra_deg_cont': 'deg',
    'dec_deg_cont': 'deg',
    'ra_err': 'arcsec',
    'dec_err': 'arcsec',
    'freq': 'MHz ',
    'flux_peak': 'mJy/beam',
    'flux_peak_err': 'mJy/beam',
    'flux_int': 'mJy',
    'flux_int_err': 'mJy',
    'maj_axis': 'arcsec',
    'min_axis': 'arcsec',
    'pos_ang': 'deg',
    'maj_axis_err': 'arcsec',
    'min_axis_err': 'arcsec',
    'pos_ang_err': 'deg',
    'maj_axis_deconv': 'arcsec',
    'min_axis_deconv': 'arcsec',
    'pos_ang_deconv': 'deg',
    'maj_axis_deconv_err': 'arcsec',
    'min_axis_deconv_err': 'arcsec',
    'pos_ang_deconv_err': 'deg',
    'chi_squared_fit': '',
    'rms_fit_gauss': 'mJy/beam',
    'spectral_index': '',
    'spectral_curvature': '',
    'spectral_index_err': '',
    'spectral_curvature_err': '',
    'rms_image': 'mJy/beam',
    'has_siblings': '',
    'fit_is_estimate': '',
    'spectral_index_from_TT': '',
    'flag_c4': '',
    'comment': ''
}


# selavy default template for creating selavy.in file.
# notes: 
#    [1] https://www.atnf.csiro.au/computing/software/askapsoft/sdp/docs/current/analysis/selavy.html
#selavy_pars = {
#    'Selavy': {
#        'image': None,
#        'imagetype': 'fits',
#        'flagLog': True,
#        'flagDS9': True,
#        'Fitter': {
#            'doFit': True,
#            'fitTypes': '[full]',
#            'numGaussFromGuess': True,
#        },
#        'snrCut': 4.0,
#        'flagGrowth': True,
#        'growthCut': 3.0,
#     },
#}
selavy_pars = {
    'Selavy': {
        'image': None,
        'imagetype': 'fits',
        'flagLog': True,
        'flagDS9': True,
        'Fitter': {
            'doFit': True,
            'fitTypes': '[full]',
            'numGaussFromGuess': True,
            'writeComponentMap': True,
        },
        'searchType': 'spatial',
        'VariableThreshold': True,
        'flagRobustStats': True,
        #'flagRobustStats': False,
        'snrCut': 4.0,
        'flagGrowth': True,
        'growthCut': 3.0,
     },
}

# sets selavy_pars
def set_pars(image,snrCut,growthCut,box_size=None,step_size=None):
    selavy_pars['Selavy']['image'] = f'"{image}"'
    if not snrCut is None:
        selavy_pars['Selavy']['snrCut'] = snrCut
    if not growthCut is None:
        selavy_pars['Selavy']['growthCut'] = growthCut
    if not box_size is None: # nb: apperently selavy has no step_size option
        if  'VariableThreshold' in selavy_pars['Selavy']:
           selavy_pars['Selavy']['__VariableThreshold'] = selavy_pars['Selavy']['VariableThreshold']
        selavy_pars['Selavy']['VariableThreshold'] = {'boxSize': int((float(box_size)-1.0)/2.0)}
    return selavy_pars

# sets and flattens selavy_pars into string for writing to selavy.in config file
def get_config(image,snrCut,growthCut,box_size=None,step_size=None):
    def flatten(template,r_key=""):
        prefix = f"{'' if r_key == '' else r_key+'.'}"
        flattened = ""
        for key in template:
            value = template[key]
            if isinstance(value,dict):
                flattened += "\n".join([f"{prefix if e != '' else ''}{e}" for e in flatten(value,key).split('\n')])
            else:
                if isinstance(value,bool):
                    value = f"{value}".lower()
                flattened += f"{prefix}{key} = {value}\n"
        return flattened if r_key != '' else re.sub("\n$","",flattened)
    return re.sub("__","",flatten(set_pars(image,snrCut,growthCut,box_size,step_size)))


# txt_file to csv_file copy-converter
def selavy_txt_to_csv(txt_file):
    csv_file = re.sub(r"\.txt$",".csv",txt_file)
    if os.path.isfile(txt_file):
        # NB: selavy bug work around: i.e., column seperators in componets/islands .txt files removed if len(table)==0 -- fix only for component files
        is_kludge = False if re.sub(r"^(.*?/)+","",txt_file) == "selavy-results.components.txt" and sum(1 for line in open(txt_file)) > 2 else True
        file_contents = list()
        if not is_kludge:
            with open(txt_file,"r") as f:
                is_first = True
                for line in f:
                    line = re.sub(r"^\s+","",line)
                    if re.match(r"^\s*?#",line):
                        if is_first and re.search(r"(island_id|ObjID)",line):
                            line = re.sub(r"^#\s*","",line)
                            is_first = False
                        else:
                            continue
                    line = re.sub(r" +",",",line)
                    file_contents.append(line)
        else:
            file_contents.append(",".join(catalogue_header))
        with open(csv_file,"w") as f:
            for line in file_contents:
                f.write(line)
    return csv_file


def cleenup_fits_file(fits_file):
    # NB: This routine only handles fits component catalogues
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
                if re.match(r"text",line):
                    continue
                file_contents.append(line)
            f.seek(0)
            f.truncate(0)
            for line in file_contents:
                f.write(line)

# debug: docker shell commands
# python3 selavy.py data/ processing/ results/ J165209+212528_s3arcmin_VLASS.fits --fits --snrCut 3.8823529411764706 --growthCut 3.878470588235294 --residual
# python3 selavy.py data/ processing/ results/ emu_simulated_04.fits --fits --snrCut 4.0 --growthCut 3.996 --residual
# python3 selavy.py data processing results J165209+212528_s3arcmin_VLASS.profound.residual.fits --snrCut 4.0 --growthCut 2.0 --residual --fits
# python3 selavy.py data processing results J165209+212528_s3arcmin_VLASS.aegean.residual.shallow.fits --snrCut 4.0 --growthCut 2.0 --residual --fits
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
    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".selavy.residual.model.fits",f"{output_dir}/{image_filename}")
    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".selavy.residual.fits",f"{output_dir}/{image_filename}")
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
        # compute (ra,dec) in pixel units
        ra  = (row['ra_deg_cont']*eval(f"u.{catalogue_header['ra_deg_cont']}")).to(u.deg).value
        dec = (row['dec_deg_cont']*eval(f"u.{catalogue_header['dec_deg_cont']}")).to(u.deg).value
        r_pix = wcs.all_world2pix([[ra,dec]],0,ra_dec_order=True)[0]
        ra_0  = r_pix[0]
        dec_0 = r_pix[1]

        # compute semi-major/minor axes in pixel units
        pa = u.Quantity(row['pos_ang'],catalogue_header['pos_ang']).to(u.rad).value
        a  = u.Quantity(row['maj_axis'],catalogue_header['maj_axis']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        b  = u.Quantity(row['min_axis'],catalogue_header['min_axis']).to(u.deg).value/(2.0*np.sqrt(2.0*np.log(2.0)))
        smaj = a*np.sqrt((np.sin(pa)/CDELT1)**2+(np.cos(pa)/CDELT2)**2)
        smin = b*np.sqrt((np.cos(pa)/CDELT1)**2+(np.sin(pa)/CDELT2)**2)

        # get peak flux
        pf = u.Quantity(row['flux_peak'],catalogue_header['flux_peak']).to(u.Jy/u.beam).value
        
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

@click.command()
@click.argument('input_dir',nargs=1)
@click.argument('processing_dir',nargs=1)
@click.argument('output_dir',nargs=1)
@click.argument('fits_image_file',nargs=1)
@click.option(
    '--snrCut',
    default = selavy_pars['Selavy']['snrCut'],
    type = float,
    show_default = True,
    help = "Threshold value (in \u03C3's) above \u00B5-noise."
)
@click.option(
    '--growthCut',
    default = selavy_pars['Selavy']['growthCut'],
    type = float,
    show_default = True,
    help = "SNR to grow detections down to."
)
@click.option('--box-size', type=int,default=None,help="Grid RMS Box Size.")
@click.option('--step-size',type=int,default=None,help="Grid RMS Step Size (dummy: not used).")
@click.option('--fits',is_flag=True,default=False,help="Output FITS catalogue. [default: CSV]")
@click.option('--residual',is_flag=True,default=False,help="Output residual and model FITS files.")
@click.option('--dump',is_flag=True,default=False, help="Dump out all processing files.")
def process(
    input_dir,
    processing_dir,
    output_dir,
    fits_image_file,
    snrcut,
    growthcut,
    box_size,
    step_size,
    fits,
    residual,
    dump
):
    """\b
       Selavy image processing tool.

       inputs:

          \b
          INPUT_DIR: location of image_filename.fits
          PROCESSING_DIR: location of scratch directory
          OUTPUT_DIR: location to place results
          FITS_IMAGE_FILE: image_filename.fits (without path)

       outputs:

          \b
          OUTPUT_DIR/image_filename.aegean.fits
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
    print(f">    --snrCut: {snrcut}")
    print(f">    --growthCut: {growthcut}")

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

    # setup selavy configurtion
    cfg =  get_config(fits_image_file,snrcut,growthcut,box_size,step_size)
    #if re.search('Selavy.VariableThreshold.boxSize',cfg):
    #    cfg=re.sub("(Selavy.VariableThreshold.boxSize)", r"Selavy.VariableThreshold = true\n\1",cfg)
    cfg_file = f"{processing_dir}/selavy.in"
    print(f"Saving Configuration: {cfg_file}")
    print("> "+"\n> ".join(re.sub(r"\n$","",cfg).split("\n")))
    with open(cfg_file,"w") as f:
        f.write(cfg)

    # run selavy
    selavy = f"selavy -c selavy.in"
    print(f"> Executing selavy:")
    print(f"> $ {selavy}")
    sys.stdout.flush()
    current_dir = os.getcwd()
    os.chdir(processing_dir)
    os.system(selavy)
    os.chdir(current_dir)
    
    # convert txt to csv files
    print(f"> Creating .csv files...")
    for f in glob.glob(f"{processing_dir}/*.txt"):
        if not re.search(r"Logfile",f):
            csv_file = selavy_txt_to_csv(f)
            if re.search(".components.",csv_file):
                catalogue_file = csv_file
                if fits:
                    fits_file = re.sub(r"\.csv$",".fits",csv_file)
                    # NB: below bails when handling selavy-results.csv: i.e., error message,
                    #        ValueError: Arguments "names" and "dtype" must match number of columns
                    qt = Table.read(csv_file)
                    Table.read(csv_file).write(fits_file,format='fits',overwrite=True)
                    cleenup_fits_file(fits_file)

    print(f"> selavy files: {processing_dir}/")
    for f in glob.glob(f"{processing_dir}/*"):
        f = re.sub(r'^(.*?/)*','',f)
        if f != fits_image_file:
            print(f"> o {f}")

    # copy result from process to data dir
    print(f"> Transfering results:")
    fits_prefix = re.sub(r"\.([Ff][Ii][Tt]([Ss]|))$","",fits_image_file)
    for f in glob.glob(f"{processing_dir}/*"):
        s_file = re.sub(r"(.*?/)*","",f)
        f0 = s_file
        if s_file != fits_image_file:
           src = f"{processing_dir}/{s_file}"
           if s_file == 'selavy.in':
               s_file = f"{fits_prefix}.selavy.in"
           elif re.search(r"selavy-results\.(...)$",s_file):
               s_file = re.sub(r"selavy-results\.(...)$",r"%s.selavy.duchamp.\1" % fits_prefix,s_file)
           elif re.search(r"component[MR]",s_file):
               s_file = re.sub(r"^(component[MR].*?)_(.*?)(\.([Ff][Ii][Tt]([Ss]|)))$",r"\2.selavy.\1\3",s_file)
           elif s_file  == 'selavy-Logfile.txt':
               s_file = f"{fits_prefix}.selavy.log"
           else:
               s_file = re.sub(r"^selavy-results",f"{fits_prefix}.selavy",s_file)
               f0 = s_file
           s_file = re.sub(r"(selavy)\.components\.(...(|.))$",r"\1.\2",s_file)
           if not dump and not re.search(r"\.selavy\.(reg|%s)$" % ('fits' if fits else 'csv'),s_file):
               continue
           if re.search(r"\.reg$",src):
               cleanup_reg_file(src)
           dst = f"{output_dir}/{s_file}"
           print(f"> $ cp {src} {dst}")
           copy(src,dst)


    if residual:
        #create_residual_image(
        #    re.sub(r"^(.*?/)+","",fits_image_file),
        #    re.sub(r"^(.*?/)+","",fits_file if fits else catalogue_file),
        #    processing_dir,
        #    output_dir
        #)
        print(f"> Fetching resdual/model images.")
        image_filename = re.sub(r"^(.*?/)+","",fits_image_file)
        model_src = f"{processing_dir}/componentMap_{image_filename}"
        model_dst = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".selavy.residual.model.fits",f"{output_dir}/{image_filename}")
        print(f"> $ cp {model_src} {model_dst}")
        copy(model_src,model_dst)
        residual_src = f"{processing_dir}/componentResidual_{image_filename}"
        residual_dst = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".selavy.residual.fits",f"{output_dir}/{image_filename}")
        print(f"> $ cp {residual_src} {residual_dst}")
        copy(residual_src,residual_dst)

    print(f"[Done]")


if __name__ == "__main__":
    process()
