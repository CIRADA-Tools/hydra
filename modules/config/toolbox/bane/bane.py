###############################################################################################
# #    * * *   B A N E   C O N T A I N E R   P R O C E S S I N G   S C R I P T    * * *
#

import re
import os
import sys
import glob
import click
import numpy as np
from pathlib import PosixPath
from shutil import copy2 as copy
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

class FITS:
    def __init__(self,fits_image_file):
        try:
            if isinstance(fits_image_file,str) or isinstance(fits_image_file,PosixPath):
                def get_calling_scipt_path_reference(filepath):
                    fp = filepath if os.path.isfile(filepath) else os.path.abspath(f"{os.getcwd()}/{filepath}")
                    if not os.path.isfile(fp):
                        print(f"WARNING: File '{filepath}' not found!")
                    return fp
                def repair_header(header):
                    # this is to overcome a problem with bane, which sometimes
                    # puts out flawed headers...
                    if WCS(header).naxis == 4 and header['NAXIS'] != 4:
                        header.update({
                            'NAXIS': 4,
                            'NAXIS3': 1,
                            'NAXIS4': 1
                        })
                    return header
                self.file = get_calling_scipt_path_reference(fits_image_file)
                hdul = fits.open(self.file)
                self.header = repair_header(hdul[0].header)
                self.shape  = hdul[0].data.shape
                self.data   = np.squeeze(hdul[0].data)
            else: # assume hdu
                self.header = fits_image_file.header
                self.shape  = fits_image_file.data.shape
                self.data   = np.squeeze(fits_image_file.data)
        except Exception as e:
            print(f"ERROR: {e}")
            print(f"Bye!")
            exit()

    def get_header(self):
        return self.header

    def get_wcs(self):
        return get_wcs(self.get_header())

    def get_data(self,is_flatten_clean=False):
        if not is_flatten_clean:
            return self.data
        data = np.ndarray.flatten(self.data)
        return data[~np.isnan(data)]

    def get_mean(self):
        return np.nanmean(self.data)

    def get_mean_over_rms(self):
        return np.nanmean(self.data)/np.nanstd(self.data)

    def get_rms(self):
        return np.nanstd(self.data)

    def get_median(self):
        data = self.get_data(is_flatten_clean=True)
        return np.nanmedian(data)

    def get_median_over_madfm(self):
        data = self.get_data(is_flatten_clean=True)
        median = np.nanmedian(data)
        meddev = np.abs(data-median)
        madfm = np.median(meddev)/0.6744888
        return median/madfm

    def get_madfm(self):
        data = self.get_data(is_flatten_clean=True)
        median = np.nanmedian(data)
        meddev = np.abs(data-median)
        return np.median(meddev)/0.6744888

    def get_sumsq(self):
        data = self.get_data(is_flatten_clean=True)
        return np.sum(data**2)

    def get_field(self,field):
        return self.header[field] if self.has_field(field) else None

    def has_field(self,field):
        return field in self.header
    
    def has_beam_shape(self):
        return self.has_field('BMAJ') and self.has_field('BMIN') and self.has_field('BPA')

def optimize_box_pars(fits_file):
    # file i/o names
    rms_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$","_rms.fits",fits_file)
    print(f"Optimizing: {fits_file}")
    print(f"> rms_file: {rms_file}")

    # get initial pars
    fits = FITS(fits_file)
    if fits.has_beam_shape():
        bmaj = fits.get_field('BMAJ')/np.abs(fits.get_field('CDELT2'))
        bmin = fits.get_field('BMIN')/np.abs(fits.get_field('CDELT2'))
        bpa  = fits.get_field('BPA')
    else:
        print(f"ERROR: Require Beam Parameters: BMAJ, BMIN, PA")
        print(f"Bye!")
        exit()
    print(f"> Beam Pars:")
    print(f"> o BMAJ: {bmaj}")
    print(f"> o BMIN: {bmin}")
    print(f"> o   PA: {bpa}")
    step_size = 4.0*(bmaj+bmin)/2.0
    # 3s <= box_size <= 6s seems ok
    box_min = 3.0*step_size
    box_max = 6.0*step_size
    ##box_max = 8.0*step_size
    #box_max = 10.0*step_size

    # run over pars
    # n_boxes=n_sizes=3 seems ok
    stats = list()
    n_sizes = 3
    n_boxes = 4
    ##n_boxes = 6
    #n_boxes = 8
    d_box = (box_max-box_min)/np.float(n_boxes-1)
    for i in range(n_boxes):
        box_size = np.int64(i*d_box+box_min)
        step_min = box_size/4.0
        step_max = box_size/2.0
        d_step = (step_max-step_min)/(n_sizes-1.0)
        for j in range(n_sizes):
            step_size = np.int64(j*d_step+step_min)
            print(f"> Doing: Box_Size={box_size}, Step_Size={step_size}")
            bane = f"BANE {fits_file} --cores 1"
            bane += f" --grid {step_size} {step_size}" # 16 (step size)
            bane += f" --box {box_size} {box_size}" # 96
            print(f"> $ {bane}")
            sys.stdout.flush()
            os.system(bane)
            fits = FITS(rms_file)
            mean            = fits.get_mean()
            median          = fits.get_median()
            rms             = fits.get_rms()
            madfm           = fits.get_madfm()
            sumsq           = fits.get_sumsq()
            mean_over_rms   = fits.get_mean_over_rms()
            median_over_madfm = fits.get_median_over_madfm()
            print(f"> \u03BC={mean}, Median={median}, \u03C3={rms}, MADFM={madfm}, \u03A3I\u00B2={sumsq}, \u03BC/\u03C3={mean_over_rms}, Median/MADFM={median_over_madfm}")
            stats.append({
                'box_size':        box_size,
                'step_size':       step_size,
                'mean':            mean,
                'median':          median,
                'rms':             rms,
                'madfm':           madfm,
                'sumsq':           sumsq,
                'mean_over_rms':   mean_over_rms,
                'median_over_madfm': median_over_madfm,
            })

    # get opitimal pars
    columns = {
        'id':                {'dtype': np.int64},
        'is_opt':            {'dtype': np.int64},
        'box_size':          {'dtype': np.int64},
        'step_size':         {'dtype': np.int64},
        'mean':              {'dtype': np.float64, 'opt_id': 1, 'label': '\u03BC'},
        'median':            {'dtype': np.float64, 'opt_id': 2, 'label': 'Median'},
        'rms':               {'dtype': np.float64, 'opt_id': 3, 'label': '\u03C3'},
        'madfm':             {'dtype': np.float64, 'opt_id': 4, 'label': 'MADFM'},
        'sumsq':             {'dtype': np.float64, 'opt_id': 5, 'label': '\u03A3I\u00B2'},
        'mean_over_rms':     {'dtype': np.float64, 'opt_id': 6, 'label': '\u03BC/\u03C3'},
        'median_over_madfm': {'dtype': np.float64, 'opt_id': 7, 'label': 'Median/MADFM'},
    }
    metric = 'mean'
    print(f"> Stats:")
    qt = Table(
        names=columns,
        dtype=[columns[col]['dtype'] for col in columns]
    )
    idx = 1
    opt_value = None
    pars_opt  = None
    opt_idx   = None
    for datum in stats:
        print("> "+", ".join([f"{item.upper()}={datum[item]}" for item in datum]))
        if opt_value is None:
            opt_value  = datum[metric]
            pars_opt = datum
            opt_idx = idx
        elif opt_value > datum[metric]:
            opt_value = datum[metric]
            pars_opt = datum
            opt_idx = idx
        qt.add_row([idx if c=='id' else (0 if c=='is_opt' else datum[c]) for c in columns])
        idx += 1
    qt['is_opt'][(qt['id']==opt_idx)] = columns[metric]['opt_id']
    print(qt)
    print(f"> Optimal Pars:")
    print(f"> METRIC: {columns[metric]['label']}")
    for par in pars_opt:
        print(f"> {par.upper()}: {pars_opt[par]}")
    print(f"[Done]")

    return {
        'metric': columns[metric]['label'],
        'opt_pars': pars_opt,
        'stats': qt,
    }


@click.command()
@click.argument('input_dir',nargs=1)
@click.argument('processing_dir',nargs=1)
@click.argument('output_dir',nargs=1)
@click.argument('fits_image_file',nargs=1)
@click.option('--optimize',is_flag=True,default=False,help="Optimize RMS Box (overrides box/step-size).")
@click.option('--box-size', type=int,default=-1,help="Grid RMS Box Size [default: ~4*beam_size]")
@click.option('--step-size',type=int,default=-1,help="Grid RMS Step Size [default: 5*box_size]")
def process(
    input_dir,
    processing_dir,
    output_dir,
    fits_image_file,
    optimize,
    box_size,
    step_size
):
    """\b
       BANE image background (bkg) and RMS (rms) generation tool.

       inputs:

          \b
          INPUT_DIR: location of image_filename.fits
          PROCESSING_DIR: location of scratch directory
          OUTPUT_DIR: location to place results
          FITS_IMAGE_FILE: image_filename.fits (without path)

       outputs:

          \b
          OUTPUT_DIR/image_filename.bane.bkg.fits
          OUTPUT_DIR/image_filename.bane.rms.fits
    """

    src = f"{input_dir}/{fits_image_file}"
    print(f"Processing: {src}")
    print(f"> Local Context: ")
    print(f">   INPUT_DIR:  {input_dir}")
    print(f">   PROCESSING_DIR: {processing_dir}")
    print(f">   OUTPUT_DIR: {output_dir}")

    # check if image file exists
    if not os.path.isfile(src):
        print(f"ERROR: Image file '{src}' not found!")
        print(f"Bye!")
        exit()

    # link image file from data to processing dir
    dst = f"{processing_dir}/{fits_image_file}"
    print(f"> Linking image file:")
    print(f"> $ ln -s {dst} {src}")
    # this gaurd condition is usefull for debugging.
    if not (os.path.islink(dst) or os.path.exists(dst)):
        os.symlink(src,dst)

    # setup bane cmd
    bane = f"BANE {dst} --cores 1"

    # set/optimize bane box/step-size if required
    if optimize:
        # compute stats and opt rms box pars
        opt_stats = optimize_box_pars(dst)

        # save stats
        stats_file = output_dir+"/"+re.sub("\.[Ff][Ii][Tt]([Ss]|)$",".rms_box_statistcs.csv",fits_image_file)
        print(f"Saving RMS Box Statistics: {stats_file}")
        opt_stats['stats'].write(stats_file,format='csv',overwrite=True)

        # set box/step-size
        step_size = opt_stats['opt_pars']['step_size']
        box_size = opt_stats['opt_pars']['box_size']
    elif step_size > 0 and box_size > 0 and step_size > box_size:
        print(f"ERROR: step_size > box_size")
        print(f">  box_size: {box_size}")
        print(f"> step_size: {step_size}")
        print(f"Bye!")
        exit()
    if step_size > 0:
        bane += f" --grid {step_size} {step_size}" # 16 
    if box_size > 0:
        bane += f" --box {box_size} {box_size}" # 96

    # execute bane
    print(f"> Executing BANE:")
    print(f"> $ {bane}")
    sys.stdout.flush()
    os.system(bane)
    print(f"> BANE files: {processing_dir}/")
    files = [re.sub(r'^(.*?/)*','',f) for f in glob.glob(f"{processing_dir}/*")]
    for f in files:
        if f != fits_image_file:
            print(f">   o {f}")

    rms_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$","_rms.fits",fits_image_file)
    bkg_file = re.sub("\.[Ff][Ii][Tt]([Ss]|)$","_bkg.fits",fits_image_file)
    src = f"{processing_dir}/{rms_file}"
    if optimize:
        fits = FITS(src)
        print(f" metric: {opt_stats['metric']}")
        print(f" --------{'-'*len(opt_stats['metric'])}")
        print(f"  BOX_SIZE: {box_size}")
        print(f" STEP_SIZE: {step_size}")
        print(f"      MEAN: {fits.get_mean()}")
        print(f"       RMS: {fits.get_rms()}")
        print(f"     MADFM: {fits.get_madfm()}")
        print(f"     SumSq: {fits.get_sumsq()}")
    dst = "{}/{}".format(output_dir,re.sub(r"_(rms\.fits)$",r".bane.\1",rms_file))
    print(f"> $ cp {src} {dst}")
    copy(src,dst)
    src = f"{processing_dir}/{bkg_file}"
    dst = "{}/{}".format(output_dir,re.sub(r"_(bkg\.fits)$",r".bane.\1",bkg_file))
    print(f"> $ cp {src} {dst}")
    copy(src,dst)

    # done
    print(f"[Done]")

if __name__ == "__main__":
    process()
