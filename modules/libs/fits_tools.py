import os
import re
import sys
from pathlib import Path
from pathlib import PosixPath
from shutil import copy2 as copy
sys.path.insert(0,"../")
from libs.exceptions import print_warning
from libs.exceptions import print_ok
from libs.exceptions import HydraContainerError
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.convolution import convolve
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from astropy import units as u
from warnings import simplefilter as warning_filter
from astropy.io.fits.verify import VerifyWarning
warning_filter('ignore',category=VerifyWarning)

# ticky business this!
def _get_calling_scipt_path_reference(filepath):
    fp = filepath if os.path.isfile(filepath) else os.path.abspath(f"{os.getcwd()}/{filepath}")
    if not os.path.isfile(fp):
        print(f"WARNING: File '{filepath}' not found!")
    return fp


# gets fits image data from file, np-array, hdul, or hdu and reduces dimension
def get_fits_image_data(fits_object):
    data = None
    if isinstance(fits_object, str) and re.search(r"\.fit(s|)(|\.gz)$",fits_object.lower()) and os.path.isfile(fits_object):
        data = fits.open(_get_calling_scipt_path_reference(fits_object))[0].data
    elif isinstance(fits_object, np.ndarray):
        data = fits_object
    elif isinstance(fits_object,fits.HDUList):
        data = fits_object[0].data
    elif isinstance(fits_object,fits.PrimaryHDU):
        data = fits_object.data
    return None if data is None else np.squeeze(data)


# general fits file handler
class FITS:
    def __init__(self,fits_image_file):
        self.file = None
        hdul = None
        try:
            if isinstance(fits_image_file,str) or isinstance(fits_image_file,PosixPath):
                def repair_header(header):
                    try:
                        # this is to overcome a problem with bane, which sometimes
                        # puts out flawed headers...
                        if WCS(header).naxis == 4 and header['NAXIS'] != 4:
                            header.update({
                                'NAXIS': 4,
                                'NAXIS3': 1,
                                'NAXIS4': 1
                            })
                    except Exception as e:
                        # kludge: this exception block was put in as there appears to be
                        # random problems with caesar's residual image header output,
                        # re. WCS: i.e.,
                        #    WARNING: FITSFixedWarning: 'celfix' made the change 'Unmatched celestial axes'. [astropy.wcs.wcs]
                        #    ERROR: ERROR 4 in wcs_types() at line 2707 of file cextern/wcslib/C/wcs.c:
                        print_warning(f"ERROR: {e}")
                        print_warning(f"Issues with header: {fits_image_file}")
                        print(f" * * *   H E A D E R   * * *")
                        print(header)
                        print_warning(f"WARNING: May cause problems downstream.")
                        print_warning(f"Continuing...")
                    return header
                self.file = _get_calling_scipt_path_reference(fits_image_file)
                hdul = fits.open(self.file)
                self.header = repair_header(hdul[0].header)
                self.shape  = hdul[0].data.shape
                self.data   = np.squeeze(hdul[0].data)
            else: # assume hdu
                self.header = fits_image_file.header
                self.shape  = fits_image_file.data.shape
                self.data   = np.squeeze(fits_image_file.data)
        except Exception as e:
            print_warning(f"ERROR: {e}")
            if not hdul is None:
                print(f" * * *   H E A D E R   * * *")
                print(hdul[0].header)
                print(f" * * *   D A T A   * * *")
                print(hdul[0].data)
            else:
                print(f"hdul: {hdul}")
            self.header = None
            self.data   = None
            if not self.file is None:
                print(f"Saving copy: {self.file}")
                print(f"> $ cp {self.file} {self.file}.error")
                copy(self.file,f"{self.file}.error")
            print_warning(f"Continuing...")


    def get_header(self):
        return self.header

    def get_wcs(self):
        return get_wcs(self.get_header())

    def get_data(self,is_flatten_clean=False):
        if not is_flatten_clean:
            return self.data
        data = np.ndarray.flatten(self.data)
        return data[~np.isnan(data)]

    def get_rms(self):
        return np.nanstd(self.data)

    def get_madfm(self):
        data = self.get_data(is_flatten_clean=True)
        median = np.nanmedian(data)
        meddev = np.abs(data-median)
        return np.median(meddev)/0.6744888

    def get_sumsq(self):
        data = self.get_data(is_flatten_clean=True)
        return np.sum(data**2)

    def get_header_and_data(self,is_flatten_clean=False):
        return [self.get_header(),self.get_data(is_flatten_clean)]

    def get_field(self,field):
        return self.header[field] if self.has_field(field) else None

    def has_field(self,field):
        return field in self.header
    
    def has_beam_shape(self):
        return self.has_field('BMAJ') and self.has_field('BMIN') and self.has_field('BPA')

    def get_survey(self):
        return self.header['SURVEY'] if self.has_field('SURVEY') else None

    def has_frequency(self):
        naxis = self.get_field('NAXIS')
        if not naxis is None and naxis > 2:
            for i in range(naxis-2):
                ctype = f"CTYPE{i+3}"
                if self.has_field(ctype):
                    value = re.sub("\s+","",self.get_field(ctype))
                    if re.search(r"^fr(e(q(u(e(n(c(y|)|)|)|)|)|)|)$",value.lower()):
                            return True
        return False

    def get_frequency(self):
        naxis = self.get_field('NAXIS')
        if not naxis is None and naxis > 2:
            for i in range(naxis-2):
                ctype = f"CTYPE{i+3}"
                if self.has_field(ctype):
                    value = re.sub("\s+","",self.get_field(ctype))
                    if re.search(r"^fr(e(q(u(e(n(c(y|)|)|)|)|)|)|)$",value.lower()):
                        frequency = self.get_field(f"CRVAL{i+3}")
                        try:
                            units = eval("u."+re.sub(r"\s+","",self.get_field(f"CUNIT{i+3}")))
                        except:
                            units = u.Hz
                        return frequency * units if not frequency is None else None
        return None

    def write(self,fits_filename):
        image = self.data
        image.shape = self.shape
        fits.PrimaryHDU(image,header=self.header).writeto(fits_filename,overwrite=True)


# fits image stats tool
# notes: https://www.gnu.org/software/gnuastro/manual/html_node/Sigma-clipping.html
def sigma_clipper(fits_file,alpha=3.25,iters=5,tol=.01):
    image = get_fits_image_data(fits_file)
    if image is None:
        return None
    image = image.flatten()
    image = image[~np.isnan(image)] # get rid of nan's
    initial_max = np.max(image)
    initial_min = np.min(image)
    sigmas = list()
    for i in range(iters):
       median = np.median(image)
       sigma =  np.std(image)
       clip_plus  = median+alpha*sigma
       clip_minus = median-alpha*sigma
       image = image[(clip_minus < image) & (image < clip_plus)]
       sigmas.append(sigma)
       if np.abs(sigma-np.std(image))/np.std(image) < tol:
           break
    final_max = np.max(image)
    final_min = np.min(image)
    sigmas.append(np.std(image))
    pdfs = [(sigmas[i]-sigmas[i+1])/sigmas[i+1] for i in range(len(sigmas)-1)]
    return {
        'mean':   np.mean(image),
        'median': np.median(image),
        'rms':    np.std(image),
        'initial_min': initial_min,
        'initial_max': initial_max,
        'final_min': final_min,
        'final_max': final_max,
        'sigmas': sigmas,
        'pdfs':   pdfs,
        'iters':  i+1
    }


# shallow image creator
class Homados:
    def __init__(self,fits_file,output_dir=None,verbose=False):
        self.deep_file    = fits_file
        self.noise_file   = re.sub(r"(\.fits)(|\.gz)$",r".noise\1",fits_file)
        self.shallow_file = re.sub(r"(\.fits)(|\.gz)$",r".shallow\1",fits_file)
        if not output_dir is None:
            if os.path.isdir(output_dir):
                self.noise_file   = os.path.abspath(output_dir+"/"+re.sub(r"^(.*?/)*","",self.noise_file  ))
                self.shallow_file = os.path.abspath(output_dir+"/"+re.sub(r"^(.*?/)*","",self.shallow_file)) 
            else:
                print_warning("WARNING: Dir. '{}' invalid: Using dir. '{}'".format(
                    output_dir,
                    re.sub(r'^((.*?/)*).*$',r'\1',self.deep_file)
                ))

        self.verbose = verbose

    def __print(self,msg):
        if self.verbose:
            print(msg)

    def print_summary(self):
        if os.path.isfile(self.deep_file):
            rms_deep = sigma_clipper(self.deep_file)['rms']*10**6
            print(f"rms_deep: {rms_deep:.0f} \u03bcJy")
        else:
            rms_deep = None
    
        if os.path.isfile(self.noise_file):
            rms_noise = sigma_clipper(self.noise_file)['rms']*10**6
            print(f"rms_noise: {rms_noise:.0f} \u03bcJy")
            if not rms_deep is None:
                print(f"rms_noise/rms_deep: {rms_noise/rms_deep:.1f}")
    
        if os.path.isfile(self.shallow_file):
            rms_shallow = sigma_clipper(self.shallow_file)['rms']*10**6
            print(f"rms_shallow: {rms_shallow:.0f} \u03bcJy")
            if not rms_deep is None:
                print(f"rms_shallow/rms_deep: {rms_shallow/rms_deep:.1f}")

        return self

    def make_noise(self,thresh=5,beam_pars=None,is_save_noise_file=True):
        # compute image mean and rms
        self.__print(f"Processing {self.deep_file} ...")
        self.__print(f"> Computing \u03bc\u00B1\u03C3")
        clip_stats = sigma_clipper(self.deep_file)
        if clip_stats is None:
            print_warning("Whoops! This can't be a fits file.")
            print_ok("Bye!")
            return self
        mean = clip_stats['mean']
        noise_floor = thresh*clip_stats['rms']
        self.__print(f"> \u03bc\u00B1\u03C3 = {mean*10**6:.5}\u00B1{clip_stats['rms']*10**6:.5} \u03bcJy")
        
        # determine beam paramaters
        hdu   = fits.open(self.deep_file)
        image = hdu[0].data
        header = hdu[0].header
        for field in beam_pars:
            if not beam_pars[field] is None:
                header.update({field: beam_pars[field]})
        
        for i in range(2):
            # compute gaussian noise
            self.__print(f"> Iter{i+1}/2: Noise floor: n\u03C3 = {noise_floor*10**6:.5} \u03bcJy")
            self.__print(f"> Iter{i+1}/2: Computing noise image...")
            noise = np.random.normal(mean,noise_floor,np.squeeze(image).shape)
            
            # convolve noise with gaussian kernal
            cdeln = 2.0*np.abs(header['CDELT1'])
            self.__print("> Iter{0}/2: Convolving with Gaussian 2D Kernel: {{BMIN: {1:.3}\", DMAJ: {2:.3}\", BPA: {3:.3}\u00B0}}".format(
                i+1,
                (header['BMIN']*u.deg).to(u.arcsec).value,
                (header['BMAJ']*u.deg).to(u.arcsec).value,
                header['BPA']
            ))
            kernel = Gaussian2DKernel(
                x_stddev=header['BMIN']/cdeln,
                y_stddev=header['BMAJ']/cdeln,
                theta=np.deg2rad(header['BPA'])
            )
            convd = convolve(noise,kernel)
            convd.shape = image.shape
            
            if i == 0:
                # magic: new noise_floor on 2nd iteration corrects noise-level to thresh
                noise_thresh = sigma_clipper(noise)['rms'] # noise rms
                convd_thresh = sigma_clipper(convd)['rms'] # convolved rms
                noise_floor = noise_thresh*thresh*clip_stats['rms']/convd_thresh
    
        # save noise image
        if is_save_noise_file:
            self.__print(f"> Saving noise image: {self.noise_file}")
            noise_image = fits.PrimaryHDU(convd,header=hdu[0].header)
            noise_image.writeto(self.noise_file,overwrite=True)
    
        # create and save shallow image
        self.__print(f"> Saving shallow image: {self.shallow_file}")
        shallow_image = fits.PrimaryHDU(image+convd,header=hdu[0].header)
        shallow_image.writeto(self.shallow_file,overwrite=True)
        self.__print("[Done]")

        # print summary
        self.print_summary()

        return self

# fits image inverter
def invert_fits_image(fits_file,output_dir=None):
    inverted_fits_file = None
    if os.path.isfile(fits_file):
        print(f"Inverting: {fits_file}")
        try:
            hdul = fits.open(fits_file)
            data = -np.squeeze(hdul[0].data)
            data.shape = hdul[0].data.shape
            inverted_image = fits.PrimaryHDU(data,header=hdul[0].header)
            inverted_fits_file = re.sub(r"\.fits",".inverted.fits",fits_file)
            if not output_dir is None:
                if os.path.isdir(output_dir):
                    inverted_fits_file = os.path.abspath(output_dir+"/"+re.sub(r"^(.*?/)*","",inverted_fits_file))
                else:
                    print_warning(f"ERROR: Invalid output dir: {output_dir}")
                    eixt()
            print(f"> Writing: {inverted_fits_file}")
            inverted_image.writeto(inverted_fits_file,overwrite=True)
        except:
            print_warning("Whoops! This can't be a fits file.")
            exit()
    else:
        print_warning(f"ERROR: File '{fits_file}' does not exist!")
        print_ok("Bye!")
        exit() 
    print(f"[Done]")
    return inverted_fits_file

# header cleaner to wcs func
def get_wcs(header):
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
    w = WCS(header)
    naxis = w.naxis
    while naxis > 2:
        w = w.dropaxis(2)
        naxis -= 1
    return w

# gets sample cutout from image /w size at center of image
def create_sample(input_file,output_file,size=(0.5*u.deg,0.5*u.deg)):

    # tile trimmer func
    # notes: https://docs.astropy.org/en/stable/nddata/utils.html
    def trim_tile(hdul, position, size):
        hdu = hdul[0]
        w = get_wcs(hdu.header)
        img_data = np.squeeze(hdu.data)
    
        # get the cutout
        stamp = Cutout2D(img_data, position, (size[1],size[0]) if isinstance(size,tuple) else size, wcs=w, mode='trim', copy=True)
    
        # create fits image
        hdu.header.update(stamp.wcs.to_header())
        hdu.header.update({'WCSAXES': 4})
        stamp.data.shape = (1,1,stamp.data.shape[0],stamp.data.shape[1]) # kludge
        trimmed = fits.PrimaryHDU(stamp.data, header=hdu.header)
    
        return trimmed

    print(f"Extracting sample image...")

    # get input target parameters
    hdul = fits.open(input_file)
    hdr = hdul[0].header
    width  = hdr['NAXIS1']*abs(hdr['CDELT1'])*u.deg
    height = hdr['NAXIS2']*abs(hdr['CDELT2'])*u.deg
    input_area = width*height
    print(f"> Input Target: {input_file}")
    print(f">> WIDTH: {width}")
    print(f">> HEIGHT: {height}")
    print(f">> AREA: {input_area}")
    
    # get output target parameters
    w = get_wcs(hdr)
    output_area = size[0].to(u.deg)*size[1].to(u.deg)
    if input_area > output_area:
        if size[0] > width:
            size = (width,output_area/width)
        elif size[1] > height:
            size = (output_area/height,height)
    print(f"> Output Target: {output_file}")
    print(f">> WIDTH: {size[0]}")
    print(f">> HEIGHT: {size[1]}")
    print(f">> AREA: {output_area}")
    
    # extract sample image
    skirt = (1.2*u.arcmin).to(u.deg) # add small cutout selection threshold
    output_skirt_area = (size[0].to(u.deg)+skirt)*(size[1].to(u.deg)+skirt)
    if input_area > output_skirt_area: # create cutout
        print(f"> Input Target Area > Output Target Area:")
        ra,dec = w.all_pix2world(hdr['NAXIS1']/2.0,hdr['NAXIS2']/2.0,1)
        center = SkyCoord(ra,dec,unit=(u.deg,u.deg,))
        print(f">> Trimming Input Target: {input_file}")
        trimmed = trim_tile(hdul,center,size)
        print(f">> Creating Output Target: {output_file}")
        trimmed.writeto(output_file,overwrite=True)
    else: # whoops input target too small --> just copy
        print(f"> Input Target Area < Output Target Area:")
        print(f"> $ cp {input_file} {output_file}")
        copy(input_file,output_file)
    print("[Done]")
    
    return output_file

# create bane rms image from fits_image_file
def load_bane_rms_image(fits_image_file,cache,rms_box=None,is_optimize=False):
    print_ok(f"Computing Background Noise Image...")
    input_file = re.sub(r"^(.*?/)*","",fits_image_file)
    input_path = re.sub(r"^((.*?/)*).*$",r"\1",fits_image_file)
    if input_path == "":
        input_path = "."
    input_path  = Path(input_path).resolve()
    dcr_cmd = "docker run --rm -t -v {}:{} -v {}:{} {} {} {} {} {}".format(
        input_path, '/home/bane/data', # data in mount point
        cache,'/home/bane/results',    # data out mount point
        'bane',                        # service
        '/home/bane/data',             # input dir
        '/home/bane/processing',       # processing dir
        '/home/bane/results',          # output dir
        input_file                     # image file
    )
    if not rms_box is None and isinstance(rms_box,dict) and 'boxsize' in rms_box and 'stepsize' in rms_box:
        dcr_cmd += f" --box-size {rms_box['boxsize']}"
        dcr_cmd += f" --step-size {rms_box['stepsize']}"
    elif is_optimize:
        dcr_cmd += " --optimize"
    print_ok(f"> $ {dcr_cmd}")
    if os.system(dcr_cmd) != 0:
        print_warning(f"BANE Container Error")
        raise HydraContainerError
    bkg_file = Path(cache+"/"+re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".bane.bkg.fits",input_file)).resolve()
    rms_file = Path(cache+"/"+re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".bane.rms.fits",input_file)).resolve()
    print(f"RMS_IMG: {rms_file}")
    rms_img = FITS(rms_file)
    os.remove(bkg_file)
    os.remove(rms_file)
    print_ok(f"[Done]")
    return rms_img
