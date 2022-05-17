import re
import sys
import time
import alphashape
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy import units as u

def get_boundary(x,y,alpha=1.0):
    # guard condition
    if len(x) < 4: 
        return {'x': x, 'y': y}

    # compute the boundry points
    is_hard = False
    pts = [(a,b) for a,b in zip(x,y)]
    try: # suggested alpha value
        x,y = alphashape.alphashape(pts,alpha).exterior.coords.xy
    except: # whoops.. use brute force
        x,y = alphashape.alphashape(pts).exterior.coords.xy
        is_hard = True

    return {'x': x.tolist(), 'y': y.tolist(), 'is_hard': is_hard}

class Pix2World:
    def __init__(self,fits_image_file):
        hdul = fits.open(fits_image_file)
        self.header = hdul[0].header

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
            if key in self.header:
                self.header.insert(key,(rotation_matrix_map[key],self.header[key]),after=True)
                self.header.remove(key)

        # wcs image header
        self.wcs = WCS(self.header)

        # trim to 2d from nd
        naxis = self.wcs.naxis
        while naxis > 2:
            self.wcs = self.wcs.dropaxis(2)
            naxis -= 1

        # image data
        self.data = np.squeeze(hdul[0].data)


    def pix2world(self,p_x,p_y):
        ra  = list()
        dec = list()
        for x,y in zip(p_x,p_y):
            alpha, delta = self.wcs.all_pix2world(x,y,0,ra_dec_order=True)
            ra.append(alpha)
            dec.append(delta)
        return pd.DataFrame({'ra': ra, 'dec': dec})


class PixSegmentImage(Pix2World):
    def __init__(self,fits_image_file,segment_image,alphashape_alpha=1.0):
        # init wcs framework
        super().__init__(fits_image_file)

        # define main pars
        self.seg_ids = np.unique(segment_image[(segment_image>0)])
        self.seg_img = segment_image
        self.alpha = alphashape_alpha

        # define lazy-load parameters
        self.pix_segments   = None
        self.pix_boundaries = None
        self.time_segment_extractions  = None
        self.time_boundary_extractions = list()

    def get_extraction_stats(self):
        if len(self.time_boundary_extractions) > 0:
            units = self.time_boundary_extractions[0].unit
            dts = [dt.value for dt in self.time_boundary_extractions]
            stats = {
                'segment': self.time_segment_extractions,
                'poly_tot':     np.sum(dts ) * units,
                'poly_avg': np.mean(dts) * units,
                'poly_std': np.std(dts ) * units
            }
        else:
            stats = {
                'segment': self.time_segment_extractions,
                'poly_tot': None,
                'poly_avg': None,
                'poly_std': None 
            }
        return stats

    def __get_pix_segments(self):
        if self.pix_segments is None:
            self.pix_segments = {id: {'x': list(), 'y': list(), 'flux': list()} for id in self.seg_ids}
            t_0 = time.time()
            print(f">> Extracting Image Segments:")
            # nb: we must associate x <=> row(i) and y <=> col(i) as profound flips the image mtx
            for i in range(self.seg_img.shape[0]):
                for j in range(self.seg_img.shape[1]):
                    if self.seg_img[i,j] > 0:
                        self.pix_segments[self.seg_img[i,j]]['x'].append(i)
                        self.pix_segments[self.seg_img[i,j]]['y'].append(j)
                        self.pix_segments[self.seg_img[i,j]]['flux'].append(self.data[j,i])
            self.time_segment_extractions = (time.time()-t_0)*u.s
            print(f">>> Extracted {len(self.seg_ids)} Segments in {self.time_segment_extractions}.")
            print(f">> [Done]")
            sys.stdout.flush()
        return self.pix_segments

    def __get_pix_boundaries(self):
        if self.pix_boundaries is None:
            segments = self.__get_pix_segments()
            def get_boundary(segment):
                x = segment['x']
                y = segment['y']
                is_hard = False
                if len(x) > 3: # compute the boundry points
                    pts = [(a,b) for a,b in zip(x,y)]
                    try: # suggested alpha value
                        x,y = alphashape.alphashape(pts,self.alpha).exterior.coords.xy
                    except: # whoops.. use brute force
                        x,y = alphashape.alphashape(pts).exterior.coords.xy
                        is_hard = True
                return {'x': x.tolist(), 'y': y.tolist(), 'is_hard': is_hard}
            cnt = 1
            n_ids = len(self.seg_ids)
            print(f">> Extracting Segment Boundaries:")
            print(f">>> Segments: {n_ids}")
            self.pix_boundaries = dict()
            dts = list()
            for seg_id in self.seg_ids:
                t_0 = time.time()
                self.pix_boundaries[seg_id] = get_boundary(segments[seg_id])
                pts = len(self.pix_boundaries[seg_id]['y'])
                dt = (time.time()-t_0)*u.s
                print(f">>> Processed: {cnt} / {n_ids} ({pts} pts, {dt})")
                sys.stdout.flush()
                self.time_boundary_extractions.append(dt)
                dts.append(dt)
                cnt += 1
                stats = self.get_extraction_stats()
            print(f">>> Time: Total={stats['poly_tot']}, AVG={stats['poly_avg'].value}\u00B1{stats['poly_std']}")
            print(f">> [Done]")
            sys.stdout.flush()
        return self.pix_boundaries

    def get_segment_ids(self):
        return self.seg_ids

    def get_segments(self):
        return self.__get_pix_segments()

    def get_segment(self,segment_id):
        segments = self.__get_pix_segments()
        if segment_id in segments:
            return segments[segment_id]
        return {'x': list(), 'y': list()}

    def get_boundaries(self):
        return self.__get_pix_boundaries()

    def get_boundary(self,segment_id):
        boundaries = self.__get_pix_boundaries()
        if segment_id in boundaries:
            return boundaries[segment_id]
        return {'x': list(), 'y': list(), 'is_hard': False}


class WorldSegmentImage(PixSegmentImage):
    def __init__(self,fits_image_file,segment_image,alphashape_alpha=1.0):
        super().__init__(fits_image_file,segment_image,alphashape_alpha)

        # define lazy-load parameters
        self.world_segments   = None
        self.world_boundaries = None

    def get_segments(self):
        if self.world_segments is None:
            self.world_segments = dict()
            for seg_id,pts in super().get_segments().items():
                world = self.pix2world(pts['x'],pts['y'])
                self.world_segments[seg_id] = {
                    'ra': [float(r) for r in world['ra']],
                    'dec': [float(d) for d in world['dec']],
                    'flux': pts['flux']
                }
        return self.world_segments

    def get_segment(self,segment_id):
        segments = self.get_segments()
        if segment_id in segments:
            return segments[segment_id]
        return {'x': list(), 'y': list()}

    def get_boundaries(self):
        if self.world_boundaries is None:
            self.world_boundaries = dict()
            for seg_id,pts in super().get_boundaries().items():
                world = self.pix2world(pts['x'],pts['y'])
                self.world_boundaries[seg_id] = {
                    'ra': [float(r) for r in world['ra']],
                    'dec': [float(d) for d in world['dec']],
                    'is_hard': pts['is_hard']
                }
        return self.world_boundaries

    def get_boundary(self,segment_id):
        boundaries = self.get_boundaries()
        if segment_id in boundaries:
            return boundaries[segment_id]
        return {'x': list(), 'y': list(), 'is_hard': False}


class ProFoundSegments(WorldSegmentImage):
    def __init__(self,fits_image_file,segment_image,alphashape_alpha=1.0):
        super().__init__(fits_image_file,segment_image,alphashape_alpha)

        # define lazy-load parameters
        self.segments   = None
        self.boundaries = None

    def get_segments(self):
        if self.segments is None:
            self.segments = list()
            segments = super().get_segments()
            for seg_id in segments:
                self.segments.append({
                    'segid': seg_id,
                    'ra': segments[seg_id]['ra'],
                    'dec': segments[seg_id]['dec'],
                    'flux': segments[seg_id]['flux']
                })
        return self.segments

    def get_boundaries(self):
        if self.boundaries is None:
            self.boundaries = list()
            boundaries = super().get_boundaries()
            for seg_id in boundaries:
                self.boundaries.append({
                    'segid': seg_id,
                    'ra': boundaries[seg_id]['ra'],
                    'dec': boundaries[seg_id]['dec'],
                    'is_hard': boundaries[seg_id]['is_hard']
                })
        return self.boundaries


def csv_to_fits(csv_file,units,segments=list(),polygons=list()):
    def pythonize(csv_file):
        with open(csv_file,"r+") as fd:
            contents = list()
            for line in fd:
                line = re.sub('NA','NaN',line)
                contents.append(line)
            fd.seek(0)
            fd.truncate(0)
            for line in contents:
                fd.write(line)
    t_0 = time.time()
    # first let's make ouR csv file Python compatible...
    pythonize(csv_file)
    # create and open base fits file
    fits_file = re.sub(r"\.csv$",".fits",csv_file)
    Table.read(csv_file).write(fits_file,format='fits',overwrite=True)
    hdul = fits.open(fits_file)
    
    # add units
    cols = hdul[1].columns
    for key,value in units.items():
        if len(value) > 0:
            cols.change_attrib(key,'unit',value)
    hdul.writeto(fits_file,overwrite=True)

    # add polygons
    if len(polygons) > 0:
        # Notes: https://docs.astropy.org/en/stable/io/fits/usage/unfamiliar.html#variable-length-array-tables
        data = hdul[1].data

        # add boudary info
        #polygons = {p['segid']:{'ra':[c[0] for c in p['ra']],'dec':[c[0] for c in p['dec']]} for p in polygons}
        polygons = {p['segid']:{'ra':p['ra'],'dec':p['dec']} for p in polygons}
        ra  = [polygons[d['segid']]['ra']  for d in data] 
        dec = [polygons[d['segid']]['dec'] for d in data] 
        max_len = max([len(polygons[segid]['ra' ]) for segid in polygons])
        cols.add_col(fits.Column(name='polygon_ra', unit='deg',format=f"PD({max_len})",array=np.array(ra, dtype=np.object_)))
        cols.add_col(fits.Column(name='polygon_dec',unit='deg',format=f"PD({max_len})",array=np.array(dec,dtype=np.object_)))

        # add segment info
        segments = {s['segid']:{'ra':s['ra'],'dec':s['dec'],'flux':s['flux']} for s in segments}
        ra   = [segments[d['segid']]['ra']   for d in data] 
        dec  = [segments[d['segid']]['dec']  for d in data] 
        flux = [segments[d['segid']]['flux'] for d in data] 
        max_len = max([len(segments[segid]['ra' ]) for segid in segments])
        cols.add_col(fits.Column(name='segments_ra',  unit='deg',format=f"PD({max_len})",array=np.array(ra,  dtype=np.object_)))
        cols.add_col(fits.Column(name='segments_dec', unit='deg',format=f"PD({max_len})",array=np.array(dec, dtype=np.object_)))
        cols.add_col(fits.Column(name='segments_flux',unit='Jy', format=f"PD({max_len})",array=np.array(flux,dtype=np.object_)))
    
        # update base fits file
        fits.BinTableHDU.from_columns(cols).writeto(fits_file,overwrite=True)
    print(f"CSV2FITS Conversion Time: {(time.time()-t_0)*u.s}")


def create_residual_image(segment_image,image_filename,input_dir,output_dir):
    t_0 = time.time()
    image_file     = f"{input_dir}/{image_filename}"
    model_file     = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".profound.residual.model.fits",f"{output_dir}/{image_filename}")
    residual_file  = re.sub(r"\.[Ff][Ii][Tt](|[Ss])$",".profound.residual.fits",f"{output_dir}/{image_filename}")
    print(f"> Creating Residual Image:")
    print(f">>    INPUT_IMAGE: {image_file}")

    # open image
    hdul = fits.open(image_file)
    header = hdul[0].header
    shape = hdul[0].data.shape
    image = np.squeeze(hdul[0].data)

    # create residual image
    model = segment_image.transpose() if not segment_image is None else np.zeros(image.shape)
    image[(model>0)] = 0.0

    # output model image
    print(f">>    OUTPUT_MODEL: {model_file}")
    model.shape = shape
    fits.PrimaryHDU(model,header=header).writeto(model_file,overwrite=True)

    # output residual image
    print(f">>    OUTPUT_IMAGE: {residual_file}")
    image.shape = shape
    fits.PrimaryHDU(image,header=header).writeto(residual_file,overwrite=True)
    print(f"Residual Image Creation Time: {(time.time()-t_0)*u.s}")
    print(f"> [Done]")





