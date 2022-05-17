from cython cimport view

import re
import aplpy
import numpy as np
import matplotlib.pyplot as plt

# NOTES: The following is a try/except code snippet for error handling, by reporting,
# adding a dummy image, and continuing a la multiprocessing... for now, it has been
# decided to simply let it crash...
#
# [1] https://code-maven.com/create-images-with-python-pil-pillow
#from PIL import Image, ImageDraw, ImageFont
#    except Exception as e:
#        print_warning(f"   * * *   R E N D E R I N G   E R R O R   * * *")
#        print_warning(f"FNAME: {datum['fname']}")
#        print_warning(f"PNAME: {datum['pname']}")
#        print_warning(f"TASK: {datum}")
#        print_warning(f"Error: {e}")
#        img = Image.new('RGB',(1200,1200), color = (255,0,0))
#        d = ImageDraw.Draw(img)
#        fnt = ImageFont.truetype('/Library/Fonts/Arial.ttf', 45)
#        d.text((500,600), "Whoops!", font=fnt, fill=(255,255,255))
#        img.save(datum['pname'])
#        #raise HydraRenderingError 

def renderer(datum):
    print(f"> Rendering[{datum['counter']['cnt']}/{datum['counter']['max']}]: {datum['pname']}")
    fig = plt.figure(frameon=False,figsize=(12,12))
    if 'image' in datum:
        gc = aplpy.FITSFigure(datum['image'],figure=fig)
    else:
        gc = aplpy.FITSFigure(datum['fname'],figure=fig)
    gc.show_colorscale(cmap="Greys_r",stretch='linear')
    gc.set_title(datum['plt_title'])
    if 'annotation' in datum:
        annotation = datum['annotation']
        clump = annotation['clump']
        source_finders = list()
        for source_finder in np.unique(clump['source_finder']):
            df = clump[(clump['source_finder']==source_finder)]
            if len(df)>0:
                source_finders.append(source_finder)
                gc.show_ellipses(
                    df['ra'],df['dec'],2.0*df['extent_semiminor'],2.0*df['extent_semimajor'],df['extent_angle'],
                    edgecolor=datum['plt_colors'][source_finder],facecolor='None',linewidths=3,coords_frame='world',alpha=0.75
                )
        for bpar in annotation['bpars']:
            match_id = bpar['match_id']
            bbox = bpar['bbox']
            gc.add_label(
                bbox['postn'].ra.value,
                bbox['postn'].dec.value,
                f"{match_id}",
                size = 36
            )
            gc.show_rectangles(
                bbox['postn'].ra.value,
                bbox['postn'].dec.value,
                bbox['width'].value,
                bbox['height'].value,
                edgecolor='cyan',
                facecolor='None',
                linewidths=2,
                alpha = 0.75
            )
        df = annotation['compl']
        if len(df)>0:
            gc.show_ellipses(
                df['ra'],df['dec'],2.0*df['extent_semiminor'],2.0*df['extent_semimajor'],df['extent_angle'],
                edgecolor='grey',facecolor='None',linewidths=2,coords_frame='world',alpha=0.75
            )
        pos = annotation['positions']
        gc.show_markers(pos.ra,pos.dec,edgecolor='red',facecolor='None',linewidths=1,marker='P',s=100)
        plt.legend(handles=annotation['legend_hdls'],framealpha=0.5,loc='upper right')
    if datum['module'] == '__main__':
        title_bar = f"{datum['depth'].capitalize()} Image"
    else:
        rname = re.sub(r".*?\.residual\.(\w+?)\.fits$",r"\1",datum['fname'])
        title_bar = f"{datum['depth'].capitalize()} {datum['plt_labels'][rname]} Residual Image"
    title_bar = f"ID {datum['clump_id']:03d}: {title_bar}"
    gc.add_label(
        0.5,0.95,title_bar,relative=True,weight="bold",size="x-large",
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1', alpha=0.5)
    )
    plt.savefig(datum['pname'],transparent=True)
    plt.close(fig)
