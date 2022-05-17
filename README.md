# Needed data files for running and testing hydra
Need to download the following file data.tar from  https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/cirada/catalogues/github_files/hydra_data.tar

Untar it the directory 'hydra'.
## Ubuntu Installation:

### Install Docker:
https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository

### Add user privileges:
https://docs.docker.com/engine/install/linux-postinstall/

### Add Docker compose:
https://docs.docker.com/compose/install/#install-compose

### Building Containers:
```
$ cd modules/config
$ docker-compose build
```

### Restarting Docker:
https://docs.docker.com/config/daemon/systemd/
```
sudo systemctl restart docker
```

---

Code Count:
```
$ pygount --suffix=py,pyx,py-jn2,dcr,R,yml --format=summary --folders-to-skip=data,.git,.depricated .

 Language    Files    %     Code    %     Comment    %   
-----------  -----  ------  ----  ------  -------  ------
Python          41   70.69  9534   83.30     1944   89.67
YAML             7   12.07  1126    9.84       35    1.61
S                1    1.72   450    3.93       76    3.51
Cython           2    3.45   336    2.94      113    5.21
__unknown__      7   12.07     0    0.00        0    0.00
-----------  -----  ------  ----  ------  -------  ------
Sum total       58          11446             2168
$ date
Fri 13 Aug 2021 12:40:37 CDT
$
```
## Commands:
The following diagram shows the Hydra sofware workflow.

<img src="docs/pics/EMU_Pipeline.png" width=500px>

The Hydra software suite consists of the following commands.
<table>
<tr><th>Tool</th><th>Command</th><th>Function</th></tr>
<tr><th>Homados</th><td>homados.py</td><td>Provides functionality for creating shallow images with n&sigma; noise (default, n=5), inverting images, and computing RMS, &mu;, <i>etc.</i> statistics.</td></tr>
<tr><th>Cerberus</th><td>cerberus.py</td><td>Alows the user to run Hydra's source finders collectively or individualy.</td></tr>
<tr><th>Typhon</th><td>typhon.py</td><td>Optimizes source-finder RMS and Island parameters to a standard baseline that can be used for comparison purposes. This tool uses Homados and Cerberus.</td></tr>
<tr><th>Hydra</th><td>hydra.py</td><td>The main tool which glues everthing together, by running Typhon for deep and shallow images, to produce various data products.</td></tr>
</table>

All commands have builtin help. The Homados, Cerberus, and Typhon command-line interfaces are ancillary tools used for diagnostic purposes. These tools also have library interfaces and are pipelined together, so as to be utilized by Hydra. Hydra is the main tool that glues everything together, by running Typhon for deep and shallow images, and producing the following data products:
<ul>
    <li> Typhon Metrics
    <ul>
        <li> Deep/Shallow Diagnostic Plots of
        <ul>
            <li> PRD
            <li> PRD CPU Times
            <li> Residual RMS
            <li> Residual MADFM
            <li> Residual &Sigma;I<sup>2</sup>
        </ul>
        <li> Table of deep and shallow optimized RMS and island parameters
    </ul>
    <li> Deep/Shallow Catalogues
    <ul>
        <li> Source finder Catalogues
        <li> Cluster Catalogue
        <li> Clump Catalogue
    </ul>
    <li> Optional Simulated Input Catalogue
    <li> Deep/Shallow Cutouts
    <ul>
        <li> Un/annotated Images
        <li> Un/annotated Residual Images
    </ul>
    <li> Diagnostic Plots
    <ul>
        <li> Clump Size Distributions
        <li> Detections <i>vs.</i> S/N
        <li> Completeness <i>vs.</i> S/N
        <li> Reliability <i>vs.</i> S/N
        <li> S<sub>out</sub>/S<sub>in</sub> <i>vs.</i> S/N
        <li> False-Positive <i>vs.</i> S/N
    </ul>
    <li> wrt injected and deep sources
    <li> Local Web-browser Tool
</ul>
All of this information is stored in a tarball. With the exception of the source finder, simulated input, and clump catalogues, the local web-browser tool allows the user to explore all of these data products, and is accessible through an ``index.html`` file in the main ``tar`` directory.

The following is the help output of the Hydra tool.
```
$ python hydra.py --help
Usage: hydra.py [OPTIONS] FITS_FILE

  Performs deep-shallow analaysis.

Options:
  --catalogues-yml TEXT  YML file pointing to catalogoues to match.
  --bypass-archived      Use FITS_FILE_PREFIX.hydra.tar.gz archive.
  --use TEXT             Use FITS_FILE_PREFIX.typhon.tar.gz optimization.
  --help                 Show this message and exit.
$
```
The ``--catalogues-yml`` option is for referencing a &lt;catalogues&gt;.yml configuration file, which points to catalogues to merge with typhon's catalogue output. This is typically used for injecting simulated sources from a simulated image for completeness and reliability studies, <i>etc.</i> For example, the following Hydra command will process (in the background) a simulated 2x2 simulated FITS image (``emu_simulated_2x2.fits``) and collate the results with the corresponding simulated catalogue (``art_src_list_2x2deg.fits``), defined in the ``art_src_list_2x2deg.yml`` configuration file (logging the results in ``sim_2x2.log``).
```
$ nohup python hydra.py data/simulated/2x2/emu_simulated_2x2.fits --catalogues-yml data/simulated/2x2/art_src_list_2x2deg.yml > sim_2x2.log 2>&1 & 
```
This will generate tarballs from Typhon, for deep (``emu_simulated_2x2.fits``) and shallow (``emu_simulated_2x2.shallow.fits``) images, and Hydra, which contains final results: <i>i.e.</i>,
```
$ ls -1 data/simulated/2x2/*.tar.gz
data/simulated/2x2/emu_simulated_2x2.hydra.tar.gz
data/simulated/2x2/emu_simulated_2x2.shallow.typhon.tar.gz
data/simulated/2x2/emu_simulated_2x2.typhon.tar.gz
$
``` 
Untarring ``emu_simulated_2x2.hydra.tar.gz`` (<i>a la</i>, ``tar xf emu_simulated_2x2.hydra.tar.gz``), reveals an ``index.html`` in the top directory of ``emu_simulated_2x2.hydra_dir`` (created after untarring), from which one can navigate to launch the Hydra Viewer.

Running the Hydra command will simply load the ``emu_simulated_2x2.hydra.tar.gz`` again, if it already exists; otherwiser, it will rebuild it from the Typhon tarballs (failing that, it'll redo the whole run). The ``--bypass-archived`` forces Hydra to reprocesses the Typhon tarballs. This is intended for software development purposes, and future development.

Finally, the ``--use`` option is for running Hydra on an image using the Typhon optimiztion results of a reference image. For example, the following command processing the VLASS image ``vlass_qle1_t12t18_2x2deg_sample.fit`` using the optimization results of the ``emu_simulated_2x2.fits`` simulated (which will use using its Typhon tarball results, or generated them).
```
nohup python hydra.py data/real/vlass/vlass_qle1_t12t18_2x2deg_sample.fits --use data/simulated/2x2/emu_simulated_2x2.fits > vlass.log 2>&1 &
```

All of these examples can be run using the testbed data in the ``data/`` directory at the top level of this repo. For the purposes of rapid software development and testing (<i>i.e.</i>, not quality control) it is recommend one uses one of the smaller test images: <i>e.g.</i>,
```
$ python hydra.py data/real/vlass/J165209+212528_s3arcmin_VLASS.fits
```
----
### Citations

Please include the following material in acknowlegements to research papers:

Hydra \citep{boyce_2022a,boyce_2022b} is written primarily in Python, with some elements of Cython \citep{behnel_2011} and R, along with their standard libraries. Hydra uses alphashape \citep{bellock_2021}, APLpy \citep{robitaille_2012}, Astropy \citep{astropy_2013, astropy_2018}, Matplotlib \citep{hunter_2007}, NumPy \citep{harris_2020}, and pandas \citep{mckinney_2010,reback_2020} Python libraries commonly used in astronomy. Hydra utilizes click, gzip, Jinja, tarfile, and YAML Python libraries as part its overall architectural infrastructure. Hydra encapsulates Aegean \citep{hancock_2012,hancock_2018}, Caesar \citep{riggi_2016,riggi_2019}, ProFound \citep{robotham_2018,hale_2019}, PyBDSF \citep{mohan_2015}, and Selavy \citep{whiting_2012} source finder software using Docker. We acknowledge the authors of the aforementioned software applications, languages, and libraries. Hydra was created as a joint CIRADA -- ASKAP intitiative through CFI funding (Project 35999).

with BibTeX items,

\@article{astropy_2013,
Adsnote = {Provided by the SAO/NASA Astrophysics Data System},
Adsurl = {http://adsabs.harvard.edu/abs/2013A%26A...558A..33A},
Archiveprefix = {arXiv},
Author = {{Astropy Collaboration} and {Robitaille}, T.~P. and {Tollerud}, E.~J. and {Greenfield}, P. and {Droettboom}, M. and {Bray}, E. and {Aldcroft}, T. and {Davis}, M. and {Ginsburg}, A. and {Price-Whelan}, A.~M. and {Kerzendorf}, W.~E. and {Conley}, A. and {Crighton}, N. and {Barbary}, K. and {Muna}, D. and {Ferguson}, H. and {Grollier}, F. and {Parikh}, M.~M. and {Nair}, P.~H. and {Unther}, H.~M. and {Deil}, C. and {Woillez}, J. and {Conseil}, S. and {Kramer}, R. and {Turner}, J.~E.~H. and {Singer}, L. and {Fox}, R. and {Weaver}, B.~A. and {Zabalza}, V. and {Edwards}, Z.~I. and {Azalee Bostroem}, K. and {Burke}, D.~J. and {Casey}, A.~R. and {Crawford}, S.~M. and {Dencheva}, N. and {Ely}, J. and {Jenness}, T. and {Labrie}, K. and {Lim}, P.~L. and {Pierfederici}, F. and {Pontzen}, A. and {Ptak}, A. and {Refsdal}, B. and {Servillat}, M. and {Streicher}, O.},
Doi = {10.1051/0004-6361/201322068},
Eid = {A33},
Eprint = {1307.6212},
Journal = {\aap},
Keywords = {methods: data analysis, methods: miscellaneous, virtual observatory tools},
Month = oct,
Pages = {A33},
Primaryclass = {astro-ph.IM},
Title = {{Astropy: A community Python package for astronomy}},
Volume = 558,
Year = 2013,
Bdsk-Url-1 = {https://dx.doi.org/10.1051/0004-6361/201322068}}

\@ARTICLE{astropy_2018,
       author = {{Astropy Collaboration} and {Price-Whelan}, A.~M. and
         {Sip{\H{o}}cz}, B.~M. and {G{\"u}nther}, H.~M. and {Lim}, P.~L. and
         {Crawford}, S.~M. and {Conseil}, S. and {Shupe}, D.~L. and
         {Craig}, M.~W. and {Dencheva}, N. and {Ginsburg}, A. and {Vand
        erPlas}, J.~T. and {Bradley}, L.~D. and {P{\'e}rez-Su{\'a}rez}, D. and
         {de Val-Borro}, M. and {Aldcroft}, T.~L. and {Cruz}, K.~L. and
         {Robitaille}, T.~P. and {Tollerud}, E.~J. and {Ardelean}, C. and
         {Babej}, T. and {Bach}, Y.~P. and {Bachetti}, M. and {Bakanov}, A.~V. and
         {Bamford}, S.~P. and {Barentsen}, G. and {Barmby}, P. and
         {Baumbach}, A. and {Berry}, K.~L. and {Biscani}, F. and {Boquien}, M. and
         {Bostroem}, K.~A. and {Bouma}, L.~G. and {Brammer}, G.~B. and
         {Bray}, E.~M. and {Breytenbach}, H. and {Buddelmeijer}, H. and
         {Burke}, D.~J. and {Calderone}, G. and {Cano Rodr{\'\i}guez}, J.~L. and
         {Cara}, M. and {Cardoso}, J.~V.~M. and {Cheedella}, S. and {Copin}, Y. and
         {Corrales}, L. and {Crichton}, D. and {D'Avella}, D. and {Deil}, C. and
         {Depagne}, {\'E}. and {Dietrich}, J.~P. and {Donath}, A. and
         {Droettboom}, M. and {Earl}, N. and {Erben}, T. and {Fabbro}, S. and
         {Ferreira}, L.~A. and {Finethy}, T. and {Fox}, R.~T. and
         {Garrison}, L.~H. and {Gibbons}, S.~L.~J. and {Goldstein}, D.~A. and
         {Gommers}, R. and {Greco}, J.~P. and {Greenfield}, P. and
         {Groener}, A.~M. and {Grollier}, F. and {Hagen}, A. and {Hirst}, P. and
         {Homeier}, D. and {Horton}, A.~J. and {Hosseinzadeh}, G. and {Hu}, L. and
         {Hunkeler}, J.~S. and {Ivezi{\'c}}, {\v{Z}}. and {Jain}, A. and
         {Jenness}, T. and {Kanarek}, G. and {Kendrew}, S. and {Kern}, N.~S. and
         {Kerzendorf}, W.~E. and {Khvalko}, A. and {King}, J. and {Kirkby}, D. and
         {Kulkarni}, A.~M. and {Kumar}, A. and {Lee}, A. and {Lenz}, D. and
         {Littlefair}, S.~P. and {Ma}, Z. and {Macleod}, D.~M. and
         {Mastropietro}, M. and {McCully}, C. and {Montagnac}, S. and
         {Morris}, B.~M. and {Mueller}, M. and {Mumford}, S.~J. and {Muna}, D. and
         {Murphy}, N.~A. and {Nelson}, S. and {Nguyen}, G.~H. and
         {Ninan}, J.~P. and {N{\"o}the}, M. and {Ogaz}, S. and {Oh}, S. and
         {Parejko}, J.~K. and {Parley}, N. and {Pascual}, S. and {Patil}, R. and
         {Patil}, A.~A. and {Plunkett}, A.~L. and {Prochaska}, J.~X. and
         {Rastogi}, T. and {Reddy Janga}, V. and {Sabater}, J. and
         {Sakurikar}, P. and {Seifert}, M. and {Sherbert}, L.~E. and
         {Sherwood-Taylor}, H. and {Shih}, A.~Y. and {Sick}, J. and
         {Silbiger}, M.~T. and {Singanamalla}, S. and {Singer}, L.~P. and
         {Sladen}, P.~H. and {Sooley}, K.~A. and {Sornarajah}, S. and
         {Streicher}, O. and {Teuben}, P. and {Thomas}, S.~W. and
         {Tremblay}, G.~R. and {Turner}, J.~E.~H. and {Terr{\'o}n}, V. and
         {van Kerkwijk}, M.~H. and {de la Vega}, A. and {Watkins}, L.~L. and
         {Weaver}, B.~A. and {Whitmore}, J.~B. and {Woillez}, J. and
         {Zabalza}, V. and {Astropy Contributors}},
        title = "{The Astropy Project: Building an Open-science Project and Status of the v2.0 Core Package}",
      journal = {\aj},
     keywords = {methods: data analysis, methods: miscellaneous, methods: statistical, reference systems, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2018,
        month = sep,
       volume = {156},
       number = {3},
          eid = {123},
        pages = {123},
          doi = {10.3847/1538-3881/aabc4f},
archivePrefix = {arXiv},
       eprint = {1801.02634},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}}



\@misc{boyce_2022a,
  AUTHOR = {Boyce, M. M. and Hopkins, A. M. and Riggi, S. and Rudnick, L. and Ramsay, M. and Hale, C. L. and  Marvil, J. and Whiting, M. and Venkataraman, P. and O'Dea, C. P. and Baum, S. A. and Gordon, Y. A. and Vantyghem, A. N. and Dionyssiou, M. and Andernach, H. and Collier, J. D. and English, J. and Koribalski, B. S. and Leahy, D. and Michałowski, M. J. and Safi-Harb, S. and Vaccari, M..},
  TITLE = {Hydra I: An extensible multi-source-finder comparison and cataloguing tool -- In preparation (PASA).},
  YEAR  = {2022}}

\@misc{boyce_2022b,
  AUTHOR = {Boyce, M. M. and Hopkins, A. M. and Riggi, S. and Rudnick, L. and Ramsay, M. and Hale, C. L. and  Marvil, J. and Whiting, M. and Venkataraman, P. and O'Dea, C. P. and Baum, S. A. and Gordon, Y. A. and Vantyghem, A. N. and Dionyssiou, M. and Andernach, H. and Collier, J. D. and English, J. and Koribalski, B. S. and Leahy, D. and Michałowski, M. J. and Safi-Harb, S. and Vaccari, M..},
  TITLE = {Hydra II: Characterization of Aegean, Caesar, ProFound, PyBDSF, and Selavy -- In preparation (PASA).},
  YEAR  = {2022}}

\@article{behnel_2011,
  title={Cython: The best of both worlds},
  author={Behnel, Stefan and Bradshaw, Robert and Citro, Craig and Dalcin, Lisandro and Seljebotn, Dag Sverre and Smith, Kurt},
  journal={Computing in Science \& Engineering},
  volume={13},
  number={2},
  pages={31--39},
  year={2011},
  publisher={IEEE}}

\@article{bellock_2021,
    author       = {Bellock, Kenneth E.},
    journal      = {PyPI},
    year         = 2021,
    version      = {latest},
    doi          = {https://pypi.org/project/alphashape/}}

\@article{hale_2019,
  title={Radio source extraction with ProFound},
  author={{Hale}, C.~L. and {Robotham}, A.~S.~G. and {Davies}, L.~J.~M. and {Jarvis}, M.~J. and {Driver}, S.~P. and {Heywood}, I.},
  journal={MNRAS},
  volume={487},
  year={2019},
  doi={https://doi.org/10.1093/mnras/stw982}}

\@article{hancock_2012,
  title={Compact continuum source finding for next generation radio surveys},
  author={{Hancock}, P.~J. and {Murphy}, T. and {Gaensler}, B.~M. and {Hopkins}, A. and {Curran}, J.~R.},
  journal={MNRAS},
  volume={422},
  year={2012},
  doi={https://doi.org/10.1111/j.1365-2966.2012.20768.x}}

\@article{hancock_2018,
  title={Source Finding in the Era of the SKA (Precursors): AEGEAN 2.0},
  author={{Hancock}, P.~J. and {Cathryn}, M.~T. and {Hurley-Walker}, N.},
  journal={PASA},
  volume={35},
  year={2018},
  doi={https://doi.org/10.1017/pasa.2018.3}}

\@Article{harris_2020,
 title         = {Array programming with {NumPy}},
 author        = {Charles R. Harris and K. Jarrod Millman and St{\'{e}}fan J.
                 van der Walt and Ralf Gommers and Pauli Virtanen and David
                 Cournapeau and Eric Wieser and Julian Taylor and Sebastian
                 Berg and Nathaniel J. Smith and Robert Kern and Matti Picus
                 and Stephan Hoyer and Marten H. van Kerkwijk and Matthew
                 Brett and Allan Haldane and Jaime Fern{\'{a}}ndez del
                 R{\'{i}}o and Mark Wiebe and Pearu Peterson and Pierre
                 G{\'{e}}rard-Marchant and Kevin Sheppard and Tyler Reddy and
                 Warren Weckesser and Hameer Abbasi and Christoph Gohlke and
                 Travis E. Oliphant},
 year          = {2020},
 month         = sep,
 journal       = {Nature},
 volume        = {585},
 number        = {7825},
 pages         = {357--362},
 doi           = {10.1038/s41586-020-2649-2},
 publisher     = {Springer Science and Business Media {LLC}},
 url           = {https://doi.org/10.1038/s41586-020-2649-2}}

\@Article{hunter_2007,
  Author    = {Hunter, J. D.},
  Title     = {Matplotlib: A 2D graphics environment},
  Journal   = {Computing in Science \& Engineering},
  Volume    = {9},
  Number    = {3},
  Pages     = {90--95},
  abstract  = {Matplotlib is a 2D graphics package used for Python for
  application development, interactive scripting, and publication-quality
  image generation across user interfaces and operating systems.},
  publisher = {IEEE COMPUTER SOC},
  doi       = {10.1109/MCSE.2007.55},
  year      = 2007}

\@InProceedings{mckinney_2010,
  author    = {{W}es {M}c{K}inney},
  title     = {{D}ata {S}tructures for {S}tatistical {C}omputing in {P}ython},
  booktitle = {{P}roceedings of the 9th {P}ython in {S}cience {C}onference},
  pages     = {56 - 61},
  year      = {2010},
  editor    = {{S}t\'efan van der {W}alt and {J}arrod {M}illman},
  doi       = {10.25080/Majora-92bf1922-00a}}

\@article{mohan_2015,
    author={Mohan, N. and Rafferty, D.},
    title={PyBDSF: Python Blob Detection and Source Finder},
    journal={Astrophysics Source Code Library},
    volume={ascl:1107.013},
    year={2015},
    doi={https://ui.adsabs.harvard.edu/abs/2015ascl.soft02007M/abstract}}

\@article{reback_2020,
    author       = {{Pandas Development Team}},
    journal      = {pandas-dev/pandas: Pandas},
    month        = feb,
    year         = 2020,
    publisher    = {Zenodo},
    version      = {latest},
    doi          = {https://doi.org/10.5281/zenodo.3509134}}

\@article{riggi_2016,
  title={Automated detection of extended sources in radio maps: progress from the SCORPIO survey},
  author={{Riggi}, S. and {Ingallinera}, A. and {Leto}, P. and {Cavallaro}, F. and {Bufano}, F. and {Schillir\'{o}}, F. and {Trigilio}, C. and {Umana}, G. and {Buemi}, C.~S. and {Norris}, R.~P.},
  journal={MNRAS},
  volume={460},
  year={2016},
  doi={https://doi.org/10.1093/mnras/stw982}}

\@article{riggi_2019,
  title={CAESAR source finder: Recent developments and testing},
  author={{Riggi}, S. and {Vitello}, F. and {Becciani}, U. and {Buemi}, F. and {Calanducci}, A. and {Cavallaro}, F. and {Costa}, A. and {Ingallinera}, A. and {Leto}, P. and {Norris}, R.~P. and {Schillir\`{o}}, F. and {Sciacca}, E. and {Trigilio}, C. and {Umana}, G.},
  journal={PASA},
  volume={36},
  year={2019},
  doi={https://doi.org/10.1017/pasa.2019.29}}

\@article{robitaille_2012,
  title={Search and detection of low frequency radio transients},
  author={Robitaille, Thomas and Bressert, Eli},
  journal={Astrophysics Source Code Library, record ascl:1208.017},
  year={2012},
  doi={https://ui.adsabs.harvard.edu/abs/2012ascl.soft08017R/abstract}}

\@article{robotham_2018,
  title={ProFound: Source Extraction and Application to Modern Survey Data},
  author={{Robotham}, A.~S.~G. and {Davies}, L.~J.~M. and {Driver}, S.~P. and {Koushan}, S.~P. and {Taranu}, D.~S. and {Casura}, S. and {Liske}, J.},
  journal={MNRAS},
  volume={476},
  year={2018},
  doi={https://doi.org/10.1093/mnras/sty440}}

\@article{whiting_2012,
  title={Source-Finding for the Australian Square Kilometre Array Pathfinder},
  author={Whiting, M. and Humphreys, B.},
  journal={PASA},
  volume={29},
  pages={371-381},
  year={2012},
  doi = {https://doi.org/10.1071/AS12028}}


