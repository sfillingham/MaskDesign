# mask-design
Useful scripts for multi object spectroscopy mask design and general
spectroscopic observation planning.
Specifically, these scripts were used to design Keck/DEIMOS slitmasks via the dsimulator package.

## Keck/DEIMOS
- determine priority.py
- make_dsim_input.py
- make_dsim_input_mask2.py
- make_dsim_input_mask3.py
- dsim_output_check.py

## General Observing Tools
- catalog_combing.py
- determine_offset.py
- master_catalog.py
- update_targets.py
- make_regionfile.py

## Overview of the Step-by-step Process
This describes the routine and order by which I compiled the object
catalog for the DEIMOS runs related to the z ~ 1 study of satellite
quenching in group environments. 

Begin by downloading all the necessary catalogs of data
- Photometric science catalog (e.g. NEWFIRM_MBS, 3D-HST, etc.)
- Spectroscopic catalog for any objects in the same field as the photometric catalog(e.g. DEEP2)
- SDSS photometry over the same field for both alignment and guide stars

Combine all the catalogs into one master FITS file (use 'newfirm_mbs' as working dir)
- Run 'master_catalog.get_catalog()' to generate the NEWFIRM catalog
- Run 'catalog_combine.specmatch()' to add DEEP spectra to the catalog
- Run 'catalog_combine.sdss_match()' to add SDSS alignment, guide, and
  astrometry objects in the NEWFIRM field

Determine completeness levels in each photometric band at desired redshift
- Plot each band vs stellar mass in a given redshift window; verifying the magnitude zeropoint from literature
- From survey literature, determine where the completeness limit is in survey selection band
- Verify completeness is consistent with science goals
- Determine magnitude cut in selection band

Determine/Apply the relative astrometric offset between the science targets and the guide/alignment stars
- Remove any faint objects below completeness levels, this helps ensure accurate astrometric matching
- Determine and apply the astrometric offsets to the guide stars,
    alignment stars, and test objects in SDSS 'determine_offset.apply_offset()'
- At this point all the science, alignment, and guide stars should have the same relative astrometry
- Verify the offsets by plotting the ra/dec offset for science and test objects

Create final observing catalog
- Make completeness cut in science band, and propogate that cut into observing bands (i.e. selection cut in K, observing in R)
- Establish a priority method and run the following scripts
        'determine_priority.pass1()'
    	'determine_priority.pass2()'
- Create a DSIM input file using 'make_dsim_input.input()'
- Pass this file into DSIMULATOR and generate an output FITS file

Rinse and repeat process as needed until satisfied.

