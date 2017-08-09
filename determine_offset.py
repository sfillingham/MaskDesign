from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

###########################################################################
##
## Determines the index and distance for each object in the guide catalog
## compared to the science catalog.
## The coordinates, (ra1,dec1), need to correspond to the "guide" catalog which
## requires an offset in order to align with the "science" catalog
## The coordinates, (ra2, dec2), are the "science" objects
##
############################################################################

def match(ra1,dec1,ra2,dec2):

    c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)

    idx, d2d, d3d = c.match_to_catalog_sky(catalog)

    return idx, d2d/u.degree

##############################################
##
## Plots the ra and dec offset for the science and guide catalogs
## If test == True, then plots will be generated; Else, no plots
## Threshold is the maximum distance two coords can be away from each other in
## order for a match to be valid, 'default' will use 1" radius
##
##############################################

def plot_offset(test, threshold):

    hdulist = fits.open('newfirm_mbs_deepz_sdss_catalog.fits')
    #hdulist2 = fits.open('newfirm_mbs_deepz_sdss_offset.fits')
    nf = hdulist[1].data
    sdss = hdulist[2].data

    #Make color cut below science cut to help with astrometric fitting
    #Make galaxy cut in sdss data for astrometric comparison
    colorcut = nf['R'] > 15.849#This is FLUX!!!! and a mag limit < 22
    sdsscut = sdss['R'] < 22#THIS IS Magnitude!!!

    #Clean ra and dec for astrometric matching
    nf_ra = nf['ra']#[colorcut]
    nf_dec = nf['dec']#[colorcut]

    sdss_ra = sdss['ra']#[sdsscut]
    sdss_dec = sdss['dec']#[sdsscut]

    print len(sdss_dec)
    print len(nf_dec)

    #use the match function to determine the closest objects
    idx, d2d = match(sdss_ra,sdss_dec,nf_ra,nf_dec)
    print len(idx)

    #print np.median(d2d)

    if threshold == 'default':

        limit = 0.0002778

    else:

        limit = 0.0002778*threshold

    hitcut = d2d < limit
    locations = idx[hitcut]

    ra_diff = np.empty(len(locations))
    dec_diff = np.empty(len(locations))

    for i in range(len(locations)):

        ra_diff[i] = nf_ra[locations[i]] - sdss_ra[hitcut][i]
        dec_diff[i] = nf_dec[locations[i]] - sdss_dec[hitcut][i]
    print len(ra_diff)

    if test == 'plot':
        plt.figure(1)
        n, bins, patches = plt.hist(ra_diff, 30, histtype = 'step', color = 'b')
        n, bins, patches = plt.hist(dec_diff, 30, histtype = 'step', color  = 'g')

        return np.median(ra_diff)*3600, np.median(dec_diff)*3600

    else:

        return np.median(ra_diff)*3600, np.median(dec_diff)*3600


def apply_offset():

    hdulist = fits.open('newfirm_mbs_deepz_sdss_catalog.fits')
    nf = hdulist[1].data
    sdss = hdulist[2].data
    
    ra_corr, dec_corr = plot_offset('plot',3)

    sdsscut = sdss['type_r'] == 6

    sdss_ra = sdss['ra']#[sdsscut]
    sdss_dec = sdss['dec']#[sdsscut]

    print ra_corr, dec_corr

    new_ra = sdss_ra + (ra_corr/3600)
    new_dec = sdss_dec + (dec_corr/3600)

    ########################################################
    #check the new ra and dec to ensure offset is < 0.001"
    ########################################################
    colorcut = nf['R'] > 15.849 # This is a flux limit corresponding to a mag limit < 22
    sdsscut = sdss['R'] < 22

    nf_ra = nf['ra'][colorcut]
    nf_dec = nf['dec'][colorcut]

    clean_ra = new_ra[sdsscut]
    clean_dec = new_dec[sdsscut]

    idx, d2d = match(clean_ra,clean_dec,nf_ra,nf_dec)

    limit = 0.0002778*3

    hitcut = d2d < limit
    locations = idx[hitcut]

    ra_diff = np.empty(len(locations))
    dec_diff = np.empty(len(locations))

    for i in range(len(locations)):

        ra_diff[i] = nf_ra[locations[i]] - clean_ra[hitcut][i]
        dec_diff[i] = nf_dec[locations[i]] - clean_dec[hitcut][i]
    

    if (np.median(ra_diff)*3600 < 0.1) & (np.median(dec_diff)*3600 < 0.1):

        plt.figure(2)
        n, bins, patches = plt.hist(ra_diff, 30, histtype = 'step', color = 'b')
        n, bins, patches = plt.hist(dec_diff, 30, histtype = 'step', color  = 'g')
        plt.show()

        nfcut = nf['K'] > 6.309
        new_nf = nf[nfcut]
        hdulist[1].data = new_nf
        
        sdss['ra'] = new_ra
        sdss['dec'] = new_dec

        hdulist.writeto('newfirm_mbs_deepz_sdss_offset.fits')

        print np.median(ra_diff)*3600, np.median(dec_diff)*3600

    else:

        print 'Check routine, offset sucks!!!'

    

    

    

    

    
    
    
