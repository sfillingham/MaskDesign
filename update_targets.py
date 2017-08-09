import numpy as np
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)

def update1(mask1, mask2, originalfile, passnum):

    hdulist = fits.open(originalfile+'.fits')
    targetdata = hdulist[1].data
    stardata = hdulist[2].data

    targetID = targetdata['id']
    targetINFO = targetdata['target']
    print len(targetINFO[np.where(targetINFO == 1)[0]])

    hdum1 = fits.open(mask1+'.fits')
    hdum2 = fits.open(mask2+'.fits')

    m1data = hdum1[1].data
    m1ID = m1data['OBJECT']
    m2data = hdum2[1].data
    m2ID = m2data['OBJECT']

    #######################################################
    ## Change the 'target' of objects already observed on Mask 1
    #######################################################
    targetIDtest = np.array([])

    for i in range(len(targetID)):
        targetIDtest = np.append(targetIDtest, np.str(targetID[i]))

    for i in range(len(m1ID)):
        targetINFO[np.where(targetIDtest == m1ID[i])[0]] = 1

    for i in range(len(m2ID)):
        targetINFO[np.where(targetIDtest == m2ID[i])[0]] = 1

    if np.max(targetINFO) != 1:
        print 'Holy FUCK'

    else:
        print 'Good to go'
        print len(targetINFO[np.where(targetINFO == 1)[0]])

    print 'saving...'
    hdulist.writeto('newfirm_mbs_deepz_sdss_priority_pass'+np.str(passnum)+'.fits')

    #return hdulist[1].data

def update2(mask1, mask2, mask3, originalfile, passnum):

    hdulist = fits.open(originalfile+'.fits')
    targetdata = hdulist[1].data
    stardata = hdulist[2].data

    targetID = targetdata['id']
    targetINFO = targetdata['target']
    print len(targetINFO[np.where(targetINFO == 1)[0]])

    hdum1 = fits.open(mask1+'.fits')
    hdum2 = fits.open(mask2+'.fits')
    hdum3 = fits.open(mask3+'.fits')

    m1data = hdum1[1].data
    m1ID = m1data['OBJECT']
    m2data = hdum2[1].data
    m2ID = m2data['OBJECT']
    m3data = hdum3[1].data
    m3ID = m3data['OBJECT']

    #######################################################
    ## Change the 'target' of objects already observed on Mask 1
    #######################################################
    targetIDtest = np.array([])

    for i in range(len(targetID)):
        targetIDtest = np.append(targetIDtest, np.str(targetID[i]))

    for i in range(len(m1ID)):
        targetINFO[np.where(targetIDtest == m1ID[i])[0]] = 1

    for i in range(len(m2ID)):
        targetINFO[np.where(targetIDtest == m2ID[i])[0]] = 1

    for i in range(len(m3ID)):
        targetINFO[np.where(targetIDtest == m3ID[i])[0]] = 1

    if np.max(targetINFO) != 1:
        print 'Holy FUCK'

    else:
        print 'Good to go'
        print len(targetINFO[np.where(targetINFO == 1)[0]])

    print 'saving...'
    hdulist.writeto('newfirm_mbs_deepz_sdss_priority_pass'+np.str(passnum)+'.fits')


    
