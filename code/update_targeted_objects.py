import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.coordinates import SkyCoord
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import pdb


def update(semester,oldmask,masknum=2):

    #target data
    if masknum == 1:
        print('try again...')
    elif masknum == 2:
        if semester == '2017B':
            hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority.fits')
        else:
            hdulist = fits.open(semester+'/newfirm_mbs_deepz_sdss_priority.fits')
    elif masknum == 100:
        hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority_mask4.fits')
    else:
        if semester == '2017B':
            hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority_mask'+np.str(masknum-1)+'.fits')
        else:
            hdulist = fits.open(semester+'/newfirm_mbs_deepz_sdss_priority_mask'+np.str(masknum-1)+'.fits')
    
    nf = hdulist[1].data
    nf_id = nf['id']
    nf_pflag = nf['Pflag']
    nf_targ = nf['target']

    hdu1 = fits.open(semester+'/maskdesign/finishedmasks/'+oldmask+'.fits')
    m1 = hdu1[1].data
    m1_ID = m1['OBJECT']
    print(len(m1_ID))


    nfID_string = []
    for ii in range(len(nf_id)):
        nfID_string = np.append(nfID_string,np.str(nf_id[ii]))

    for i in range(len(m1_ID)):
        targcut = np.where(m1_ID[i] == nfID_string)[0]
        nf_targ[targcut] = 1

    print(len(nf_targ[np.where(nf_targ==1)[0]]))

    if semester == '2017B':
        hdulist.writeto(semester+'/newfirm_mbs_specz_sdss_priority_mask'+np.str(masknum)+'.fits',overwrite=True)
    else:
        hdulist.writeto(semester+'/newfirm_mbs_deepz_sdss_priority_mask'+np.str(masknum)+'.fits',overwrite=True)
