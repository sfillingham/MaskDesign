from astropy import wcs
import astropy.io.fits as fits
import numpy as np
import skimage.draw as draw
from astropy.table import Table


def region():

    imagehdu = fits.open('/Users/Sean/research/egs/AEGIS-N2_K_sci.fits')
    m1hdu = fits.open('masks_round2/mask10c.fits')
    m2hdu = fits.open('masks_round2/mask10c.fits')
    #m30hdu = fits.open('masks_round2/mask30c.fits')
    hdulist = fits.open('newfirm_mbs_deepz_sdss_priority_pass1.fits')
    groupdata = Table.read('/Users/Sean/research/deep/groupcat_deimos16a.txt', format = 'ascii')
    
    w = wcs.WCS(imagehdu[0].header)

    group_ra = np.array(groupdata['RA'])
    group_dec = np.array(groupdata['DEC'])
    group_priority = np.array(groupdata['PRIORITY'])

    m1_data = m1hdu[1].data
    m2_data = m2hdu[1].data
    #m30_data = m30hdu[1].data

    m1_ra = m1_data['RA_OBJ']
    m2_ra = m2_data['RA_OBJ']
    #m30_ra = m30_data['RA_OBJ']
    m1_dec = m1_data['DEC_OBJ']
    m2_dec = m2_data['DEC_OBJ']
    #m30_dec = m30_data['DEC_OBJ']
    m1_pa = m1_data['MajAxPA']
    m2_pa = m2_data['MajAxPA']

    nf_data = hdulist[1].data
    nf_ra = nf_data['ra']
    nf_dec = nf_data['dec']
    nf_pflag = nf_data['Pflag']
    print np.str(np.max(nf_pflag))
    nf_deep = nf_data['deepz']
    nf_ellip = nf_data['K_ellip']
    nf_theta = nf_data['K_theta_J2000']

    nf_theta[np.where(nf_ellip < 0.3)[0]] = 0
    print nf_ellip[:20]
    print nf_theta[:20]
    
    #nf_rflux = nf_data['R']
    #nf_rmag = 25.-2.5*np.log10(nf_rflux)
    bright = (nf_pflag == 10000)#(nf_rmag <= 24) & (nf_rmag > 0)
    dull = (nf_pflag == 100)#(nf_rmag > 24) & (nf_rmag < 25.5)
    
    sdss_data = hdulist[2].data
    sdss_ra = sdss_data['ra']
    sdss_dec = sdss_data['dec']
    sdss_pflag = sdss_data['Pflag']

    guidecut = sdss_pflag == -1
    aligncut = sdss_pflag == -2
    deeponlycut = (nf_deep >= 3) & (nf_pflag == 1)
    deepcut = (nf_deep >= 3) 
    scicut1 = nf_pflag == 1
    scicut10 = nf_pflag == 10
    lowcut = nf_pflag == 1000
    highcut = nf_pflag == 10000
    membercut = (nf_pflag == 10000) ^ (nf_pflag == 100)
    
    Ghigh = group_priority == 1
    Glow = group_priority == 2
    
    
    sdssxx,sdssyy = w.wcs_world2pix(sdss_ra,sdss_dec,0,ra_dec_order=True)
    nfxx,nfyy = w.wcs_world2pix(nf_ra,nf_dec,0,ra_dec_order=True)
    m1xx,m1yy = w.wcs_world2pix(m1_ra,m1_dec,0,ra_dec_order=True)
    m2xx,m2yy = w.wcs_world2pix(m2_ra,m2_dec,0,ra_dec_order=True)
    gxx,gyy = w.wcs_world2pix(group_ra,group_dec,0,ra_dec_order=True)


    #f = open('2016a_NMBS_DEEP_region.reg','w')
    #f = open('2016a_NMBS_region_pass1.reg','w')
    #f = open('2016a_NMBS_region_pass2.reg','w')
    #f = open('2016a_NMBS_region_stars.reg','w')
    #f = open('2016a_NMBS_region_mask1c.reg','w')
    #f = open('2016a_NMBS_region_mask2c.reg','w')
    #f = open('2016a_NMBS_region_mask3.reg','w')
    #f = open('2016a_NMBS_region_mask4.reg','w')
    #f = open('2016a_NMBS_region_mask10.reg','w')
    f = open('2016a_NMBS_region_groupmembers.reg','w')
    #f = open('2016a_NMBS_region_PAcheck.reg','w')

    hxx = nfxx[highcut]
    hyy = nfyy[highcut]
    htheta = nf_theta[highcut]
    deephxx = nfxx[membercut] #& deepcut]
    deephyy = nfyy[membercut] #& deepcut]

    lxx = nfxx[lowcut]
    lyy = nfyy[lowcut]
    ltheta = nf_theta[lowcut]
    deeplxx = nfxx[lowcut & deepcut]
    deeplyy = nfyy[lowcut & deepcut]

    ghxx = gxx[Ghigh]
    ghyy = gyy[Ghigh]

    glxx = gxx[Glow]
    glyy = gyy[Glow]

    sci1xx = nfxx[scicut1]
    sci1yy = nfyy[scicut1]
    sci10xx = nfxx[scicut10]
    sci10yy = nfyy[scicut10]

    deepxx = nfxx[deeponlycut]
    deepyy = nfyy[deeponlycut]

    guidexx = sdssxx[guidecut]
    guideyy = sdssyy[guidecut]
    alignxx = sdssxx[aligncut]
    alignyy = sdssyy[aligncut]

    print len(hxx)
    #print len(ghxx)

    print len(lxx)
    #print len(glxx)
    
    #for i in range(len(hxx)):
        #print >>f, 'circle '+np.str(hxx[i])+' '+np.str(hyy[i])+' 15 #color=magenta'

    #for i in range(len(m1xx)):
        #print >>f, 'circle '+np.str(m1xx[i])+' '+np.str(m1yy[i])+' 15 #color=yellow'

    #for i in range(len(m2xx)):
        #print >>f, 'circle '+np.str(m2xx[i])+' '+np.str(m2yy[i])+' 20 #color=red width=2'

    for i in range(len(deephxx)):
        print >>f, 'box '+np.str(deephxx[i])+' '+np.str(deephyy[i])+' 25 25 45#color=magenta width=2'

    #for i in range(len(lxx)):
        #print >>f, 'circle '+np.str(lxx[i])+' '+np.str(lyy[i])+' 10 #color=cyan'

    #for i in range(len(deeplxx)):
        #print >>f, 'box '+np.str(deeplxx[i])+' '+np.str(deeplyy[i])+' 25 25 45#color=cyan width=2'

    for i in range(len(ghxx)):
        print >>f, 'circle '+np.str(ghxx[i])+' '+np.str(ghyy[i])+' 25 #color=red width=2'

    for i in range(len(glxx)):
        print >>f, 'circle '+np.str(glxx[i])+' '+np.str(glyy[i])+' 25 #color=blue width=2'

    #for i in range(len(deepxx)):
        #print >>f, 'box '+np.str(deepxx[i])+' '+np.str(deepyy[i])+' 25 25 45#color=green width=2'

    #for i in range(len(guidexx)):
        #print >>f, 'box '+np.str(guidexx[i])+' '+np.str(guideyy[i])+' 25 25 #color=blue'

    #for i in range(len(alignxx)):
        #print >>f, 'box '+np.str(alignxx[i])+' '+np.str(alignyy[i])+' 20 20 45 #color=red'

    #for i in range(len(nfxx[bright])):
        #print >>f, 'box '+np.str(nfxx[bright][i])+' '+np.str(nfyy[bright][i])+' 25 25 0 #color=red'

    #for i in range(len(nfxx[dull])):
        #print >>f, 'box '+np.str(nfxx[dull][i])+' '+np.str(nfyy[dull][i])+' 25 25 0 #color=blue'

    #for i in range(len(m1xx)):
        #print >>f, 'box '+np.str(m1xx[i])+' '+np.str(m1yy[i])+' 20 4 '+np.str(m1_pa[i]+90.)+' #color=yellow'
        #print >>f, 'box '+np.str(m1xx[i])+' '+np.str(m1yy[i])+' 20 4 20 #color=yellow'

    #for i in range(len(m2xx)):
        #print >>f, 'box '+np.str(m2xx[i])+' '+np.str(m2yy[i])+' 20 4 '+np.str(m2_pa[i]+90.)+' #color=red'
        #print >>f, 'box '+np.str(m2xx[i])+' '+np.str(m2yy[i])+' 20 4 20 #color=yellow'

    #for i in range(len(lxx)):
        #print >>f, 'box '+np.str(lxx[i])+' '+np.str(lyy[i])+' 20 4 '+np.str(ltheta[i])+' #color=green'

    f.close()

    return nf_ellip[lowcut], nf_theta[lowcut]

        
