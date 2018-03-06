from astropy import wcs
import astropy.io.fits as fits
import numpy as np
import skimage.draw as draw
from astropy.table import Table


def region(field,semester,maskname):

    if field=='egs':
        imagehdu = fits.open('/Users/Sean/research/egs/AEGIS-N2_K_sci.fits')
    elif field=='cosmos':
        imagehdu = fits.open('/Users/Sean/research/cosmos/COSMOS-1_K_sci.fits')
    else:
        print('wrong field')
        
    m1hdu = fits.open(semester+'/maskdesign/finishedmasks/m17B70.fits')
    m2hdu = fits.open(semester+'/maskdesign/finishedmasks/m17B71.fits')
    m3hdu = fits.open(semester+'/maskdesign/finishedmasks/m17B80.fits')
    hdulist = fits.open(semester+'/newfirm_mbs_specz_sdss_priority.fits')
    
    if field=='egs':
        groupdata = Table.read('/Users/Sean/research/deep/groupcat_deimos16a.txt', format = 'ascii')
    elif field=='cosmos':
        grouphdu = fits.open('/Users/Sean/geec/geec_groupdata.fits')
        groupdata = grouphdu[1].data
    else:
        print('wrong group catalog')
        
    w = wcs.WCS(imagehdu[0].header)

    group_ra = np.array(groupdata.RAJ2000)
    group_dec = np.array(groupdata.DEJ2000)
    #group_priority = np.array(groupdata['PRIORITY'])

    m1_data = m1hdu[1].data
    m2_data = m2hdu[1].data
    m3_data = m3hdu[1].data

    m1_ra = m1_data['RA_OBJ']
    m2_ra = m2_data['RA_OBJ']
    m3_ra = m3_data['RA_OBJ']
    m1_dec = m1_data['DEC_OBJ']
    m2_dec = m2_data['DEC_OBJ']
    m3_dec = m3_data['DEC_OBJ']
    m1_pa = m1_data['MajAxPA']
    m2_pa = m2_data['MajAxPA']
    m3_pa = m3_data['MajAxPA']

    nf_data = hdulist[1].data
    nf_ra = nf_data['ra']
    nf_dec = nf_data['dec']
    nf_pflag = nf_data['Pflag']
    print np.str(np.max(nf_pflag))
    nf_target = nf_data['target']
    nf_deep = nf_data['zcosmos']
    nf_ellip = nf_data['K_ellip']
    nf_theta = nf_data['K_theta_J2000']

    nf_theta[np.where(nf_ellip < 0.3)[0]] = 0
    print nf_ellip[:20]
    print nf_theta[:20]
    
    #nf_rflux = nf_data['R']
    #nf_rmag = 25.-2.5*np.log10(nf_rflux)
    bright = (nf_pflag == 10000)^(nf_pflag == 5000)#(nf_rmag <= 24) & (nf_rmag > 0)
    dull = (nf_pflag == 100)^(nf_pflag == 50)#(nf_rmag > 24) & (nf_rmag < 25.5)
    
    sdss_data = hdulist[2].data
    sdss_ra = sdss_data['ra']
    sdss_dec = sdss_data['dec']
    sdss_pflag = sdss_data['Pflag']

    guidecut = sdss_pflag == -1
    aligncut = sdss_pflag == -2
    deeponlycut = (nf_deep >= 0) & (nf_pflag == 1)
    deepcut = (nf_deep >= 0) 
    scicut1 = nf_pflag == 1
    scicut10 = nf_pflag == 10
    lowcut = nf_pflag == 1000
    highcut = (nf_pflag == 10000)^(nf_pflag == 5000)
    membercut = ((nf_pflag == 10000)) ^ ((nf_pflag == 100)) #^ ((nf_pflag == 5000)&(nf_target==1)) ^ ((nf_pflag == 50)&(nf_target==1))
    membercut2 = ((nf_pflag == 5000)) ^ ((nf_pflag == 50))
    
    #Ghigh = group_priority == 1
    #Glow = group_priority == 2
    
    
    sdssxx,sdssyy = w.wcs_world2pix(sdss_ra,sdss_dec,0,ra_dec_order=True)
    nfxx,nfyy = w.wcs_world2pix(nf_ra,nf_dec,0,ra_dec_order=True)
    m1xx,m1yy = w.wcs_world2pix(m1_ra,m1_dec,0,ra_dec_order=True)
    m2xx,m2yy = w.wcs_world2pix(m2_ra,m2_dec,0,ra_dec_order=True)
    m3xx,m3yy = w.wcs_world2pix(m3_ra,m3_dec,0,ra_dec_order=True)
    gxx,gyy = w.wcs_world2pix(group_ra,group_dec,0,ra_dec_order=True)


    ### Possible file names...total shit show...add new name to the bottom
    
    #f = open('2016a_NMBS_DEEP_region.reg','w')
    #f = open('2016a_NMBS_region_pass1.reg','w')
    #f = open('2016a_NMBS_region_pass2.reg','w')
    #f = open('2016a_NMBS_region_stars.reg','w')
    #f = open('2016a_NMBS_region_mask1c.reg','w')
    #f = open('2016a_NMBS_region_mask2c.reg','w')
    #f = open('2016a_NMBS_region_mask3.reg','w')
    #f = open('2016a_NMBS_region_mask4.reg','w')
    #f = open('2016a_NMBS_region_mask10.reg','w')
    #f = open('2016a_NMBS_region_groupmembers.reg','w')
    #f = open('2016a_NMBS_region_PAcheck.reg','w')
    ### 2017B COSMOS Keck/DEIMOS run ###
    #f = open(semester+'/maskdesign/finishedmasks/2017B_COSMOS_regionfile.reg','w')
    #f = open(semester+'/maskdesign/finishedmasks/2017B_COSMOS_regionfile_m70.reg','w')
    #f = open(semester+'/maskdesign/finishedmasks/2017B_COSMOS_regionfile_m71.reg','w')
    #f = open(semester+'/maskdesign/finishedmasks/2017B_COSMOS_regionfile_m80.reg','w')
    f = open(semester+'/maskdesign/finishedmasks/2017B_COSMOS_regionfile_groups.reg','w')

    

    hxx = nfxx[highcut]
    hyy = nfyy[highcut]
    htheta = nf_theta[highcut]
    memberxx = nfxx[membercut] #& deepcut]
    memberyy = nfyy[membercut] #& deepcut]
    memberxx2 = nfxx[membercut2] #& deepcut]
    memberyy2 = nfyy[membercut2] #& deepcut]
    print(len(memberxx))

    lxx = nfxx[lowcut]
    lyy = nfyy[lowcut]
    ltheta = nf_theta[lowcut]
    deeplxx = nfxx[lowcut & deepcut]
    deeplyy = nfyy[lowcut & deepcut]

    #ghxx = gxx[Ghigh]
    #ghyy = gyy[Ghigh]

    #glxx = gxx[Glow]
    #glyy = gyy[Glow]

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

    #print len(hxx)
    #print len(ghxx)

    #print len(lxx)
    #print len(glxx)

    print(len(m1xx))
    print(len(m2xx))
    print(len(m3xx))
    
    #for i in range(len(hxx)):
        #print >>f, 'circle '+np.str(hxx[i])+' '+np.str(hyy[i])+' 15 #color=magenta'

    ### region print statements for each DEIMOS mask ###

    #for i in range(len(m1xx)):
        #print >>f, 'circle '+np.str(m1xx[i])+' '+np.str(m1yy[i])+' 15 #color=yellow width=2'

    #for i in range(len(m2xx)):
        #print >>f, 'circle '+np.str(m2xx[i])+' '+np.str(m2yy[i])+' 20 #color=magenta width=2'

    #for i in range(len(m3xx)):
        #print >>f, 'circle '+np.str(m3xx[i])+' '+np.str(m3yy[i])+' 25 #color=cyan width=2'

    for i in range(len(memberxx)):
        print >>f, 'box '+np.str(memberxx[i])+' '+np.str(memberyy[i])+' 25 25 45#color=green width=2'

    for i in range(len(memberxx2)):
        print >>f, 'box '+np.str(memberxx2[i])+' '+np.str(memberyy2[i])+' 30 30 45#color=red width=2'

    #for i in range(len(lxx)):
        #print >>f, 'circle '+np.str(lxx[i])+' '+np.str(lyy[i])+' 10 #color=cyan'

    #for i in range(len(deeplxx)):
        #print >>f, 'box '+np.str(deeplxx[i])+' '+np.str(deeplyy[i])+' 25 25 45#color=cyan width=2'

    #for i in range(len(ghxx)):
        #print >>f, 'circle '+np.str(ghxx[i])+' '+np.str(ghyy[i])+' 25 #color=red width=2'

    #for i in range(len(glxx)):
        #print >>f, 'circle '+np.str(glxx[i])+' '+np.str(glyy[i])+' 25 #color=blue width=2'

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

        
