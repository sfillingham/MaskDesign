from astropy import units as u
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
import pdb


##################################################################
##
## This will take as input any mask from DSIMULATOR and check its
## accuracy against the master catalog. 
##
##
##################################################################

def check(filename, mastername, num):

    masterhdu = fits.open(mastername)
    masterdata = masterhdu[1].data
    sdss_data = masterhdu[2].data
    sdss_cols = masterhdu[2].columns
    sID = sdss_data['sID']
    sdss_ID = sdss_data['id']
    sdss_ra = sdss_data['ra']
    sdss_dec = sdss_data['dec']
    
    mcols = masterhdu[1].columns
    #print mcols
    m_ID = masterdata['ID']
    m_pflag = masterdata['Pflag']
    m_rband = masterdata['R']
    m_kband = masterdata['K']
    m_mstar = masterdata['mstar']
    
    hdulist = fits.open(filename)
    data = hdulist[1].data
    cols = hdulist[1].columns
    dsim_ID = data['OBJECT']
    pband = data['pBand']

    print "SDSS Guide and Alignment Stars"
    sdsstargets = data['OBJECT'][np.where(np.core.defchararray.find(data['OBJECT'],'s')==0)[0]]
    print sdsstargets

    for i in range(len(sdsstargets)):
        ID = sID[np.where(sID==sdsstargets[i])[0]]
        ra = sdss_ra[np.where(sID==sdsstargets[i])[0]]
        dec = sdss_dec[np.where(sID==sdsstargets[i])[0]]

        print np.str(ID[0])+' '+np.str(ra[0])+' '+np.str(dec[0])
    
    Pflag10000 = np.array([])
    Pflag1000 = np.array([])
    Pflag100 = np.array([])
    Pflag10 = np.array([])
    Pflag5 = np.array([])
    Pflag1 = np.array([])
    other = np.array([])
    Rflux = np.array([])
    Kflux = np.array([])
    mstar = np.array([])
    

    for i in range(len(dsim_ID)-len(sdsstargets)):

        cut = np.where(m_ID == np.int(dsim_ID[i]))[0]
        #print cut
        #cut = cut[0]
        dsim_flag = m_pflag[cut]

        #pdb.set_trace()

        Rflux = np.append(Rflux, m_rband[cut])
        Kflux = np.append(Kflux, m_kband[cut])
        mstar = np.append(mstar, m_mstar[cut])

        if dsim_flag == 10000:
            Pflag10000 = np.append(dsim_flag,Pflag10000)

        elif dsim_flag == 1000:
            Pflag1000 = np.append(dsim_flag,Pflag1000)

        elif dsim_flag == 100:
            Pflag100 = np.append(dsim_flag,Pflag100)

        elif dsim_flag == 10:
            Pflag10 = np.append(dsim_flag,Pflag10)

        elif dsim_flag == 5:
            Pflag5 = np.append(dsim_flag,Pflag5)

        elif dsim_flag == 1:
            Pflag1 = np.append(dsim_flag,Pflag1)

        else:
            other = np.append(dsim_flag,other)

    print 'tier 1 = '+np.str(len(Pflag10000))
    print 'tier 2 = '+np.str(len(Pflag1000))
    print 'tier 3 = '+np.str(len(Pflag100))
    print 'tier 4 = '+np.str(len(Pflag10))
    print 'tier 5 = '+np.str(len(Pflag5))
    print 'tier 6 = '+np.str(len(Pflag1))
    print 'other = '+np.str(len(other))
    


    ########################################################################
    ## Plot the R band magnitude to check brightness of typical object
    ########################################################################
    Rband = -2.5*np.log10(Rflux) + 25.0
    Rband = Rband[Rflux > 0]
    Kband = -2.5*np.log10(Kflux) + 25.0
    Kband = Kband[Kflux > 0]
    
    plt.figure(num)
    n, bins, patches = plt.hist(Rband[np.where(Rband < 50)[0]], 15, facecolor='red', alpha=0.75)
    #plt.plot(Rband,mstar,'ro')
    plt.savefig('Rmag_distribution_mask'+np.str(np.int(num/10.))+'.pdf')
    plt.show()

    plt.figure(num*10)
    n, bins, patches = plt.hist(Kband[np.where(Rband < 50)[0]], 15, facecolor='green', alpha=0.75)
    #plt.plot(Rband,mstar,'ro')
    plt.savefig('Kmag_distribution_mask'+np.str(np.int(num/10.))+'.pdf')
    plt.show()
    


