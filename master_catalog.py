import numpy as np
import astropy.io.fits as fits
from astropy.table import Table

def get_catalog():

    cat_data = Table.read('aegis/aegis-n2.deblend.v5.1.cat', format = 'ascii')
    mstar_data = Table.read('aegis/aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout', format = 'ascii')
    z_data = Table.read('aegis/aegis-n2.deblend.redshifts/aegis-n2.deblend.v5.1.zout', format = 'ascii')

    cat_paramlist = np.array(['id','x','y','ra','dec','Ks','K','H','H2','H1','J','J3','J2','J1','z','I','R','G','U','z_spec','star_flag'])
    mstar_paramlist = np.array(['col1','col7']) #'col1' = 'id', 'col7' = 'lmass'
    z_paramlist = np.array(['id','z_spec','z_peak'])

    #data from the general catalog
    ID = np.array(cat_data['id'])
    x = np.array(cat_data['x'])
    y = np.array(cat_data['y'])
    ra = np.array(cat_data['ra'])
    dec = np.array(cat_data['dec'])
    Ks = np.array(cat_data['Ks'])
    K = np.array(cat_data['K'])
    H = np.array(cat_data['H'])
    H2 = np.array(cat_data['H2'])
    H1 = np.array(cat_data['H1'])
    J = np.array(cat_data['J'])
    J3 = np.array(cat_data['J3'])
    J2 = np.array(cat_data['J2'])
    J1 = np.array(cat_data['J1'])
    z = np.array(cat_data['z'])
    I = np.array(cat_data['I'])
    R = np.array(cat_data['R'])
    G = np.array(cat_data['G'])
    U = np.array(cat_data['U'])
    z_spec = np.array(cat_data['z_spec'])
    star_flag = np.array(cat_data['star_flag'])
    K_ellip = np.array(cat_data['K_ellip'])
    K_theta = np.array(cat_data['K_theta_J2000'])
    N_star = np.array(cat_data['Near_Star'])
    
    pflag = np.empty(len(ID))
    for i in range(len(pflag)):
        if cat_data['star_flag'][i] < 1:
            pflag[i] = 1
        else:
            pflag[i] = -3

    deepz = np.empty(len(ID))
    deepz[:] = -1
    sID = np.empty(len(ID))
    sID[:] = -1
    targ = np.empty(len(ID))
    targ[:] = 0

    #data from the SPS catalog
    msID = np.array(mstar_data['col1'])
    mstar = np.array(mstar_data['col7'])

    #data from the photometric spec catalog
    redID = np.array(z_data['id'])
    redspec = np.array(z_data['z_spec'])
    redpeak = np.array(z_data['z_peak'])


    ###########
    ## Put the above data into a single FITS file
    ###########

    col1 = fits.Column(name='id',format='K',array=ID)
    col2 = fits.Column(name='ra',format='D',array=ra)
    col3 = fits.Column(name='dec',format='D',array=dec)
    col4 = fits.Column(name='Ks',format='E',array=Ks)
    col5 = fits.Column(name='K',format='E',array=K)
    col6 = fits.Column(name='H',format='E',array=H)
    col7 = fits.Column(name='H2',format='E',array=H2)
    col8 = fits.Column(name='H1',format='E',array=H1)
    col9 = fits.Column(name='J',format='E',array=J)
    col10 = fits.Column(name='J3',format='E',array=J3)
    col11 = fits.Column(name='J2',format='E',array=J2)
    col12 = fits.Column(name='J1',format='E',array=J1)
    col13 = fits.Column(name='z',format='E',array=z)
    col14 = fits.Column(name='I',format='E',array=I)
    col15 = fits.Column(name='R',format='E',array=R)
    col16 = fits.Column(name='G',format='E',array=G)
    col17 = fits.Column(name='U',format='E',array=U)
    col18 = fits.Column(name='mstar',format='E',array=mstar)
    col19 = fits.Column(name='z_spec',format='E',array=redspec)
    col20 = fits.Column(name='z_peak',format='E',array=redpeak)
    col21 = fits.Column(name='star_flag',format='I',array=star_flag)
    col22 = fits.Column(name='K_ellip',format='E',array=K_ellip)
    col23 = fits.Column(name='K_theta_J2000',format='E',array=K_theta)
    col24 = fits.Column(name='Near_Star',format='I',array=N_star)
    col25 = fits.Column(name='Pflag',format='I',array=pflag)
    col26 = fits.Column(name='deepz',format='E',array=deepz)
    col27 = fits.Column(name='sID',format='E',array=sID)
    col28 = fits.Column(name='target',format='I',array=targ)

    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28])

    tbhdu = fits.BinTableHDU.from_columns(cols)

    tbhdu.writeto('newfirm_mbs_catalog.fits')
    
