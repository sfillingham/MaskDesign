import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u


def color(idlist):

    data1 = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.153-154.rf', format='ascii')
    data2 = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.153-155.rf', format='ascii')
    data3 = Table.read('/Users/spf/newfirmMBS/aegis/aegis-n2.deblend.rfcolors/aegis-n2.deblend.v5.1.155-161.rf', format='ascii')

    id1 = np.array(data1['id'])
    id2 = np.array(data2['id'])
    id3 = np.array(data3['id'])

    U1 = np.array(data1['L153'])
    B1 = np.array(data1['L154'])
    U2 = np.array(data2['L153'])
    V2 = np.array(data2['L155'])
    V3 = np.array(data3['L155'])
    J3 = np.array(data3['L161'])
    
    UBoutput = np.empty(len(idlist))
    UVoutput = np.empty(len(idlist))
    VJoutput = np.empty(len(idlist))
    


    ###########################################################################
    ## Determine the color for every input ID from NMBS
    ###########################################################################

    for ii in range(len(idlist)):
        id = idlist[ii]

        cut1 = np.where(id==id1)[0]
        cut2 = np.where(id==id2)[0]
        cut3 = np.where(id==id3)[0]

        if len(cut1)==0:
            return "UB color missing, object number " + np.str(ii)
        else:
            Umag = U1[cut1]
            Bmag = B1[cut1]

            UBoutput[ii] = -2.5*np.log10(Umag/Bmag)

        if len(cut2)==0:
            return "UV color missing, object number " + np.str(ii)
        else:
            Umag = U2[cut2]
            Vmag = V2[cut2]

            UVoutput[ii] = -2.5*np.log10(Umag/Vmag)

        if len(cut3)==0:
            return "VJ color missing, object number " + np.str(ii)
        else:
            Vmag = V3[cut3]
            Jmag = J3[cut3]

            VJoutput[ii] = -2.5*np.log10(Vmag/Jmag)


    return UBoutput, UVoutput, VJoutput
        
  
            



