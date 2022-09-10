
import os, sys

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.table import join, Table
# from astroquery.xmatch import XMatch
from astropy import units as u
from scipy.interpolate import interp1d


cmap2 = cm.get_cmap("viridis",7)
colors = {"IC_2391": cmap2(0),
         "IC_2602": cmap2(4),
         "NGC_2547": cmap2(3),
         "NGC_2451A": cmap2(2),
         "Collinder_135": cmap2(1)}


shapes= {"IC_2391": "o",
         "IC_2602": "d",
         "NGC_2547": "v",
         "NGC_2451A": "^",
         "Collinder_135": "s"}


if __name__=="__main__":

    # Interpolation to get masses
    mfile = at.read(os.path.expanduser("~/Dropbox/Models/mamajek_colors.dat.txt"),
                    fill_values=[('...', '0')])
    print(mfile.dtype)

    bv_mass = interp1d(mfile["B-V"],mfile["Msun"],bounds_error=False)
    vi_mass = interp1d(mfile["V-Ic"],mfile["Msun"],bounds_error=False)
    vk_mass = interp1d(mfile["V-Ks"],mfile["Msun"],bounds_error=False)


    # h per
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/hper_rotation_moraux2013.tsv")
    per = at.read(pfile,delimiter="|",data_start=3)
    plt.plot(per["Mass"],per["Per"],'.',color="grey",label="h Per (13 Myr; Moraux+ 2013)")
    plt.legend(loc=2)


    plt.ylim(0.1,50)
    plt.xlim(1.3,0.1)
    plt.yscale("log")

    # plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
    plt.xlabel(r"Mass (M$_{Sun}$)")
    plt.ylabel("Period (d)")
    plt.savefig("0013Myr_hper.png",bbox_inches="tight")
    plt.close()

    # Literature ZAMS clusters
    # IC 2391
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2391_rotation_patten1996.csv")
    per = at.read(pfile,delimiter=",")
    mass = vi_mass(per["V-I"])
    plt.plot(mass,per["Period"],shapes["IC_2391"],color=colors["IC_2391"],label="IC 2391 (Patten & Simon 1996)")

    # IC 2602
    # Barnes
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_barnes1999.csv")
    per = at.read(pfile,delimiter=",")
    mass = bv_mass(per["B-V"])
    plt.plot(mass,per["Prot"],shapes["IC_2602"],color=colors["IC_2602"],label="IC 2602 (Barnes+ 1999)")

    # Tschape & Rudiger
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_tschape2001.csv")
    per = at.read(pfile,delimiter=",")
    mass = bv_mass(per["B-V"])
    plt.plot(mass,per["Prot"],'D',color="#174f34",label="IC 2602 (Tschape & Rudiger 2001)")

    # NGC 2547
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/ngc2547_rotation_irwin2008b.tsv")
    per = at.read(pfile,delimiter="|",data_start=3)
    plt.plot(per["Mass"],per["Per"],shapes["NGC_2547"],color=colors["NGC_2547"],
             label="NGC 2547 (Irwin+ 2008)")

    plt.legend(loc=1)#,ncol=2)
    plt.ylim(0.1,50)
    plt.xlim(1.3,0.1)
    plt.yscale("log")

    # plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
    plt.xlabel(r"Mass (M$_{Sun}$)")
    plt.ylabel("Period (d)")
    plt.savefig("0035Myr_ZAMS.png",bbox_inches="tight")
    plt.close()


    # Pleiades
    pfile = os.path.expanduser("~/Dropbox/data/catalogs/pleiades_rotation_rebull2016.csv")
    per = at.read(pfile,delimiter=",")
    mass = vk_mass(per["(V-K)0"])
    plt.plot(mass,per["Per1"],".",color="grey",label="Pleiades (125 Myr; Rebull+ 2016)")
    plt.legend(loc=2)


    plt.ylim(0.1,50)
    plt.xlim(1.3,0.1)
    plt.yscale("log")

    # plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
    plt.xlabel(r"Mass (M$_{Sun}$)")
    plt.ylabel("Period (d)")
    plt.savefig("0125Myr_pleiades.png",bbox_inches="tight")
    plt.close()
