import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmasher as cmr
import cartopy.crs as ccrs
from specfabpy import specfabpy as sf
from numpy import genfromtxt
import pandas as pd
import math
import pickle
from ai_correlator_tg import ai_correlator
import csv

lm, nlm_len = sf.init(4) # L=4 truncation is sufficient in this case
nlm = np.zeros((nlm_len), dtype=np.complex64) # expansion coefficients

# ----------------------------------------------------
# ----------------------------------------------------

def plot_ODF(nlm, lm, ax, geo, cmap='Greys', cblabel='$\psi/N$ (ODF)', lvls=np.linspace(0.0,0.5,9), tickintvl=4, latres=60, plotAxes=False):
   
    # Discretize over S^2
    theta = np.linspace(0,   np.pi,   latres) # co-lat
    phi   = np.linspace(0, 2*np.pi, 2*latres) # lon
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0))
    F[F<0] = 0 # hide numerical/truncation errors
    
    # Plot    
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white for visibility
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max', cmap=cmap, nchunk=5) # "nchunk" argument must be larger than 0 for constant-ODF (isotropy) to be plotted correctly.

    # Add grid lines
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb = plt.colorbar(h, ax=ax, fraction=0.075, aspect=9,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    cb.set_label(cblabel)
    cb.ax.xaxis.set_ticks(lvls, minor=True)
    
    if plotAxes:
        ax.plot([0],[90], marker=r'$z$', ms=9, c='tab:red',   transform=geo) # z axis
        ax.plot([90],[0], marker=r'$y$', ms=9, c='tab:blue',  transform=geo) # y axis
        ax.plot([0],[0],  marker=r'$x$', ms=9, c='tab:green', transform=geo) # x axis

    return h, cb


def discretize_ODF(nlm, lm, latres=60):
    # latres = 60 # latitude resolution on S^2
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi = np.linspace(0, 2 * np.pi, 2 * latres)  # LON
    phi, theta = np.meshgrid(phi, theta)  # gridded
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    _, nlm_len = lm.shape
    F = np.real(np.sum([nlm[ii] * sp.sph_harm(lm[1][ii], lm[0][ii], phi, theta) for ii in np.arange(nlm_len)], axis=0))

    return (F, lon, lat)

def calc_enhancement(Etensor, nlm, i):
    ### Coordinate basis vectors for enhancement-factor matrix
    #if rot==0:
    (e1, e2, e3) = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]))  # x,y,z cartesian basis
    #else:
    #(e1,e2,e3, eigvals) = sf.frame(nlm, 'e') # a2 eigen basis

    ### Transversely isotropic monocrystal parameters for ice
    #  Here the linear model by Rathmann & Lilien (2021) is used
    nprime = 1  # n'=1 => linear grain rheology
    Ecc = 1  # Enhancement factor for compression along c-axis
    Eca = 1e3  # Enhancement factor for shear parallel to basal plane
    alpha = 0.0125  # Taylor--Sachs weight

    ### Calculate enhancement-factor matrix in the basis (e1,e2,e3)
    Eij = sf.Eeiej(nlm, e1, e2, e3, Ecc, Eca, alpha, nprime)  # 3x3 matrix
    Etensor[i, 0] = Eij[0, 0]  # Exx
    Etensor[i, 1] = Eij[0, 1]  # Exy
    Etensor[i, 2] = Eij[0, 2]  # Exz
    Etensor[i, 3] = Eij[1, 1]  # Eyy
    Etensor[i, 4] = Eij[1, 2]  # Eyz
    Etensor[i, 5] = Eij[2, 2]  # Ez

    return Etensor, Eij

########################################################################################################
# coordinate system: x=along flow, y=perpendicular to flow, z=vertical
def twtdata(rot):
    # 1) load structure tensor from NEGIS
    my_data = genfromtxt('a2_cp.csv', delimiter=',')

    Etensor = np.zeros([len(my_data), 6])
    for i in range(0,len(my_data)):
        lhmin = my_data[i, 4]
        lhmax = my_data[i, 5]
        lvert = my_data[i, 6]

        if rot == 1:  # rotate coordinate system from eigenframe to iceflow frame
            theta = my_data[i, 7]
        else:
            theta = 0

        # define a2, assuming that smaller horizontal eigenvalue is in flow-direction
        a2 = np.diag([lhmin, lhmax, lvert])
        aicorr = ai_correlator()
        # calculate array of complex expansion coefficients
        nlm = aicorr.lami_to_nlm(a2, theta)

        # calculate enhancement tensor
        Etensor, Eij = calc_enhancement(Etensor, nlm, i)

    # 4) calculate equivalent dT
    dT_cold, dT_temp = calc_deltaT(Etensor)

    # saving
    Etensor_cp = np.zeros([len(my_data), 25])
    Etensor_cp[:, 0:7] = my_data
    Etensor_cp[:, 7:13] = Etensor
    Etensor_cp[:, 13:19] = dT_cold
    Etensor_cp[:, 19:25] = dT_temp
    headerList = ['lat', 'lon', 'x', 'y', 'lmin', 'lmax', 'lvert', 'Exx', 'Exy', 'Exz', 'Eyy', 'Eyz', 'Ezz',
                  'dTc_exx', 'dTc_exy', 'dTc_exz', 'dTc_eyy', 'dTc_eyz', 'dTc_ezz',
                  'dTt_exx', 'dTt_exy', 'dTt_exz', 'dTt_eyy', 'dTt_eyz', 'dTt_ezz']
    np.savetxt("Etensor_cp.csv", Etensor_cp, delimiter=",")
    file = pd.read_csv("Etensor_cp.csv")
    file.to_csv("Etensor_cp.csv", header=headerList, index=False)
    print(Eij)

def modeldata(rot):
    lm, nlm_len = sf.init(4)  # L=4 truncation is sufficient in this case
    nlm = np.zeros((nlm_len), dtype=np.complex64)

    my_data =genfromtxt('a2_model.csv', delimiter=',')
    Etensor = np.zeros([len(my_data), 6])
    for i in range(0,len(my_data)):
        lhmax = my_data[i, 5]
        lhmin = my_data[i, 4]
        lvert = my_data[i, 6]

        # ensure that eigenvalues are positive
        if lhmin < 0:
            lhmin = abs(lhmin)
            lhmax = 1-lvert-lhmin

        if rot == 1:    # rotate coordinate system from eigenframe to iceflow frame
            theta = my_data[i, 7]
        else:           # calculate eigenhancements
            theta = 0

        a2 = np.diag([lhmin, lhmax, lvert])
        # without ai_correlation (Smooth results)
        a4 = sf.a4_IBOF(a2) # using the a4_IBOF function to reconstruct a4 tensor from Elmer/Ice

        # 2) construct nlm with ai_correlation (this results in jumps in the output because of the definition of single max vs girdle)
        #aicorr = ai_correlator()
        #nlm = aicorr.lami_to_nlm(a2, theta)

        nlm[:6] = sf.a2_to_nlm(a2)  # using a2 only gives different behaviour
        nlm[:15] = sf.a4_to_nlm(a4)

        nlm = sf.rotate_nlm4(nlm, 0, math.radians(theta))  # rotate by "theta" to change from eigenframe to flow frame.

        # 3) calculate Enhancement tensor
        Etensor, Eij = calc_enhancement(Etensor, nlm, i)

    # 4) calculate equivalent dT
    dT_cold, dT_temp = calc_deltaT(Etensor)

    # saving
    Etensor_model = np.zeros([len(my_data), 26])
    Etensor_model[:, 0:8] = my_data
    Etensor_model[:, 8:14] = Etensor
    Etensor_model[:, 14:20] = dT_cold
    Etensor_model[:, 20:26] = dT_temp
    headerList = ['lat', 'lon', 'x', 'y', 'lmin', 'lmax', 'lvert','theta', 'Exx', 'Exy', 'Exz', 'Eyy', 'Eyz', 'Ezz',
                  'dTc_exx', 'dTc_exy', 'dTc_exz', 'dTc_eyy', 'dTc_eyz', 'dTc_ezz',
                  'dTt_exx', 'dTt_exy', 'dTt_exz', 'dTt_eyy', 'dTt_eyz', 'dTt_ezz']

    if rot == 1:
        np.savetxt("Etensor_model_IBOF_rot.csv", Etensor_model, delimiter=",")
        file = pd.read_csv("Etensor_model_IBOF_rot.csv")
        file.to_csv("Etensor_model_IBOF_rot.csv", header=headerList, index=False)
    else:
        np.savetxt("Etensor_model_IBOF.csv", Etensor_model, delimiter=",")
        file = pd.read_csv("Etensor_model_IBOF.csv")
        file.to_csv("Etensor_model_IBOF.csv", header=headerList, index=False)

    print(Eij)

def birefringencedata(rot,type):
    if type == 1:       # Querprofil
        my_data = genfromtxt('a2_birQ.csv', delimiter=',')
    elif type == 2:      # Diagonalprofil
        my_data = genfromtxt('a2_birD.csv',delimiter=',')
    elif type == 3:     # Längsprofil
        my_data = genfromtxt('a2_birL.csv',delimiter=',')

    Etensor = np.zeros([len(my_data), 6])
    for i in range(0,len(my_data)):
        lhmin = my_data[i, 4]
        lhmax = my_data[i, 5]
        lvert = my_data[i, 6]

        if rot == 1:    # rotate coordinate system from eigenframe to iceflow frame
            theta = my_data[i, 7]
        else:
            theta = 0

        a2 = np.diag([lhmin, lhmax, lvert])
        # estimate a4 from a2:
        aicorr = ai_correlator()
        nlm = aicorr.lami_to_nlm(a2, theta)

        # 3) calculate Enhancement tensor
        Etensor, Eij = calc_enhancement(Etensor, nlm, i)

    # 4) calculate equivalent dT
    dT_cold, dT_temp = calc_deltaT(Etensor)

    # saving
    Etensor_bir = np.zeros([len(my_data), 25])
    Etensor_bir[:, 0:7] = my_data
    Etensor_bir[:, 7:13] = Etensor
    Etensor_bir[:, 13:19] = dT_cold
    Etensor_bir[:, 19:25] = dT_temp
    headerList = ['lat', 'lon', 'x', 'y', 'lmin', 'lmax', 'lvert', 'Exx', 'Exy', 'Exz', 'Eyy', 'Eyz', 'Ezz',
                  'dTc_exx', 'dTc_exy', 'dTc_exz', 'dTc_eyy', 'dTc_eyz', 'dTc_ezz',
                  'dTt_exx', 'dTt_exy', 'dTt_exz', 'dTt_eyy', 'dTt_eyz', 'dTt_ezz']
    if type == 1:
        np.savetxt("Etensor_birQ.csv", Etensor_bir, delimiter=",")
        file = pd.read_csv("Etensor_birQ.csv")
        file.to_csv("Etensor_birQ.csv", header=headerList, index=False)
    elif type == 2:
        np.savetxt("Etensor_birD.csv", Etensor_bir, delimiter=",")
        file = pd.read_csv("Etensor_birD.csv")
        file.to_csv("Etensor_birD.csv", header=headerList, index=False)
    elif type == 3:
        np.savetxt("Etensor_birL.csv", Etensor_bir, delimiter=",")
        file = pd.read_csv("Etensor_birL.csv")
        file.to_csv("Etensor_birL.csv", header=headerList, index=False)
    print(Eij)

def calc_deltaT(Etensor):
    T1 = 253.15     # assume cold ice is -20 degree
    T2 = 263.15 # -10 degree
    R = 8.314       # gas constant J/K/mol
    Qc = 60000       # valid for ice colder than -10 degrees, J/mol
    Qt = 152000     # valid for ice warmer than -10 degrees, J/mol
    dT_cold = np.zeros([len(Etensor), 6])
    dT_temp = np.zeros([len(Etensor), 6])

    dT_cold[:, 0] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 0])) / (1 - T1 * R / Qc * np.log(Etensor[:, 0]))     # dT exx
    dT_cold[:, 1] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 1])) / (1 - T1 * R / Qc * np.log(Etensor[:, 1]))     # dT exy
    dT_cold[:, 2] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 2])) / (1 - T1 * R / Qc * np.log(Etensor[:, 2]))     # dT exz
    dT_cold[:, 3] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 3])) / (1 - T1 * R / Qc * np.log(Etensor[:, 3]))     # dT eyy
    dT_cold[:, 4] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 4])) / (1 - T1 * R / Qc * np.log(Etensor[:, 4]))     # dT eyz
    dT_cold[:, 5] = (T1 ** 2 * R / Qc * np.log(Etensor[:, 5])) / (1 - T1 * R / Qc * np.log(Etensor[:, 5]))     # dT ezz

    dT_temp[:, 0] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 0])) / (1 - T2 * R / Qt * np.log(Etensor[:, 0]))    # dT exx
    dT_temp[:, 1] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 1])) / (1 - T2 * R / Qt * np.log(Etensor[:, 1]))    # dT exy
    dT_temp[:, 2] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 2])) / (1 - T2 * R / Qt * np.log(Etensor[:, 2]))    # dT exz
    dT_temp[:, 3] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 3])) / (1 - T2 * R / Qt * np.log(Etensor[:, 3]))    # dT eyy
    dT_temp[:, 4] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 4])) / (1 - T2 * R / Qt * np.log(Etensor[:, 4]))    # dT eyz
    dT_temp[:, 5] = (T2 ** 2 * R / Qt * np.log(Etensor[:, 5])) / (1 - T2 * R / Qt * np.log(Etensor[:, 5]))    # dT ezz

    return dT_cold, dT_temp

def example_plot():
    Etensor = np.zeros([1, 6])
    # shear margin core
    #lhmin = 0.0
    #lhmax = 0.85
    #lvert = 0.15

    # EGRIP:
    lhmin = 0.0
    lhmax = 0.6
    lvert = 0.4

    #downstream: 72 km
    #lhmin = 0.205
    #lhmax = 0.35
    #lvert = 0.445

    #downstream: 116 km
    #lhmin = 0.192
    #lhmax = 0.444
    #lvert = 0.364

    # vertical single max
    #lhmin = 0.05
    #lhmax = 0.05
    #lvert = 0.9

    # iso
    #lhmin = 0.33
    #lhmax = 0.33
    #lvert = 1-lhmin-lhmax

    #trns-girdle
    #lhmin = 0.25
    #lhmax = 0.38
    #lvert = 1-lhmin-lhmax

    #trns-girdle2
    #lhmin = 0.15
    #lhmax = 0.42
    #lvert = 1 - lhmin - lhmax
    #lvert = 0.42

    #trns-sm
    #lhmin = 0.28
    #lhmax = 0.44
    #lvert = 1 - lhmin - lhmax
    #lvert = 0.28

    #trns-sm2
    #lhmin = 0.21
    #lhmax = 0.57
    #lvert = 1 - lhmin - lhmax
    #lvert = 0.21


    #theta = -45-20
    theta = 0
    a2 = np.diag([lhmin, lhmax, lvert])

    aicorr = ai_correlator()
    nlm = aicorr.lami_to_nlm(a2, theta)
    nlmvec = np.zeros(len(nlm)+1, dtype=np.complex64)
    nlmvec[1:17] = nlm
    np.savetxt("nlm.csv", nlmvec, delimiter=",")
    file = pd.read_csv("nlm.csv")
    file.to_csv("nlm.csv", header="nlm", index=False)

    # 3) calculate Enhancement tensor
    Etensor, Eij = calc_enhancement(Etensor, nlm, 0)

    # 4) calculate equivalent dT
    dT_cold, dT_temp = calc_deltaT(Etensor)

    # Setup figure
    fig = plt.figure(figsize=(3, 4))
    #inclination, rot = 45, -45  # view angle
    inclination, rot = 45, -65
    prj, geo = ccrs.Orthographic(rot, 90 - inclination), ccrs.Geodetic()
    ax = plt.subplot(projection=prj)
    ax.set_global()  # show entire S^2
    # Plot
    h, cb = plot_ODF(nlm, lm, ax, geo, plotAxes=False)
    fig.tight_layout()
    plt.show(block=False)
    plt.savefig("exampleCOF.png",transparent=True)

    print(Eij)

modeldata(1)
modeldata(2)
twtdata(0)      # 1: rotate from eigenframe to flowframe, 0: no rotation
birefringencedata(0, type = 1)  # type =1: Querprofil, type=2: Diagonalprofil, type=3: Längsprofil
birefringencedata(0, type = 2)
#example_plot()
