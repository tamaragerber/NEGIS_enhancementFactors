# N. Rathmann, 2022

import numpy as np
import pickle
from specfabpy import specfabpy as sf
import scipy.special as sp
import math

class ai_correlator():

    def __init__(self):
     
        # Load the save polynomial fit for n_4^0(n_2^0)
        self.x_model, self.y_model, self.pc, self.p_model = pickle.load(open("latrot-validation-corrpoly.p", "rb"))
        self.lm, self.nlm_len = sf.init(4) # use L=4 truncation, nothing larger is needed for this purpose.

    def lami_to_nlm(self, a2, theta):

        # eigenvalues defined in coordinate system of eigenvectors
        lamx = a2[0, 0]
        lamy = a2[1, 1]
        lamz = a2[2, 2]

        # eigenvalues sorted by size
        lam1, lam2, lam3 = np.sort([lamx, lamy, lamz])  # sorted such that lam1 <= lam2 <= lam3
        #print(round(np.sum([lam1, lam2, lam3]),3))

        if round(np.sum([lam1, lam2, lam3]), 3) != 1:
            #print(round(np.sum([lam1, lam2, lam3]), 3) )
            lam3 = 1-lam1-lam2
            #raise ValueError('sum(lami) != 1')

        # single max.: lam_3 > 2/3 and  lam2 ~ lam1 < 1/6

        if lam3 >= 2/3 and lam1 <= 1/6 and lam2 <= 1/6: # conditions for single maximum
        #if lam2 - lam1 <0.2:
        #if lam3 > lam2:
        #if lam3-lam2 > 0.4:
       # if lam3 >0.5:
            a2_vert = np.diag([lam1,lam2,lam3]) # set as vertical single max.
            if lam3 == lamz:
                phi = 0     # if vertical already
            elif lam3 == lamx:
                phi = np.pi     # rotation angle for horizontal single max
            else:
                phi = np.pi/2   # rotation angle for horizontal single max

        # else assume is girdle: lam_3 ~ lam_2 > lam1
        else:
            if lam1 == lamz:    # for horizontal girdle
                a2_vert = a2
                phi = 0
            elif lam1 == lamx:  # for vertical girdle in y
                phi = np.pi
                # rotate to the vertical, accounting for non-symmetric girdles
                if lam3 == lamy:
                    a2_vert = np.diag([lam2, lam3, lam1])
                else:
                    a2_vert = np.diag([lam3, lam2, lam1])
            else:               # for vertical girdle in x
                phi = np.pi/2
                # rotate to the vertical, accounting for non-symmetric girdles
                if lam3 == lamy:
                    a2_vert = np.diag([lam3, lam2, lam1])
                else:
                    a2_vert = np.diag([lam2, lam3, lam1])


        nlm = np.zeros((self.nlm_len), dtype=np.complex64)
        nlm[:6] = sf.a2_to_nlm(a2_vert) # get spectral coefficients for l=2 modes
        n2m = nlm[3] # this is the psi_2^0 component
        n4m = self.p_model(n2m) # this is the psi_4^0 component from the data+model correlation curve
        nlm[10] = n4m # insert the coefficient into its correct position in the coefficient array.

        # if phi != 0, then rotate the fabric pattern down into the horizontal x--y plane, 
        #   followed by "phi" about the vertical z-axis.
        if phi != 0:
            nlm = sf.rotate_nlm4(nlm, -np.pi/2, 0) # rotate into z-axis into x-axis
            nlm = sf.rotate_nlm4(nlm, 0, phi) # rotate by "phi" about z-axis
            nlm = sf.rotate_nlm4(nlm, 0, math.radians(theta)) # rotate by "theta" to change from eigenframe to flow frame.
                
        return nlm
            
#----------------------------------------------------            
#----------------------------------------------------
            
# For plotting below if run as standalone script
def plot_ODF(nlm, lm, ax=None, cmap='Greys', cblabel='$\psi/N$ (ODF)', rot0=-40, lvls=np.linspace(0.0,0.6,9), tickintvl=4, latres=60):
    
    import scipy.special as sp

    import matplotlib.pyplot as plt
    from matplotlib import rcParams, rc, colors
    from matplotlib.offsetbox import AnchoredText
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as mticker

    import cmasher as cmr
    import cartopy.crs as ccrs

    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    
    pltshow = (ax is None)
    
    if ax is None:
        size = 1.5
        plt.figure(figsize=(1.6*size,2*size))
        inclination = 45 # view angle
        rot = rot0 -90 # view angle
        prj = ccrs.Orthographic(rot, 90-inclination)
        geo = ccrs.Geodetic()
        ax = plt.subplot(projection=prj)
        ax.set_global()
    
    F, lon,lat = discretize_ODF(nlm, lm, latres=latres)
    F[F<0] = 0 # fix numerical/truncation errors
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white.
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap=cmap, nchunk=5) # "nchunk" argument must be larger than 0 for constant-ODF (e.g. isotropy) is plotted correctly.
    #ax.set_facecolor(color_bad) # "bad" (masked) values, default white

    # Add grid lines
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb1 = plt.colorbar(h, ax=ax, fraction=0.075, aspect=9,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    cb1.set_label(cblabel)
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    
    if pltshow: 
        plt.tight_layout()
        plt.show()
    
    return h
    
def discretize_ODF(nlm, lm, latres=60):

    #latres = 60 # latitude resolution on S^2        
    theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
    phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    
    _,nlm_len = lm.shape
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0))
    
    return (F, lon,lat)
    

# This is how you can use the correlator class for getting "nlm" given the a2 eigenvalues
if __name__ == "__main__":

    # This part executes only if run as a (standalone) script           

   # plot_ODF(nlm, aicorr.lm)
    aicorr = ai_correlator()
    
    lami = np.flipud(np.array([1, 0, 0])) # a single max
    #lami = np.flipud(np.array([0.5, 0.5, 0])) # a girdle
    nlm = aicorr.lami_to_nlm(lami, phi=4*np.pi/4) # if phi != 0, then the resulting ODF will be rotated for you (see above code documentation).
    # ... you can now calculate enhancement factors using sf.Eij(nlm, ...)
#    print(nlm)
    
    plot_ODF(nlm, aicorr.lm)
    
    
