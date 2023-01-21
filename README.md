# Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream
## Flow enhancement factors

These scripts were used to calculate directional flow enhancement factors in the Northeast Greenland Ice Stream (NEGIS) which are published in Gerber et al., (2023), Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream, Nature Communications.

These scripts can be used to calculate enhancement factors from any fabric described by second-order structure tensors.

## Script overview 

1) Estimate 2nd order structure tensor for radar-based methods:

`eigenvalues_icestream.m`:
- read delta lambda for crosspoint analysis (**CP.mat**) and beat-signature analysis (**combined_D.csv**, **combined_Q.csv**, **combined_L.csv**)
- estimate the second-order structure tensor from horizontal anisotropy by making the following assumptions:
  - if inside ice stream: set lx = 0.
  - if less than 2km distance from shear margin set lx = lz = (1-delta lambda)/3.
  - if outside shear zone set lx = 0.1
- save files for crosspoint data (**a2_cp.csv**) and beat-signature analysis (**a2_birD.csv**, **a2_birQ.csv**, **a2_birL.csv**)
- for Elmer/Ice results (**theta_flow_model.shp**), read eigenvalues and flow-angle and save in **a2_model.csv**

2) calculate enhancement factors

`NEGIS_fabric_aicorr.py`
- derive 4th-order tensor from a2 with `nlm = aicorr.lami_to_nlm(a2,theta)`, using `ai_correlator_tg.py`
- calculate enhancement factors with `calc_enhancement(Etensor, nlm, i)`
- calculate equivalent temperature difference with `calc_deltaT(Etensor)`
- save output for beat-signature profiles (**Etensor_birQ.csv**, **Etensor_birL.csv**, **Etensor_birD.csv**) and crosspoint travel-time analysis (**Etensor_cp.csv**)
- For Elmer/Ice, instead of using estimating a4 from a2 as above we use `a4_IBOF(a2)` to obtain 4th order instead. This function reconstructs a4 directly obtained from Elmer/Ice.
- The COF evolution model in Elmer/Ice also simulates the COF orientation, i.e. we can calculate enhancement factors in the eigenframe (deformation along COF principal axes) or in the flowframe (deformation along/across flow direction)
- save output from Elmer/Ice in eigenframe (**Etensor_model_IBOF.csv**) and flowframe (**Etensor_model_IBOF_rot.csv**)

## Additional data

To reconstruct the results obtained from the crosspoint travel-time analysis requires **CP.mat** which can be downloaded from ERDA data archive (https://doi.org/?????/) under *raw_data/crosspoint_traveltime_analysis/*

