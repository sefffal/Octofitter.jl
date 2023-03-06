using CSV
using StringDistances


function HARPS_observations(target, catalog=datadep"HARPS_RVBank")

    
    rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    target_matched_i = findfirst(==(target), rvbank.target)

    if isnothing(target_matched_i)
        avail_filt = rvbank.target[rvbank.target.!=""] 
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return rvbank[target_matched_i]
   
end

function HARPS_rvs(target, catalog=datadep"HARPS_RVBank")

    rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    target_rows = findall(==(target), rvbank.target)
    table = rvbank[target_rows]

    return RadialVelocityLikelihood(Table(;
        epoch=mjd2jd.(table.BJD),
        inst_idx=ones(Int,size(table,1)),
        rv=-table.RV_mlc_nzp,
        σ_rv=table.e_RV_mlc_nzp,
    ))
end

# HARPS RV Bank detailed Header
#=

J/A+A/???/???    Radial velocities for HAPRS targets           (Trifonov+, 2020)
================================================================================
A public HARPS radial velocity database corrected for systematic errors
       T. Trifonov et al.
      <Astron. Astrophys. ???, ??? (2019)>
      =2020A&A...???..???T
================================================================================
ADC_Keywords: Radial velocities
Keywords: methods: data analysis - planetary systems

Abstract:

See https://arxiv.org/pdf/2001.05942.pdf

Description:
    Time series for radial velocities and activity indicators from HARPS 
    spectrograph are presented. See Trifonov et al. (2020) for a detailed 
    description of the parameters.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl    Records    Explanations
--------------------------------------------------------------------------------
ReadMe            80          .    This file
--------------------------------------------------------------------------------

Byte-by-byte Description of file: HARPS_RVBank.dat
--------------------------------------------------------------------------------
   Bytes Format Units      Label    Explanations
------------------------------------------------------------------------------   
   1-13  F13.5 d          BJD            Barycentric Julian date
  14-22  F8.3  m/s        RV_mlc_nzp     Radial velocity (from mlc.dat, sa, NZP and drift corrected, -pre and -post calculated separately)
  23-29  F6.3  m/s      e_RV_mlc_nzp     Radial velocity error  
  
  14-22  F8.3  m/s        RV_drs_nzp     Radial velocity (from drs.dat, sa, NZP and drift corrected -pre and -post calculated separately)
  23-29  F6.3  m/s      e_RV_drs_nzp     Radial velocity error
  
  14-22  F8.3  m/s        RV_mlc         Radial velocity (from mlc.dat, sa and drift corrected, -pre and -post calculated separately)
  23-29  F6.3  m/s      e_RV_mlc         Radial velocity error

  14-22  F8.3  m/s        RV_drs         Radial velocity (from drs.dat, sa and drift corrected, -pre and -post calculated separately)
  23-29  F6.3  m/s      e_RV_drs         Radial velocity error

  14-22  F8.3  m/s        RV_mlc_j       Radial velocity (from mlc.dat, sa and drift corrected, -pre and -post calculated jointly)
  23-29  F6.3  m/s      e_RV_mlc_j       Radial velocity error
  
  46-54  F8.3  m/s        CRX      Chromatic index
  55-63  F8.3  m/s      e_CRX      Chromatic index error
  64-72  F8.3  m/s*km/s   dLW      Differential line width
  73-81  F8.3  m/s*km/s e_dLW      Differential line width error
  91-99  F8.4  ---        Halpha   Halpha index
100-108  F8.4  ---      e_Halpha   Halpha index error
  91-99  F8.4  ---        NaD1     NaD1 index
100-108  F8.4  ---      e_NaD1     NaD1 index error
  91-99  F8.4  ---        NaD2     NaD2 index
100-108  F8.4  ---      e_NaD2     NaD2 index error
109-112  I3    ---      f_RV       [0,15] Bitwise flag (1)
130-138  F8.3  km/s       FWHM_DRS Full width half maximum
139-147  F8.3  ---        CONTRAST_DRS
148-156  F8.3  m/s        Bisector
157-166  F8.3  km/s       RVGUESS  User guess for RV
167-175  F8.3  m/s*km/s   SNR_DRS  Signal-to-noise ratio in order 55
176-189  F13.5  km/s      BJD_DRS  Barycentric Julian date from DRS 
  82-90  F8.3  km/s       BERV     Barycentric Earth radial velocity
190-198  F8.3  km/s       BERV_DRS Barycentric Earth radial velocity from DRS
199-207  F8.3  m/s        DRIFT    drift measure
208-218  F8.3  m/s      E_DRIFT    drift measure error
219-227  F8.3  m/s        SA        Contribution from secular acceleration
228-239  F8.3  m/s        NZP_mlc       NZP for BJD
240-251  F8.3  m/s        dNZP_mlc      Contribution from intra-night systmatics
252-260  F8.3  ---        TMMEAN    Flux weighted mean point
261-270  F8.3  s          EXPTIME   Exposure time  
298-306  F8.3  m/s        MLCRX     ML Chromatic index (Slope over logarithmic wavelength)
307-315  F8.3  m/s      E_MLCRX     error for MLCRX (slope error)
         A     ---        TIMEID    Identifier for file (time stamp)
         A     ---        DRIFT_LAMP  FP,ThAr?
         A     ---        MASK     DRS MASK for CCF
         A     ---        PROGID   Prog-ID
         A     ---        PROGID   Prog-PI
         F8.2  m/s*km/s   AIRMASS  Airmass
         A     ---        OBJAB
         A     ---        THAR_FP
         A     ---        DPR_TYPE

 

Comments:

1) The formating in this version is still not perfect !!! TBFixed
2) Some entries and details might change.

Note (1): The flags are bitwise where:
    0 - normal (reliable) spectra,
    1 - nosci frame, e.g. calibration files,
    2 - spectra taken in I2 mode,
    4 - spectra taken in "eggs" mode,
    16 - coordinates too much off,
    32 - spectra not within a nautical twilight (daytime),
    64 - spectra with too low S/N,
    128 - spectra with too high S/N.
--------------------------------------------------------------------------------

Acknowledgements:
  Trifon Trifnov <trifonov@mpia.de>

================================================================================
(End)              Trifon Trifnov [MPIA, Germany]     15-August-2019

=#
