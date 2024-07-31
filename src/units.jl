
# The Hipparcos epoch is J1991.25 --- it is not the decimal year 1991.25.
# The proper date corresponds to 1991-04-02 13:30:00.000 ISO time.

const julian_year = 365.25
const sec2jyear = 60*60*24*julian_year
# const siderial_year = 365.256363004
# const tropical_year = 365.2422
const hipparcos_catalog_epoch_mjd = 48348.5625 # Time("J1991.25", format="jyear_str").mjd

const IAU_pc2au = 648_000/π
const IAU_au = 149_597_870_700 # m
const IAU_pc2km = IAU_pc2au*IAU_au/1e3
