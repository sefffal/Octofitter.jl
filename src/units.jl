
# The Hipparcos epoch is J1991.25 --- it is not the decimal year 1991.25.
# The proper date corresponds to 1991-04-02 13:30:00.000 ISO time.

# TODO: these should reference the constants stored by PlanetOrbits

const julian_year = 365.25
const sec2jyear = 60*60*24*julian_year
# const siderial_year = 365.256363004
# const tropical_year = 365.2422
const hipparcos_catalog_epoch_mjd = 48348.5625 # Time("J1991.25", format="jyear_str").mjd

const IAU_pc2au = 648_000/π
const IAU_au = 149_597_870_700 # m
const IAU_pc2km = IAU_pc2au*IAU_au/1e3

# TODO: this could move to PlanetOrbits in the next major version
mjd2jd(jd) = jd + 2400000.5
jd2mjd(jd) = jd - 2400000.5

export jd2mjd, mjd2jd