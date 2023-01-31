# [Fit Radial Velocity](@id fit-pma)

You can use Octofitter as a basic tool for fitting radial velocity data, by itself, or in combination with other kinds of data.
Multiple instruments (up to five) are supported.

!!! tip
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. You'll need
    to add both packaged to continue.

The following is an example of jointly fitting the radial velocity and astrometric motion of the star $\epsilon$ Eridani using some of the radial velocity data collated in Mawet et al.

Datasets from two different radial velocity insturments are included and modelled together with separate jitters and instrumental offsets.

For the purposes of this quick example, a fixed semi-major axis, eccentricity, and m*sin(i) are used.


```julia


using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits, Plots

gaia_id = 5164707970261890560 


@named b = Planet{VisualOrbit}(
    Variables(
        e = 0,
        τ = UniformCircular(1.0),
        mass = (sys,pl) -> 0.78*sin(pl.i),
        a=3.48,
        i=Sine(),
        ω=0.0,#UniformCircular(),
        Ω=UniformCircular(),

    ),
    # No planet astrometry is included since it has not yet been directly detected
)

# Convert from JD to MJD
# Data tabulated from Mawet et al
jd(mjd) = mjd - 2400000.5
rvs = RadialVelocity(
    (;inst_idx=1, epoch=jd(2455110.97985),  rv=−6.54, σ_rv=1.30),
    (;inst_idx=1, epoch=jd(2455171.90825),  rv=−3.33, σ_rv=1.09),
    (;inst_idx=1, epoch=jd(2455188.78841),  rv=7.90, σ_rv=.11),
    (;inst_idx=1, epoch=jd(2455231.7593),  rv=−8.39, σ_rv=.13),
    (;inst_idx=1, epoch=jd(2455255.70841),  rv=1.66, σ_rv=.70),
    (;inst_idx=1, epoch=jd(2455260.71231),  rv=1.77, σ_rv=.01),
    (;inst_idx=1, epoch=jd(2455261.71825),  rv=0.75, σ_rv=.30),
    (;inst_idx=1, epoch=jd(2455413.14376),  rv=−10.67, σ_rv=0.76),
    (;inst_idx=1, epoch=jd(2455414.13849),  rv=−16.73, σ_rv=0.99),
    (;inst_idx=1, epoch=jd(2455415.14082),  rv=−20.89, σ_rv=0.78),
    (;inst_idx=1, epoch=jd(2455426.14477),  rv=−17.57, σ_rv=0.86),
    (;inst_idx=1, epoch=jd(2455427.14813),  rv=−18.05, σ_rv=0.87),
    (;inst_idx=1, epoch=jd(2455428.14758),  rv=−21.46, σ_rv=0.87),
    (;inst_idx=1, epoch=jd(2455429.14896),  rv=−18.67, σ_rv=0.90),
    (;inst_idx=1, epoch=jd(2455434.14805),  rv=7.21, σ_rv=.86),
    (;inst_idx=1, epoch=jd(2455435.14705),  rv=4.46, σ_rv=.89),
    (;inst_idx=1, epoch=jd(2455436.14535),  rv=−2.48, σ_rv=0.83),
    (;inst_idx=1, epoch=jd(2455437.15006),  rv=−5.03, σ_rv=0.94),
    (;inst_idx=1, epoch=jd(2455438.15172),  rv=−14.24, σ_rv=0.90),
    (;inst_idx=1, epoch=jd(2455439.14979),  rv=−13.17, σ_rv=0.51),
    (;inst_idx=1, epoch=jd(2455440.15188),  rv=−22.38, σ_rv=0.88),
    (;inst_idx=1, epoch=jd(2455441.15033),  rv=−19.71, σ_rv=0.99),
    (;inst_idx=1, epoch=jd(2455456.01632),  rv=4.52, σ_rv=.97),
    (;inst_idx=1, epoch=jd(2455465.07401),  rv=−12.99, σ_rv=0.98),
    (;inst_idx=1, epoch=jd(2455469.1284),  rv=7.81, σ_rv=1.01),
    (;inst_idx=1, epoch=jd(2455471.97444),  rv=−4.15, σ_rv=1.16),
    (;inst_idx=1, epoch=jd(2455487.00413),  rv=−9.44, σ_rv=0.96),
    (;inst_idx=1, epoch=jd(2455500.98687),  rv=−2.23, σ_rv=1.05),
    (;inst_idx=1, epoch=jd(2455521.89317),  rv=−11.42, σ_rv=1.05),
    (;inst_idx=1, epoch=jd(2455542.95125),  rv=−8.56, σ_rv=1.20),
    (;inst_idx=1, epoch=jd(2455613.70363),  rv=0.65, σ_rv=.01),
    (;inst_idx=1, epoch=jd(2455791.13884),  rv=1.87, σ_rv=.87),
    (;inst_idx=1, epoch=jd(2455792.13464),  rv=−9.19, σ_rv=0.90),
    (;inst_idx=1, epoch=jd(2455793.13858),  rv=−17.85, σ_rv=0.89),
    (;inst_idx=1, epoch=jd(2455795.14053),  rv=−15.43, σ_rv=0.96),
    (;inst_idx=1, epoch=jd(2455797.13828),  rv=−5.67, σ_rv=0.83),
    (;inst_idx=1, epoch=jd(2455798.14195),  rv=−5.00, σ_rv=0.84),
    (;inst_idx=1, epoch=jd(2455807.1116),  rv=−3.91, σ_rv=.99),
    (;inst_idx=1, epoch=jd(2455809.1367),  rv=−0.90, σ_rv=.99),
    (;inst_idx=1, epoch=jd(2455870.9902),  rv=1.81, σ_rv=1.20),
    (;inst_idx=1, epoch=jd(2455902.82961),  rv=4.20, σ_rv=.74),
    (;inst_idx=1, epoch=jd(2455960.69933),  rv=−8.22, σ_rv=1.21),
    (;inst_idx=1, epoch=jd(2456138.12976),  rv=−2.69, σ_rv=0.86),
    (;inst_idx=1, epoch=jd(2456149.05961),  rv=−2.49, σ_rv=0.53),
    (;inst_idx=1, epoch=jd(2456173.13157),  rv=−1.22, σ_rv=0.96),
    (;inst_idx=1, epoch=jd(2456202.99824),  rv=19.64, σ_rv=0.71),
    (;inst_idx=1, epoch=jd(2456327.70174),  rv=20.33, σ_rv=1.05),
    (;inst_idx=1, epoch=jd(2456343.7026),  rv=16.52, σ_rv=.05),
    (;inst_idx=1, epoch=jd(2456530.11763),  rv=6.76, σ_rv=.90),
    (;inst_idx=1, epoch=jd(2456532.12218),  rv=8.06, σ_rv=.85),
    (;inst_idx=1, epoch=jd(2456587.96668),  rv=14.41, σ_rv=1.03),
    (;inst_idx=1, epoch=jd(2456613.91026),  rv=15.04, σ_rv=1.02),
    (;inst_idx=1, epoch=jd(2456637.81493),  rv=23.88, σ_rv=1.02),
    (;inst_idx=1, epoch=jd(2456638.79118),  rv=32.35, σ_rv=1.07),
    (;inst_idx=1, epoch=jd(2456674.80603),  rv=11.70, σ_rv=1.03),
    (;inst_idx=1, epoch=jd(2456708.78257),  rv=2.49, σ_rv=.99),
    (;inst_idx=1, epoch=jd(2456884.13093),  rv=12.85, σ_rv=0.95),
    (;inst_idx=1, epoch=jd(2456889.14678),  rv=18.51, σ_rv=0.82),
    (;inst_idx=1, epoch=jd(2456890.14703),  rv=13.09, σ_rv=0.86),
    (;inst_idx=1, epoch=jd(2456894.13998),  rv=8.71, σ_rv=.83),
    (;inst_idx=1, epoch=jd(2456896.11131),  rv=15.09, σ_rv=0.78),
    (;inst_idx=1, epoch=jd(2456910.94964),  rv=13.84, σ_rv=0.64),
    (;inst_idx=1, epoch=jd(2457234.13834),  rv=9.97, σ_rv=.85),
    (;inst_idx=1, epoch=jd(2457240.99109),  rv=6.26, σ_rv=.52),
    (;inst_idx=1, epoch=jd(2457243.14297),  rv=3.19, σ_rv=.78),
    (;inst_idx=1, epoch=jd(2457245.14532),  rv=5.26, σ_rv=.90),
    (;inst_idx=1, epoch=jd(2457246.14242),  rv=−1.45, σ_rv=0.99),
    (;inst_idx=1, epoch=jd(2457247.14678),  rv=−5.60, σ_rv=1.01),
    (;inst_idx=1, epoch=jd(2457254.14889),  rv=8.50, σ_rv=.80),
    (;inst_idx=1, epoch=jd(2457255.15244),  rv=6.36, σ_rv=.91),
    (;inst_idx=1, epoch=jd(2457256.15168),  rv=5.80, σ_rv=.83),
    (;inst_idx=1, epoch=jd(2457265.14924),  rv=5.74, σ_rv=.88),
    (;inst_idx=1, epoch=jd(2457291.04683),  rv=6.07, σ_rv=.05),
    (;inst_idx=1, epoch=jd(2457326.9831),  rv=6.10, σ_rv=1.12),
    (;inst_idx=1, epoch=jd(2457353.88153),  rv=−0.55, σ_rv=1.09),
    (;inst_idx=1, epoch=jd(2457378.78993),  rv=2.19, σ_rv=.08),
    (;inst_idx=1, epoch=jd(2457384.78144),  rv=14.17, σ_rv=1.10),
    (;inst_idx=1, epoch=jd(2457401.75106),  rv=6.07, σ_rv=.99),
    (;inst_idx=1, epoch=jd(2457669.02614),  rv=1.91, σ_rv=.10),
    (;inst_idx=1, epoch=jd(2457672.99494),  rv=−1.33, σ_rv=1.20),
    (;inst_idx=1, epoch=jd(2457678.97973),  rv=−13.88, σ_rv=1.10),
    (;inst_idx=1, epoch=jd(2457704.03411),  rv=−14.12, σ_rv=0.67),
    (;inst_idx=1, epoch=jd(2457712.99284),  rv=−4.84, σ_rv=1.18),
    (;inst_idx=1, epoch=jd(2457789.74988),  rv=−13.12, σ_rv=1.12),
    (;inst_idx=1, epoch=jd(2457790.737),  rv=−8.09, σ_rv=1.01),
    (;inst_idx=1, epoch=jd(2457803.70407),  rv=−4.25, σ_rv=1.09),
    (;inst_idx=1, epoch=jd(2457804.70718),  rv=−6.55, σ_rv=1.09),
    (;inst_idx=1, epoch=jd(2457806.79201),  rv=−11.62, σ_rv=1.13),
    (;inst_idx=1, epoch=jd(2457828.7545),  rv=−12.69, σ_rv=1.12),
    (;inst_idx=1, epoch=jd(2457829.71875),  rv=−19.82, σ_rv=0.98),
    (;inst_idx=1, epoch=jd(2457830.71979),  rv=−12.66, σ_rv=1.10),




    (inst_idx=2, epoch=jd(2456582.93034), rv=26.64, σ_rv=2.73),
    (inst_idx=2, epoch=jd(2456597.91368), rv=6.40, σ_rv=2.36),
    (inst_idx=2, epoch=jd(2456606.68427), rv=16.52, σ_rv=0.75),
    (inst_idx=2, epoch=jd(2456608.10376), rv=4.69, σ_rv=0.78),
    (inst_idx=2, epoch=jd(2456610.7625), rv=16.04, σ_rv=1.18),
    (inst_idx=2, epoch=jd(2456618.88476), rv=−2.11, σ_rv=0.78),
    (inst_idx=2, epoch=jd(2456624.72004), rv=4.20, σ_rv=1.11),
    (inst_idx=2, epoch=jd(2456626.81421), rv=24.46, σ_rv=0.75),
    (inst_idx=2, epoch=jd(2456628.72976), rv=24.14, σ_rv=0.70),
    (inst_idx=2, epoch=jd(2456631.42746), rv=−2.26, σ_rv=0.88),
    (inst_idx=2, epoch=jd(2456632.80921), rv=14.46, σ_rv=0.62),
    (inst_idx=2, epoch=jd(2456644.75696), rv=8.20, σ_rv=2.30),
    (inst_idx=2, epoch=jd(2456647.81171), rv=14.44, σ_rv=0.63),
    (inst_idx=2, epoch=jd(2456648.59184), rv=12.62, σ_rv=1.10),
    (inst_idx=2, epoch=jd(2456662.63738), rv=9.77, σ_rv=0.73),
    (inst_idx=2, epoch=jd(2456663.75415), rv=10.43, σ_rv=1.11),
    (inst_idx=2, epoch=jd(2456667.52792), rv=18.00, σ_rv=0.78),
    (inst_idx=2, epoch=jd(2456671.68695), rv=19.96, σ_rv=1.05),
    (inst_idx=2, epoch=jd(2456675.75647), rv=7.84, σ_rv=1.12),
    (inst_idx=2, epoch=jd(2456679.83732), rv=17.70, σ_rv=1.05),
    (inst_idx=2, epoch=jd(2456682.56608), rv=17.80, σ_rv=0.82),
    (inst_idx=2, epoch=jd(2456689.76638), rv=26.34, σ_rv=0.75),
    (inst_idx=2, epoch=jd(2456875.02028), rv=7.12, σ_rv=2.18),
    (inst_idx=2, epoch=jd(2456894.88054), rv=8.28, σ_rv=1.30),
    (inst_idx=2, epoch=jd(2456901.06193), rv=9.95, σ_rv=1.54),
    (inst_idx=2, epoch=jd(2456909.10279), rv=−4.71, σ_rv=1.21),
    (inst_idx=2, epoch=jd(2456922.07953), rv=12.25, σ_rv=2.13),
    (inst_idx=2, epoch=jd(2456935.94021), rv=−2.43, σ_rv=1.27),
    (inst_idx=2, epoch=jd(2456937.92403), rv=−0.55, σ_rv=1.35),
    (inst_idx=2, epoch=jd(2456950.03798), rv=3.82, σ_rv=1.44),
    (inst_idx=2, epoch=jd(2456985.64755), rv=−1.80, σ_rv=2.28),
    (inst_idx=2, epoch=jd(2456988.63095), rv=5.93, σ_rv=1.29),
    (inst_idx=2, epoch=jd(2456999.76434), rv=8.84, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457015.72916), rv=−2.17, σ_rv=1.10),
    (inst_idx=2, epoch=jd(2457026.78021), rv=−1.44, σ_rv=1.34),
    (inst_idx=2, epoch=jd(2457058.45996), rv=−3.69, σ_rv=1.89),
    (inst_idx=2, epoch=jd(2457234.08236), rv=7.73, σ_rv=1.39),
    (inst_idx=2, epoch=jd(2457245.86234), rv=−4.19, σ_rv=1.41),
    (inst_idx=2, epoch=jd(2457249.93007), rv=−3.94, σ_rv=1.31),
    (inst_idx=2, epoch=jd(2457253.11257), rv=5.63, σ_rv=1.33),
    (inst_idx=2, epoch=jd(2457257.15719), rv=−1.02, σ_rv=1.15),
    (inst_idx=2, epoch=jd(2457258.94437), rv=−12.69, σ_rv=1.23),
    (inst_idx=2, epoch=jd(2457261.02221), rv=−2.76, σ_rv=1.32),
    (inst_idx=2, epoch=jd(2457262.94505), rv=−7.81, σ_rv=1.36),
    (inst_idx=2, epoch=jd(2457265.95783), rv=9.67, σ_rv=1.24),
    (inst_idx=2, epoch=jd(2457275.01304), rv=−1.91, σ_rv=1.23),
    (inst_idx=2, epoch=jd(2457283.96368), rv=1.88, σ_rv=1.29),
    (inst_idx=2, epoch=jd(2457287.02735), rv=−1.11, σ_rv=1.35),
    (inst_idx=2, epoch=jd(2457290.95635), rv=3.19, σ_rv=1.42),
    (inst_idx=2, epoch=jd(2457305.83659), rv=−5.63, σ_rv=1.23),
    (inst_idx=2, epoch=jd(2457308.90844), rv=13.30, σ_rv=1.26),
    (inst_idx=2, epoch=jd(2457318.83435), rv=8.72, σ_rv=1.26),
    (inst_idx=2, epoch=jd(2457321.79157), rv=6.64, σ_rv=1.36),
    (inst_idx=2, epoch=jd(2457325.84352), rv=2.87, σ_rv=1.41),
    (inst_idx=2, epoch=jd(2457331.10764), rv=9.90, σ_rv=1.36),
    (inst_idx=2, epoch=jd(2457332.78237), rv=9.64, σ_rv=1.25),
    (inst_idx=2, epoch=jd(2457334.82998), rv=5.22, σ_rv=1.30),
    (inst_idx=2, epoch=jd(2457337.7891), rv=5.41, σ_rv=1.59),
    (inst_idx=2, epoch=jd(2457340.95644), rv=−1.99, σ_rv=1.27),
    (inst_idx=2, epoch=jd(2457347.86896), rv=4.10, σ_rv=1.29),
    (inst_idx=2, epoch=jd(2457348.77993), rv=4.65, σ_rv=1.27),
    (inst_idx=2, epoch=jd(2457350.72611), rv=5.83, σ_rv=1.20),
    (inst_idx=2, epoch=jd(2457354.70613), rv=−0.88, σ_rv=1.65),
    (inst_idx=2, epoch=jd(2457361.64656), rv=17.26, σ_rv=1.43),
    (inst_idx=2, epoch=jd(2457364.77113), rv=−7.80, σ_rv=1.30),
    (inst_idx=2, epoch=jd(2457365.70544), rv=0.72, σ_rv=1.26),
    (inst_idx=2, epoch=jd(2457424.71436), rv=−1.68, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457426.63205), rv=3.62, σ_rv=1.42),
    (inst_idx=2, epoch=jd(2457427.38923), rv=3.97, σ_rv=1.17),
    (inst_idx=2, epoch=jd(2457429.72793), rv=2.42, σ_rv=0.90),
    (inst_idx=2, epoch=jd(2457432.60322), rv=6.20, σ_rv=1.25),
    (inst_idx=2, epoch=jd(2457435.69406), rv=−18.61, σ_rv=18.79),
    (inst_idx=2, epoch=jd(2457443.66061), rv=2.25, σ_rv=1.24),
    (inst_idx=2, epoch=jd(2457446.70278), rv=3.96, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457471.55712), rv=5.85, σ_rv=1.63),
    (inst_idx=2, epoch=jd(2457599.93545), rv=−5.69, σ_rv=0.85),
    (inst_idx=2, epoch=jd(2457605.99828), rv=−5.33, σ_rv=1.27),
    (inst_idx=2, epoch=jd(2457607.92844), rv=−24.97, σ_rv=1.39),
    (inst_idx=2, epoch=jd(2457611.16197), rv=−16.02, σ_rv=1.26),
    (inst_idx=2, epoch=jd(2457613.86777), rv=2.47, σ_rv=1.54),
    (inst_idx=2, epoch=jd(2457615.04307), rv=3.50, σ_rv=1.48),
    (inst_idx=2, epoch=jd(2457617.08138), rv=0.91, σ_rv=1.29),
    (inst_idx=2, epoch=jd(2457619.05397), rv=−12.30, σ_rv=1.46),
    (inst_idx=2, epoch=jd(2457621.79772), rv=−13.43, σ_rv=1.57),
    (inst_idx=2, epoch=jd(2457626.10874), rv=0.39, σ_rv=1.33),
    (inst_idx=2, epoch=jd(2457627.95628), rv=−4.92, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457633.96762), rv=−8.24, σ_rv=1.70),
    (inst_idx=2, epoch=jd(2457636.08672), rv=−1.33, σ_rv=1.18),
    (inst_idx=2, epoch=jd(2457637.95848), rv=−7.66, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457643.92459), rv=−14.39, σ_rv=1.33),
    (inst_idx=2, epoch=jd(2457668.93315), rv=−0.83, σ_rv=1.34),
    (inst_idx=2, epoch=jd(2457669.90475), rv=2.76, σ_rv=1.43),
    (inst_idx=2, epoch=jd(2457670.88203), rv=−8.82, σ_rv=1.42),
    (inst_idx=2, epoch=jd(2457674.61398), rv=−5.61, σ_rv=1.42),
    (inst_idx=2, epoch=jd(2457679.98028), rv=−12.42, σ_rv=1.78),
    (inst_idx=2, epoch=jd(2457687.77138), rv=1.17, σ_rv=1.37),
    (inst_idx=2, epoch=jd(2457694.76122), rv=−3.81, σ_rv=1.33),
    (inst_idx=2, epoch=jd(2457696.82099), rv=−5.60, σ_rv=1.32),
    (inst_idx=2, epoch=jd(2457700.96748), rv=−10.84, σ_rv=1.41),
    (inst_idx=2, epoch=jd(2457701.84849), rv=−11.69, σ_rv=1.38),
    (inst_idx=2, epoch=jd(2457702.89789), rv=−14.82, σ_rv=1.22),
    (inst_idx=2, epoch=jd(2457703.82658), rv=−19.89, σ_rv=1.25),
    (inst_idx=2, epoch=jd(2457705.73282), rv=−9.58, σ_rv=1.32),
    (inst_idx=2, epoch=jd(2457707.78376), rv=−9.03, σ_rv=1.24),
    (inst_idx=2, epoch=jd(2457717.79818), rv=−15.06, σ_rv=1.22),
    (inst_idx=2, epoch=jd(2457722.75749), rv=−12.43, σ_rv=2.05),
    (inst_idx=2, epoch=jd(2457728.81592), rv=−7.64, σ_rv=1.67),
    (inst_idx=2, epoch=jd(2457741.79955), rv=−14.52, σ_rv=1.16),
    (inst_idx=2, epoch=jd(2457743.5028), rv=−17.28, σ_rv=1.32),
    (inst_idx=2, epoch=jd(2457745.93451), rv=−17.74, σ_rv=1.31),
    (inst_idx=2, epoch=jd(2457749.71344), rv=−5.63, σ_rv=1.30),
    (inst_idx=2, epoch=jd(2457751.64976), rv=−16.16, σ_rv=1.32),
    (inst_idx=2, epoch=jd(2457753.47716), rv=−12.45, σ_rv=1.30),
    (inst_idx=2, epoch=jd(2457798.55461), rv=−18.91, σ_rv=2.25),
    (inst_idx=2, epoch=jd(2457821.65582), rv=−5.60, σ_rv=1.63),
)

@named ϵEri = System(
    Variables(
        M = 0.78,
        plx = gaia_plx(;gaia_id),
        pmra = Normal(-975, 10),
        pmdec = Normal(20,  10),

        rv0_1 = Normal(0,10),
        rv0_2 = Normal(0,10),
        jitter_1 = truncated(Normal(0,10),lower=0),
        jitter_2 = truncated(Normal(0,10),lower=0),

    ),
    ProperMotionAnomHGCA(;gaia_id),
    rvs,
    b
)

## Build model
model = Octofitter.LogDensityModel(ϵEri; autodiff=:ForwardDiff, verbosity=4) # defaults are ForwardDiff, and verbosity=0

## Sample from chains

results = Octofitter.advancedhmc(
    model, 0.65;
    adaptation =  2000,
    iterations =  5000,
    verbosity = 4,
    tree_depth = 12
)

## Display results
timeplotgrid(results, cmap=:greys, clims=(-0.5,1));
savefig("eps-eri.png")
```

![model plot with astrometry](assets/eps-eri-pma-fit-4.png)
