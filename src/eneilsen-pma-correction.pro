pro convert_pm,ra1,dec1,parallax1,pmra1,pmdec1,rv1,epoch1,$
ra2,dec2,parallax2,pmra2,pmdec2,rv2,epoch2,epoch2a=epoch2a,$
rvsystem1 = rvsystem1,rvsystem2=rvsystem2
               
;Rigorous updating of coordinates.
;RA, Dec should be in degrees
;parallax in mas
;proper motion in mas/yr
;RV in km/s
;epoch in years

;will also calculates light travel time, returning updated epochs
;(epoch2a) due to change in distance between epoch1 and epoch2.
;epoch2 will be when the light was detected, epoch2a will be the
;"emitted" time accounting for the different positions between epoch1
;and epoch 2.

;convert everything to double

ra1 = double(ra1)
dec1 = double(dec1)
parallax1 = double(parallax1)
pmra1 = double(pmra1)
pmdec1 = double(pmdec1)
rv1 = double(rv1)
epoch1 = double(epoch1)
epoch2 = double(epoch2)

mydtor = !dpi / 180d0
my206265 = 180d0 / !dpi * 60d0 * 60d0
;sec2year = 365.2421897d0 * 24d0 * 60d0 * 60d0
sec2year = 365.25d0 * 24d0 * 60d0 * 60d0
pc2km = 3.08567758149137d13

distance1 = 1000d0 / parallax1


;convert RV to pc/year, convert delta RA and delta Dec to radians/year

dra1 = pmra1 / 1000d0 / my206265 / cos(dec1 * mydtor)
ddec1 = pmdec1 / 1000d0 /my206265
ddist1 = rv1 / pc2km * sec2year

;convert first epoch to x,y,z and dx,dy,dz

x1 = cos(ra1*mydtor) * cos(dec1*mydtor) * distance1
y1 = sin(ra1*mydtor) * cos(dec1*mydtor) * distance1
z1 = sin(dec1*mydtor) * distance1

;Excellent.  Now dx,dy,dz,which are constants

dx = -1d0 * sin(ra1 * mydtor) * cos(dec1 * mydtor) * distance1 * dra1 $
          - cos(ra1 * mydtor) * sin(dec1 * mydtor) * distance1 * ddec1 $
          + cos(ra1 * mydtor) * cos(dec1 * mydtor) * ddist1 

dy = 1d0  * cos(ra1 * mydtor) * cos(dec1 * mydtor) * distance1 * dra1 $
          - sin(ra1 * mydtor) * sin(dec1 * mydtor) * distance1 * ddec1 $
          + sin(ra1 * mydtor) * cos(dec1 * mydtor) * ddist1

dz = 1d0 * cos(dec1 * mydtor) * distance1 * ddec1 + sin(dec1 * mydtor) * ddist1

;Now the simple part:

x2 = x1 + dx * (epoch2-epoch1)
y2 = y1 + dy * (epoch2-epoch1)
z2 = z1 + dz * (epoch2-epoch1)

;And done.  Now we just need to go backward.

distance2 = sqrt(x2^2d0 + y2^2d0 + z2^2d0)

parallax2 = 1000d0/distance2

ra2 = ((atan(y2,x2)/mydtor + 360d0) mod 360d0)
dec2 = asin(z2 / distance2) / mydtor

ddist2 = 1d0 / sqrt(x2^2d0 + y2^2d0 + z2^2d0) * (x2 * dx + y2 * dy + z2 * dz)
dra2 = 1d0 / (x2^2d0 + y2^2d0) * (-1d0 * y2 * dx + x2 * dy)
ddec2 = 1d0 / (distance2 * sqrt(1d0 - z2^2d0 / distance2^2d0)) * (-1d0 * z2 * ddist2 / distance2 + dz)

pmra2 = dra2  * my206265 * 1000d0 * cos(dec2 * mydtor)
pmdec2 = ddec2 * 1000d0 * my206265
rv2 = ddist2 * pc2km / sec2year


;dra1 = pmra1 / 1000d0 * distance1
;ddec1 = pmdec1 / 1000d0 * distance1
;ddist1 = rv1 * 3.24078d-14 * 3.154d7

;stop

;light travel time

delta_time = (distance2 - distance1) * 3.085677d13 / 2.99792d5 ;in seconds
epoch2a = epoch2 - delta_time/3.154d7

return
end