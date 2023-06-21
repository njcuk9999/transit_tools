# Transit tools

Simple tools for transits


## Transit Predict

Given the properties of a transiting system, predict a transit light curve closest to a future date.

Set the transit parameters and a `time_find` and the code will give you a plot of the transit closest to that date.

### Example

1. Get the transit parameters from Nasa Exoplanet archives:

https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TransitView/nph-visibletbls?dataset=transits

2. Enter the values in the code

```python
# object name
objname = 'WASP 107b'
# transit midpoint as a astropy time object
time0 = Time(2457584.329746, format='jd')
# transit duration (period) with astropy units
period = 5.721492420 * uu.day
# radius of star with astropy units
rstar = 0.657000 * uu.R_sun
# radius of planet with astropy units
rplanet = 0.924000 * uu.R_jup
# semi major axis with astropy units
semimajoraxis = 0.055300 * uu.AU
# longitude of periastron with astropy units
longofperiastron = 0.000000 * uu.deg
# eccentricity (0=circular, 1=parabolic)
eccentricity = 0.000000
# inclination with astropy units
inclination = np.arcsin(0.999971) * uu.rad
# delta t for before/after light curve
delta_t = 0.1
# give me the best transit to a date
time_find = Time('2023-06-12')
```

3. Run the code

You will get a plot similar to this:
![image](https://github.com/njcuk9999/transit_tools/assets/6008608/8b08c76f-d672-4d7c-9129-209b2ac24f75)


