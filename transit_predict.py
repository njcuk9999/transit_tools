#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given the properties of a transiting system, predict a transit light curve
closest to a future date

Created on 2023-06-21 at 13:29

@author: cook
"""
import batman
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as uu
from astropy.time import Time


# =============================================================================
# Define variables
# =============================================================================
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


# =============================================================================
# Define functions
# =============================================================================
def get_model(t0, per, rp, rs, a, inc, ecc, w, limb_dark='nonlinear',
              u=(0.5, 0.1, 0.1, -0.1), tmin=-0.025, tmax=0.025, nt=1000):
    # object to store transit parameters
    params = batman.TransitParams()
    # time of inferior conjunction (in jd)
    params.t0 = t0
    # orbital period (days)
    params.per = per.to(uu.day).value
    # planet radius (in units of stellar radii)
    params.rp = rp.to(uu.R_sun).value / rs.to(uu.R_sun).value
    # semi-major axis (in units of stellar radii)
    params.a = a.to(uu.AU).value / rs.to(uu.AU).value
    # orbital inclination (in degrees)
    params.inc = inc.to(uu.deg).value
    # eccentricity
    params.ecc = ecc
    # longitude of periastron (in degrees)
    params.w = w.to(uu.deg) .value
    # limb darkening model
    params.limb_dark = limb_dark
    # limb darkening coefficients [u1, u2, u3, u4]
    params.u = u
    # times at which to calculate light curve
    params.t = np.linspace(tmin, tmax, nt)
    # initializes model
    model = batman.TransitModel(params, params.t)

    return params, model


def mid_transit(t0, per, num):
    return t0 + num * per


def predict_transit(t0, per, num=None, tfind=None):
    # Case 1: we are given a number of transits to predict
    if num is not None:
        return mid_transit(t0, per, num)
    # Case 2: we are given a time and want to find the closest transit (tfind)
    else:
        # deal with closest being None
        if tfind is None:
            raise ValueError('Must give either num or closest')
        # Case 2a: we are given a time before the transit
        if tfind < t0:
            direction = -1
        else:
            direction = 1
        # set the dt values
        dt_best = np.inf
        num_best = 0
        num = 0
        # get the closest transit to tfind
        while True:
            # get the next transit
            tnext = mid_transit(t0, per, num + direction)
            # get the time difference between find time and next transit
            dt_next = tfind - tnext

            # if it gets worse then break
            if abs(dt_next) > dt_best:
                break
            # otherwise update the best
            dt_best = abs(dt_next)
            num_best = num + direction
            # update num
            num += 1
        # return the mid transit time
        return mid_transit(t0, per, num_best)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # find the closest transit to
    closest_time = predict_transit(time0.jd, period.to(uu.day).value,
                                   tfind=time_find.jd)
    # get the transit model
    _params, _model = get_model(0, period, rplanet, rstar, semimajoraxis,
                                inclination, eccentricity, longofperiastron,
                                tmin=-delta_t, tmax=delta_t, nt=1000)
    # get the flux for these parameters
    flux = _model.light_curve(_params)

    # push the time into astropy time
    time = Time(_params.t + closest_time, format='jd')


    fig, frame = plt.subplots(ncols=1, nrows=1)

    frame.plot(time.mjd, flux, color='r')

    pargs = [objname, Time(closest_time, format='jd').iso]
    frame.set_title('Transit Prediction for {0} (Center={1})'.format(*pargs))

    plt.ticklabel_format(useOffset=False)
    plt.show()

# =============================================================================
# End of code
# =============================================================================
