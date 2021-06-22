"""
Created by Ciara McGrath 20/03/2019

* DISCLAIMER *
  The S/W remains property of Ciara McGrath and shall not be modified,
  nor shall it be forwarded to third parties without prior written consent

"""




import numpy as np          # need to import modules to be used in every individual module that calls it
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt



def greg2jd(month, day, year):
    """ convert from gregorian calendar to julian date

    :param month:
    :param day:
    :param year:
    :return jdn: julian day number
    """
    y = year
    m = month
    b = 0
    c = 0

    if m <= 2:
        y = y - 1
        m = m + 12

    if y < 0:
        c = -.75

    # check for valid calendar date (5th - 14th Oct 1582 not a valid period)
    if year < 1582:
        pass
    elif year > 1582:
        a = np.floor(y / 100)
        b = 2 - a + np.floor(a / 4)
    elif month < 10:
        pass
    elif month > 10:
        a = np.floor(y / 100)
        b = 2 - a + np.floor(a / 4)
    elif day <= 4:
        pass
    elif day > 14:
        a = np.floor(y / 100)
        b = 2 - a + np.floor(a / 4)
    else:
        print('dates specific within 5th - 14th Oct 1582, which is not a valid date for JD conversion')
        quit() # exit simulation

    jd = np.floor(365.25 * y + c) + np.floor(30.6001 * (m + 1))
    jdn = jd + day + b + 1720994.5

    return jdn


def gast(jdate):
    """ Greenwich apparent sidereal time

    :param jdate: julian date
    :return gst: greenwich siderial time
    """
    dtr = np.pi/180  # degrees to radians
    atr = dtr/3600  # arc second to radians

    # time arguments
    t = (jdate - 2451545) / 36525 # number of julian centuries since 12:00 01 Jan 2000
    t2 = t * t
    t3 = t * t2

    # fundamental trig arguments (modulo 2pi functions)
    l = (dtr * (280.4665 + 36000.7698 * t)) % (2*np.pi)
    lp = (dtr * (218.3165 + 481267.8813 * t)) % (2*np.pi)
    lraan = (dtr * (125.04452 - 1934.136261 * t)) % (2*np.pi)

    # nutations in longitude and obliquity
    dpsi = atr * (-17.2 * np.sin(lraan) - 1.32 * np.sin(2 * l) - 0.23 * np.sin(2 * lp) + 0.21 * np.sin(2 * lraan))
    deps = atr * (9.2 * np.cos(lraan) + 0.57 * np.cos(2 * l) + 0.1 * np.cos(2 * lp) - 0.09 * np.cos(2 * lraan))

    # mean and apparent obliquity of the ecliptic
    eps0 = (dtr * (23 + 26 / 60 + 21.448 / 3600) + atr * (-46.815 * t - 0.00059 * t2 + 0.001813 * t3)) % (2*np.pi)
    obliq = eps0 + deps

    # greenwich mean and apparent sidereal time
    gstm = (dtr * (280.46061837 + 360.98564736629 * (jdate - 2451545) + 0.000387933 * t2 - t3 / 38710000)) % (2*np.pi)
    gst = (gstm + dpsi * np.cos(obliq)) % (2*np.pi)

    return gst



# equations of motion for a spacecraft experiencing J2 only and with no thrust
def eq_of_motion_noman(y, t, mu, Re, J2, incl):
    '''
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
    '''
    # assuming J2 only and no thrust

    a, raan, aol = y    # feeds in current state

    # differential equations
    da_dt = 0
    draan_dt = (-3/2)*a**(-2)*(a**(-3)*mu)**(1/2)*np.cos(incl)*J2*Re**2*(1+(3/2)*a**(-2)*(1+(-3/2)*np.sin(incl)**2)*J2*Re**2)
    daol_dt = (1/64)*a**(-4)*(a**(-3)*mu)**(1/2)*(8*a**2+3*(1+3*np.cos(2*incl))*J2*Re**2)*(8*a**2+3*(3+5*np.cos(2*incl))*J2*Re**2)

    return [da_dt, draan_dt, daol_dt]

# Fixed time step solve for ODEs
def odeFixedStepnoman(eq, y0, tmin, tmax, tstep, mu, Re, J2, incl):
    """ Fixed-step solver for ODEs

    :param eq: equations of motion
    :param y0: initial state vector
    :param tmax: end time
    :param tmin: start time
    :param tstep: step size

    :return sol: state vector at each time step
    :return t: time vector
    """


    t = np.linspace(tmin, tmax, int((tmax - tmin)/tstep)) # linspace takes args of : start, stop, and number of divisions

    sol = odeint(eq_of_motion_noman, y0, t, args = (mu, Re, J2, incl)) # call ode solver and pass in necessary arguments for function

    return sol, t

    

# calculate total time that a target is in view of the satellite
def total_time_in_view(timesteps_list, step):
    window_start_time = timesteps_list[0]    # start time for first window
    time_in_view = 0                    # initialise time in view
    for i, j in zip(timesteps_list, timesteps_list[1:]):      # loop over adjacent values in timesteps
        if j - i > step * 2:                                    # if gap between adjacent values is greater than double the timestep
            time_in_view += i - window_start_time               # add the length of the window to the total time in view
            window_start_time = j                               # set new window start time
    time_in_view += timesteps_list[-1] - window_start_time      # add on last window
    
    return time_in_view  

# calculate latitude and longitude of sub-satellite point
def lat_long_ssp_calc(state, ttotal, Re, incl, flattening, rot_rate, jdate_start):
       # jdate_start = greg2jd(month, day, year)     # calculate julian date of start date for use in gst calculations
       # jdate_end = jdate_start + ttotal / 86400  # add on manoeuvre time as fraction of day to give Julian date of time of interest NOT NEEDED - want gst of start
       
       a_total, RAAN_total, u_total = state 
       gst_start = gast(jdate_start)       # calculate gst for time of interest
       
       lat_ssp_c = np.arcsin(np.sin(incl) * np.sin(u_total))    # calculate geocentric latitude of sub-satellite point
       lat_ssp = np.arctan(np.tan(lat_ssp_c)/(1 - flattening * (2 - flattening)))   # convert geocentric latitude of ssp to geodetic
       long_ssp = np.arctan2( np.cos(incl) * np.sin(u_total) , np.cos(u_total) ) - rot_rate * ttotal + RAAN_total - gst_start        # calculate longitude of ssp (need to do gst as function of time....)

       return(lat_ssp, long_ssp)


# calculate haversince curved surface distance between subsatellite point and target
def find_dist(lat_ssp, long_ssp, lat_poi, long_poi, Re):

       # calculate curved surface distance from sub-satellite point to target (need to do at each time step for moving target...)
       dist_to_target = 2 * Re * np.arcsin(np.sqrt((np.sin((lat_ssp - lat_poi)/(2)))**2 + np.cos(lat_ssp) * np.cos(lat_poi) * (np.sin((long_ssp - long_poi)/(2)))**2))

       return(dist_to_target)


# find satellite elevatio compared to target
def find_sat_elevation(a, dist, Re):
    angle_between_ssp_and_tar = dist / Re     # assuming the Earth is a sphere, calculate and between the line of the Earth's centre and subsatellite point and the line of the target and Earth's centre. Treat the distance as an arc and the Earth's radius as the radius of the sector.
    
    earth_angular_radius = np.arcsin(Re / a)
    nadir_angle = np.arctan((np.sin(earth_angular_radius) * np.sin(angle_between_ssp_and_tar)) / (1 - np.sin(earth_angular_radius) * np.cos (angle_between_ssp_and_tar)))

    elevation = np.pi / 2 - angle_between_ssp_and_tar - nadir_angle
    return elevation
