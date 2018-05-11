from __future__ import print_function
# UCSCScheduler_V3.py
# now featuring Dynamic Priorities(tm)!

from datetime import datetime, timedelta
import ephem
import gspread
import json
#from oauth2client.client import SignedJwtAssertionCredentials
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter, getEXPTime
import getpriority

import numpy as np
import os
import pickle
import sys
import time
try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import re
from astropy.table import Table
from astropy.io import ascii

# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20
TARGET_ELEVATION_PREF_MIN = 45
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX =  3* 60 * 60 # 1 hour
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

# Maximum single exposure time in seconds
MAX_EXPTIME = 1200.
MIN_EXPTIME = 900.
MIN_SHUTTERTIME = 60
MAX_I2 = 40000
MIN_I2 = 2000
MAX_EXPMETER = 2e9


SLOWDOWN_MIN = 0.4
SLOWDOWN_MAX = 10.0

###
# arGGGGGG!!
###

last_objs_attempted = []

def computeMaxTimes(sn,exp_times):
    maxtimes = np.zeros_like(exp_times)
    maxtimes += TARGET_EXPOSURE_TIME_MAX
    for i,n in enumerate(sn):
        if re.match("K2APF",str(n)):
            maxtimes[i] *= 3.
    
    return maxtimes

def parseStarlist(starlist):
    """ Parse a scriptobs-compatible starlist for the scheduler.

    names, star_table, lines, stars = parseStarlist(starlist)
    starlist - a filename

    names - a list of stars in the starlist
    star_table - a numpy array
    lines - a list of strings that can be used for scriptobs input
    stars - a list of pyEphem objects 
    """
    names = []
    lines = []
    stars = []
    star_table = []
    try:
        f = open(starlist,'r')
    except IOError:
        apflog("Warning: Could not open %s. No target can be selected." % starlist,echo=True)
        return None
    else:
        for line in f:
            if not re.search("\A\#",line):
                ls = line.split()
                names.append(ls[0])
                row = []
                # RA value in radians
                row.append(getRARad(ls[1], ls[2], ls[3]))
                # Dec value in radians
                row.append(getDECRad(ls[4], ls[5], ls[6]))
                # PM RA
                row.append(float(ls[8].split('=')[-1]))
                # PM Dec
                row.append(float(ls[9].split('=')[-1]))
                # V mag
                row.append(float(ls[10].split('=')[-1]))
                # Exposure time
                row.append(float(ls[11].split('=')[-1]))
                # Desired Counts
                row.append(float(ls[16].split('=')[-1]))
                # Filler not used here
                row.append(0.)
                row.append(0.)
                # Number of exposures
                row.append(int(ls[19].split('=')[-1]))

                star_table.append(row)

                # Save the scriptobs line for later
                lines.append(line)

                # Generate a pyEphem object for this target
                star = ephem.FixedBody()
                star._ra = ephem.hours(":".join([ls[1], ls[2], ls[3]]))
                star._dec = ephem.degrees(":".join([ls[4], ls[5], ls[6]]))
                stars.append(star)
            
    return names, np.array(star_table), lines, stars

def get_spreadsheet(sheetn="FakeGoogledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """ Get the spreadsheet from google

    worksheet = get_spreadsheet(sheetn="FakeGoogledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')
    worksheet - the worksheet object returned by the gspread module

    sheetn - name of the google sheet, defaults to "FakeGoogledex"
    certificate - certificate used to control access to the google sheet
    
    """
    # this downloads the googledex from the Google Drive
    # the certificate must be available
    # these certificates are generated through the Google Developer Interface
    # the developer must select the correct API for access

    # the certificate has an email associated with it, that email must
    # have the document shared with it to allow access 

    certificate_path = os.path.dirname(__file__)
    
    json_key = json.load(open(os.path.join(certificate_path, certificate)))
    scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gs = gspread.authorize(credentials)

    apflog("Successfully logged in.", echo=True)
    spreadsheet = gs.open(sheetn)
    apflog("Loaded Main %s" % (sheetn),echo=True)
    worksheet = spreadsheet.sheet1
    apflog("Got spreadsheet", echo=True)

    return worksheet

def findColumns(col_names,req_cols):
    """ findColumns finds the indices for the column names in the list of required columns
    indices = findColumns(col_names, req_cols)
    
    indices - a list of indices, each index maps to where in col_names the column is found and in the order of req_cols
    col_names - list of column names to be searched
    req_cols - list of names that should be in the first list
    """
    idx = []
    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) , level="Alert",echo=True)

    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in didx.keys():
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Error",echo=True)

    return didx

def write_as_table(full_codex, outfn):
    hdr = full_codex[0]
    success =  False
    try:
        ascii.write(full_codex[1:],outfn,names=full_codex[0],delimiter=",")
        success = True
    except:
        apflog("Cannot write %s" % (outfn), level="warn")
    return success

def parseGoogledex(sheetn="FakeGoogledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="./newgoogledex.csv"):
    """ parseGoogledex parses the google sheet and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseGoogledex(sheetn="FakeGoogledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat")
    names - a list of stars in the starlist
    star_table - a numpy array
    do_flag - a list of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    try:
        alldata = ascii.read(os.path.join(os.getcwd(),outfn))
    except IOError:
        apflog( "Starting Googledex parse",echo=True)
        worksheet = get_spreadsheet(sheetn=sheetn,certificate=certificate)
        full_codex = worksheet.get_all_values()
        
        #time = (datetime.now() - start).total_seconds()
        #print "Loaded Values. Took {0:f} seconds.".format(time)
        ofn = os.path.join(os.getcwd(),"newgoogledex.csv")
        alldata = write_as_table(full_codex, ofn)
        
    

    good = ((alldata['est_K'] > 0.0) & (alldata['vmag'] < 11) & (alldata['vsini'] < 5)) | (alldata['comments'] == 'B star')
    data = alldata[good]

    names = data['starname']
    do_flag = []
    stars = []
    ras = []
    decs = []

    # Build the star table to return to 
    for ind in range(0,len(data)):
        star = ephem.FixedBody()
        star._ra = ephem.hours(":".join([str(data['rahr'][ind]),str(data['ramin'][ind]),str(data['rasec'][ind])]))
        star._dec = ephem.degrees(":".join([str(data['decdeg'][ind]),str(data['decmin'][ind]),str(data['decsec'][ind])]))
        star.name = str(data['starname'][ind])
        stars.append(star)
        do_flag.append(data['closecompanion'][ind])
        ras.append(getRARad(str(data['rahr'][ind]),str(data['ramin'][ind]),str(data['rasec'][ind])))
        decs.append(getDECRad(str(data['decdeg'][ind]),str(data['decmin'][ind]),str(data['decsec'][ind])))
    data['ra'] = np.asarray(ras)
    data['dec'] = np.asarray(decs)
    data['precision'] = data['est_K'] * 0.5
    data['precision'][data['precision'] < 1.0] = 1.0
    return (names, data, do_flag, stars)
    
def update_googledex_lastobs(filename, sheetn="FakeGoogledex",time=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
        update_googledex_lastobs(filename, sheetn="FakeGoogledex",time=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')

        filename - where the observations are logged
    """
    names, times = getObserved(filename)
    if len(names) == 0:
        return
    if time is None:
        time = datetime.utcnow()
    

    ws = get_spreadsheet(sheetn=sheetn,certificate=certificate)
    vals = ws.get_all_values()

    col = vals[0].index("lastobs") 
    
    for i, v in enumerate(vals):
        # Did we observe this target tonight?
        if v[0] in names:
            # We observed this target, so update the cell in the worksheet
            # update_cell(row, col, val) - col and row are 1 indexed
            otime = times[names.index(v[0])]
            if isinstance(otime,float):
                t = datetime.fromtimestamp(otime)
            else:
                hr, min = otime
                t = datetime(time.year, time.month, time.day, hr, min)
            jd = float(ephem.julian_date(t))
            try:
                pastdate = float(v[col])
                if jd > pastdate:
                    ws.update_cell(i+1, col+1, round(jd, 2) )
            except:
                print (v[0], v[col])
                
    apflog( "Updated Googledex",echo=True)

def update_local_googledex(time,googledex_file="newgoogledex.csv", observed_file="observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time.
        update_local_googledex(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """
    names, times = getObserved(observed_file)

    phase_edges = [ getpriority.EDGE1, getpriority.EDGE2, getpriority.EDGE3, getpriority.EDGE4, getpriority.EDGE5]

    try:
        full_codex = ascii.read(googledex_file)
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
        return names,times

    if type(time) != datetime:
        time = datetime.utcnow()

    if len(names):
        name = names[0]
        # We have observed this star, so lets update the lastobs / phasebins
        otime = times[names.index(name)]
        if isinstance(otime,float):
            t = datetime.fromtimestamp(otime)
        else:
            ohr, omin = otime
        t = datetime(time.year, time.month, time.day, ohr, omin)
        
        jd = float(ephem.julian_date(t))
        sidx, = np.where(full_codex['starname'] == name)
        full_codex['lastobs'][sidx] = jd
        if full_codex['initialphase'][sidx] > 0:
            phases = getpriority.compute_currentphase(jd,full_codex['foldperiod'][sidx],full_codex['initialphase'][sidx])
            binnedphase = np.average(phases)
            try:
                length = len(binnedphase)
                binnum = np.digitize(binnedphase, phase_edges)
            except:
                binnum = np.digitize([binnedphase], phase_edges)
            if binnum > 3:
                binnum = 0
            binname = "phase%dbin" % (binnum+1)
            
            full_codex[binname][sidx] += 1
            apflog( "Updating local googledex star %s phase bin %s to %d" % (name, binname, int(full_codex[binname][sidx])),echo=True)

        ascii.write(full_codex,googledex_file,delimiter=",")
            
    return names, times


def getLST(date, longitude):
    """Take a datetime and longitude and calculate the Local Sidereal Time."""
    # Assumes date is a datetime object, and that the longitude is formatted as in PyEphem 

    ll = [float(v) for v in longitude.split(':')]
    if ll[0] > 0:
        sign = 1
    else:
        sign = -1
    ut = date.hour + date.minute/60. + date.second/3600.
    lng = ll[0] + sign*ll[1]/60. + sign*ll[2]/3600.
    d  = ephem.julian_date() - 2451545.0
    lst = 100.46 + 0.985647 * d + lng + 15*ut
    return lst % 360.


def getRADeg(hr, mn, sec):
    ra_hr = float(hr) + float(mn)/60. + float(sec)/3600.
    return ra_hr * 15

def getRARad(hr, mn, sec):
    ra_hours = float(hr) + float(mn)/60. + float(sec)/3600.
    return ra_hours * 15 * np.pi/180.0

def getDECDeg(deg, mn, sec, neg=False):
    if float(deg) < 0:
        neg = True
    elif float(mn) < 0:
        neg = True
    elif float(sec) < 0:
        neg = True
    x = abs(float(deg)) + float(mn)/60. + float(sec)/3600.
    if neg:
        return x*-1
    else:
        return x

def getDECRad(deg, mn, sec, neg=False):
    if float(deg) < 0:
        neg = True
    elif  float(mn) < 0:
        neg = True
    elif float(sec) < 0:
        neg = True
    x = abs(float(deg)) + float(mn)/60. + float(sec)/3600.
    x = x * np.pi/180.
    if neg:
        return x*-1
    else:
        return x

def getCoordStr(floatval,isRA=False):

    neg = False
    nround = 2
    if isRA:
        floatval /= 15.
        nround = 3
    if floatval < 0:
        neg = True
    floatval = abs(floatval)
    deghrval = int(floatval)
    minval = (floatval % 1) * 60.0 
    secval = round( (minval % 1) *60.0, nround)

    if neg:
        ret = "-" + str(deghrval) + ' '
    else:
        ret = str(deghrval) + ' '
    ret += str(int(minval)) + ' ' + str(secval)
    return ret

        
def getElAz(ra, dec, lat, lng, time):
    """Given RA, DEC, Latitude, and a time, returns the corresponding elevation and azimuth angles
       Works with single values, or numpy arrays
       """
    lst = getLST(time, lng)
    ha = ((lst- np.degrees(ra)) % 360.) * np.pi/180.
    el = np.arcsin(np.sin(dec) * np.sin(lat) + \
                   np.cos(dec) * np.cos(lat) * np.cos(ha))
    az = np.arccos( (np.sin(dec) - np.sin(el)*np.sin(lat)) / \
                         (np.cos(el) * np.cos(lat)))
    return (np.degrees(el), np.degrees(az))

def makeScriptobsLine(res, do_flag, t, decker="W",I2="Y",owner='Vogt'):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y")
    name - name of star, first column in line
    row - star_table row for star that begins with name, cotains all of the data needed for the line except
    do_flag - a string for whether or not scriptob needs to do a pointing check before slewing to the target
    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    """
    focval = 0
    if res['PRI'] > 9.9:
        focval = 2
    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = str(res['NAME']) + ' '
    # Add the RA as three elements, HR, MIN, SEC
    rastr = getCoordStr(np.degrees(res['RA']),isRA=True)
    ret += rastr + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = getCoordStr(np.degrees(res['DEC']))
    ret += decstr + ' '
    # Epoch
    ret += '2000 '
    # Proper motion RA and DEC
    ret += 'pmra=' + str(res['PM_RA']) + ' '
    ret += 'pmdec=' + str(res['PM_DEC']) + ' '
    # V Mag
    ret += 'vmag=%.2f ' % (res['VMAG'])
    # T Exp
    if res['EXP_TIME'] > MAX_EXPTIME:
        ret += 'texp=%d ' % (int(MAX_EXPTIME))
    else:
        ret += 'texp=%d '  % (int(res['EXP_TIME']))
    # I2
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=%d ' % (t.hour) 
    ret += 'utm=%d ' % (t.minute) 
    # Exp Count
    if res['COUNTS'] > 3e9:
        ret += 'expcount=%.3g' % (3e9) + ' '
    else:
        ret += 'expcount=%.3g' % (res['COUNTS']) + ' '
    # Decker
    ret += 'decker=%s ' % (decker)
    # do flag
    if do_flag:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    ret += 'count=' + str(int(res['NEXP'])) 

    ret += ' foc=' + str(int(focval))
        
    return ret

def getObserved(filename):
    """ getObserved parses a file to find the object names and times
    names, times = getObserved(filename)
    names - list of names, must be first column of file called filename
    times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line

    """
    obs = []
    times = []
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename,level="warn",echo=True)
        return obs, times
    else: 
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '#' or line == "":
                    pass
                else:
                    ls = line.split()
                    obs.append(ls[0])
                    if len(ls) > 15:
                        times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
                    else:
                        times.append(float(ls[1]))
            
    obs.reverse()
    times.reverse()
    return obs, times
	
def calculate_ucsc_exposure_time(vmag, precision, elevation, seeing, bmv, decker="W"):
    """ calculate_ucsc_exposure_time uses the recipe from Burt et al. (2015) to compute the exposure time for a target.

    exp_time, exp_counts, i2counts = calculate_ucsc_exposure_time(vmag, precision, elevation, seeing, bmv, decker="W")
    vmag - numpy array of V magnitudes (Johnson filter, Vega mags)
    precision - required precision for the velocity in m/s
    elevation - elevation of the star above the horizon at the start of the exposure
    seeing - FWHM of the seeing in pixels on the guider
    bmv - (B - V) for the star (both Johnson filters, Vega zeropoint)
    decker - apeture 
    
    exp_time - a numpy array of times in seconds, are integer values
    exp_counts - values for the exposure meter, this can be a floating point value
    i2counts - the required number of median Iodine cell counts, this is calculated from the precision and color of the star, this, in effect, sets the exposure time.


    """
    vmag = np.array(vmag)
    precision = np.array(precision)
    bmv = np.array(bmv)
    precision = np.array(precision)
		
	# Now lets calculate the exposure times
	
	# Desired I2 counts for precision
    i2counts = getI2_K(precision)
    mstars = np.where(bmv > 1.2)
    if len(mstars) > 0:
        i2counts[mstars] = getI2_M(precision[mstars])

    # minimum I2 counts so exposures are not rejected by P. Butler's DRP
    mini2_idx = np.where(i2counts < MIN_I2)
    if len(mini2_idx) > 0:
        i2counts[mini2_idx] = MIN_I2
	
	# Exposure Meter counts to reach desired I2 counts
    exp_counts = getEXPMeter(i2counts, bmv)
    #	exp_counts = 1e9
	# Exposure time to reach desired I2 counts
    exp_time = getEXPTime(i2counts, vmag, bmv, elevation, seeing, decker=decker)
	
    return exp_time, exp_counts, i2counts

def calc_elevations(stars, observer):
    els = []
    for s in stars:
        observer.date = ephem.Date(observer.date)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        els.append(cur_el)
    return np.array(els)

def is_visible(stars, observer, obs_len, min_el, pref_min_el, max_el):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). This is the time frame for which visibility is checked
            min_el: The minimum body elevation to be visible ( degrees )
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    prev_horizon = observer.horizon
    cdate = observer.date
    ret = []
    fin_elevations = []
    start_elevations = []    
    observer.horizon = str(min_el)
    # Now loop over each body to check visibility
    for s, dt in zip(stars, obs_len):
        s.compute()
#        apflog( "is_visible(): Finding limits for %s %s" % (s.ra, s.dec),echo=True)

        # Is the target visible at the end of the observations?
#        apflog( "is_visible(): Is the target visible at the end of the observations?", echo=True)
        observer.date = ephem.Date(cdate + dt/86400.)
        s.compute(observer)
        fin_el = np.degrees(s.alt)
        fin_elevations.append(fin_el)

        # Is the target visible now?
#        apflog( "is_visible(): Is the target visible now?", echo=True)
        observer.date = ephem.Date(cdate)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        start_elevations.append(cur_el)
        
        if fin_el < min_el or fin_el > max_el:
            ret.append(False)
            continue

        if cur_el < min_el or cur_el > max_el:
            ret.append(False)
            continue

        # Does the target remain visible through the observation?
        # The next setting/rising functions throw an exception if the body never sets or rises
        # ex. circumpolar 
        try:
            next_set = observer.next_setting(s)
        except:
            # If it never sets, no need to worry.
            pass
        else:
            # Making the assumption that next_set is a datetime object. Might not be the case
            if next_set < dt:
                # The object will set before the observation finishes
                ret.append(False)
                continue
 #       apflog( "is_visible(): Does the target remain visible through the observation?", echo=True)
        observer.horizon = max_el
        s.compute(observer)

        try:
            next_rise = observer.next_rising(s)
        except:
            # If the body never rises above the max limit no problem
            pass
        else:
            if next_rise < dt:
                # The object rises above the max el before the observation finishes
                ret.append(False)
                continue
        #        apflog( "is_visible(): If the body never rises above the max limit no problem", echo=True)

        observer.horizon = str(pref_min_el)
        s.compute(observer)
        if not s.neverup:
            # will transit above preferred elevation and still rising
            try:
                if ((s.set_time-s.rise_time) > dt/86400.) and cur_el < pref_min_el and np.degrees(s.az) < 180:
                    # this star is currently low on the horizon but will be above the preferred elevation for the requested exposure time
                    # NEED TO HANDLE ONLY RISING STARS
                    ret.append(False)
                    continue
            except:
                pass
        # Everything seems to be fine, so the target is visible!
        ret.append(True)
#	apflog( "is_visible(): done searching targets", echo=True)
    observer.horizon = prev_horizon
    return ret, np.array(start_elevations), np.array(fin_elevations)


def smartList(starlist, time, seeing, slowdown):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    # Convert the unix timestamp into a python datetime

    # punt
    dt = datetime.utcnow()

    if type(time) == float:
        dt = datetime.utcfromtimestamp(int(time))
    elif type(time) == datetime:
        dt = time
    elif type(time) == ephem.Date:
        dt = time.datetime()

    observed, _ = getObserved(os.path.join(os.getcwd(),"observed_targets"))

    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(TARGET_ELEVATION_MIN)
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180.

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Parse the starlist
    try:
        sn, star_table, lines, stars = parseStarlist(starlist)
    except ValueError:
        # This will be raised when the starlist could not be parsed successfully.
        apflog( "No target could be selected because starlist could not be parsed.", level="warn",echo=True)
        return None
    targNum = len(sn)

    # Minimum Brightness based on conditions
    VMAX = 14

    # Distance to stay away from the moon [Between 15 and 25 degrees]

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  

    available = np.ones(targNum, dtype=bool)
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check

    # If seeing is bad, only observe bright targets ( Large VMAG is dim star )       
    brightenough = np.where(star_table[:, DS_VMAG] < VMAX,True,False)
    available = available & brightenough

    obs_length = star_table[:,DS_EXPT] * star_table[:,DS_NSHOTS] + 45 * (star_table[:,DS_NSHOTS]-1)
    vis, star_elevations, fin_els = is_visible(stars,apf_obs,obs_length, TARGET_ELEVATION_MIN, TARGET_ELEVATION_PREF_MIN, TARGET_ELEVATION_MAX)
    available = available & vis
        
    done = [ True if n in observed else False for n in sn ]
    availableandnotdone = available & np.logical_not(done)
    
    if not any(availableandnotdone):
        apflog( "All visible targets have been observed",level="warn",echo=True)
        (good,) = np.where(available)
    else:
        (good,) = np.where(availableandnotdone)

    delta = fin_els[good] - star_elevations[good]
    neg   = np.where(delta < 0)
    pos   = np.where(delta >= 0)
    inv_els = fin_els[good]
    inv_els[neg]  = fin_els[good[neg]] - 90
    inv_els[pos] = 90 - fin_els[good[pos]]
    sort_fin_els_idx = np.argsort(inv_els)
    idx = good[sort_fin_els_idx[0]]
                  
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV'] = 0.6
    res['PRI'] = 10.
    res['SCORE'] = 1.0
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NEXP'] = star_table[idx, DS_NSHOTS]
    res['NAME']   = sn[idx]
    res['SCRIPTOBS'] = lines[idx]
    return res

def format_expmeter(exp_counts, nexp):
    
    exp_counts *= 1.1 
    long_idx = np.where(exp_counts > MAX_EXPMETER, True, False)
    nexp[long_idx] = np.ceil((exp_counts[long_idx]/MAX_EXPMETER) + 1)
    exp_counts[long_idx] = MAX_EXPMETER
    return exp_counts, nexp

def format_time(total, i2counts, hitthemall=False):
    total = np.array(total)
    times = np.zeros(len(total))
    exps  = np.zeros(len(total))

    short_idx = np.where(total < MIN_EXPTIME, True, False)
    times[short_idx] = MIN_EXPTIME  # pad out to make it more likely exposure meter threshold sets actual limit
    exps[short_idx] = [ (np.ceil(MIN_EXPTIME/(t+40)) + 1) for t in total[short_idx] ] 

    middle_idx = np.where((total > MIN_EXPTIME ) &(total < MAX_EXPTIME))
    times[middle_idx] = MAX_EXPTIME # pad out to make it more likely exposure meter threshold sets actual limit
    exps[middle_idx] = 1

    max_idx = np.where(total > MAX_EXPTIME, True, False)
    if hitthemall:
        exps[max_idx] = 1
    else:
        exps[max_idx] = np.ceil(total[max_idx]/MAX_EXPTIME)
    times[max_idx] = MAX_EXPTIME

#    bright_idx = np.where((i2counts > MAX_I2) & (exps == 1), True, False)
#    exps[bright_idx] = [ np.ceil(i/MAX_I2) for i in i2counts[bright_idx] ]
#    times[bright_idx] = np.ceil(total[bright_idx]/exps[bright_idx])

    return times, exps


def getNext(time, seeing, slowdown, star_dates, bstar=False, standardstar=False, verbose=False,sheetn="FakeGoogledex",owner='Vogt',googledex_file="./newgoogledex.csv",method="inquad",observed_file="observed_targets"):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """
    # Convert the unix timestamp into a python datetime
    if type(time) == float:
        dt = datetime.utcfromtimestamp(int(time))
    elif type(time) == datetime:
        dt = time
    elif type(time) == ephem.Date:
        dt = time.datetime()
    else:
        dt = datetime.utcnow()
        # punt

    confg = dict()
    confg['I2'] = 'Y'
    confg['decker']='W'
                
    if verbose:
        apflog( "getNext(): Finding target for time %s" % (dt),echo=True)

    if slowdown > SLOWDOWN_MAX:
        apflog( "getNext(): Slowndown value of %f exceeds maximum of %f at time %s" % (slowdown,SLOWDOWN_MAX,dt),echo=True)
        return None


    try:
        apfguide = ktl.Service('apfguide')
        stamp = apfguide['midptfin'].read(binary=True)
        ptime = datetime.fromtimestamp(stamp)
    except:
        if type(dt) == datetime:
            ptime = dt
        else:
            ptime = datetime.utcnow()

    # List of targets already observed
    observed, obstimes = update_local_googledex(ptime,googledex_file=googledex_file, observed_file=observed_file)
    
    global last_objs_attempted
    try:
        lastline = ktl.read("apftask","SCRIPTOBS_LINE")
        if not bstar:             # otherwise from previous night
            lastobj = lastline.split()[0]

        if verbose:
            apflog( "getNext(): Last object attempted %s" % (lastobj),echo=True)
    except:
        lastobj = None

    if lastobj:
        if lastobj not in observed and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)
            if verbose:
                apflog( "getNext(): Last objects attempted %s" % (last_objs_attempted),echo=True)
            
            if len(last_objs_attempted) > 5:
                apflog( "getNext(): 5 failed acquisition attempts",echo=True)
                last_objs_attempted = []
                return None

    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###
         
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(TARGET_ELEVATION_MIN)
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180. 

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians
    if verbose:
        apflog("getNext(): Parsing the Googledex...",echo=True)
    sn, star_table, do_flag, stars = parseGoogledex(sheetn=sheetn,outfn=googledex_file)
    star_table['apfpri'],star_table['phases'] = getpriority.getpriority(star_table['starname'],star_table,ephem.julian_date(dt),star_dates,apf_obs.sidereal_time(),method=method,standard=standardstar)
    
    targNum = len(sn)
    if verbose:
        apflog("getNext(): Parsed the Googledex...",echo=True)

    # Note which of these are B-Stars for later.
    bstars = star_table['comments'] == 'B star'
    if verbose:
        apflog("getNext(): Finding B stars",echo=True)
        
    available = np.ones(targNum, dtype=bool)
    totexptimes = np.zeros(targNum, dtype=float)
    cur_elevations = np.zeros(targNum, dtype=float)
    i2cnts = np.zeros(targNum, dtype=float)
    star_table['counts'] = np.zeros(targNum, dtype=float)
    star_table['expt'] = np.zeros(targNum, dtype=float)
    star_table['nshots'] = np.zeros(targNum, dtype=int)

    # Is the target behind the moon?
    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  
    moonDist = np.degrees(np.sqrt((moon.ra - star_table['ra'])**2 + (moon.dec - star_table['dec'])**2))
    if verbose:
        apflog("getNext(): Culling stars behind the moon",echo=True)
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check

    # We just need a B star, so restrict our math to those
    if bstar:
        if verbose:
            apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars
        
        f = available
        if verbose:
            apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations = is_visible(fstars, apf_obs, [400]*len(bstars[f]), TARGET_ELEVATION_MIN, TARGET_ELEVATION_PREF_MIN, TARGET_ELEVATION_MAX)
        
        available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]
       		
        star_table['counts'][available] = 1e9
        star_table['expt'][available] = 900
        star_table['nshots'][available] = 2
        totexptimes[available] = 400

    # Just need a normal star for observing
    else:
        # Available and not a BStar
        if verbose:
            apflog("getNext(): Culling B stars",echo=True)
        available = np.logical_and(available, np.logical_not(bstars))
        
        if len(last_objs_attempted)>0:
            for n in last_objs_attempted:
                attempted = (sn == n)
                available = available & np.logical_not(attempted) # Available and not observed

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function
        f = available
        if verbose:
            apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
#        star_elevations=calc_elevations(fstars,apf_obs)
        vis,star_elevations,fin_star_elevations = is_visible(fstars, apf_obs, [0]*len(fstars),  TARGET_ELEVATION_MIN, TARGET_ELEVATION_PREF_MIN, TARGET_ELEVATION_MAX)
        available[f] = available[f] & vis
        f = available
        fstars = [s for s,_ in zip(stars,f) if _ ]
        if verbose:
            apflog("getNext(): Computing exposure times",echo=True)
        exp_times, exp_counts, i2counts = calculate_ucsc_exposure_time( star_table['vmag'][f], star_table['precision'][f], star_elevations[np.array(vis)], seeing, star_table['b-v'][f])
        
        exp_times = exp_times * slowdown
        maxtimes = computeMaxTimes(sn[f],exp_times)
        totexptimes[f] += exp_times
        i2cnts[f] += i2counts
        if verbose:
            apflog("getNext(): Formating exposure times",echo=True)
        star_table['expt'][f], star_table['nshots'][f] = format_time(exp_times,i2counts)
        #        exp_counts /= star_table[f, DS_NSHOTS]
        if verbose:
            apflog("getNext(): Formating exposure meter",echo=True)
        star_table['counts'][f], star_table['nshots'][f] = format_expmeter(exp_counts,star_table['nshots'][f])

        # Is the exposure time too long?
        if verbose:
            apflog("getNext(): Removing really long exposures",echo=True)
        time_check = np.where( exp_times < maxtimes, True, False)
        
        available[f] = available[f] & time_check
        f = available

        # Is the star currently visible?
        if verbose:
            apflog("getNext(): Computing stars visibility",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations = is_visible(fstars, apf_obs, exp_times,  TARGET_ELEVATION_MIN, TARGET_ELEVATION_PREF_MIN, TARGET_ELEVATION_MAX)
        if vis != []:
            available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]

    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

    cadence_check = (ephem.julian_date(dt) - star_table['lastobs']) / star_table['cadence']
    good_cadence = np.where(cadence_check >  1.0, True, False)
    good_cadence_available = available & good_cadence


    if any(good_cadence_available):
        try:
            pri = max(star_table['apfpri'][good_cadence_available])
            sort_i = np.where(star_table['apfpri'][good_cadence_available] == pri, True, False)
            available = good_cadence_available
        except:
            pri = max(star_table['apfpri'][available])
            sort_i = np.where(star_table['apfpri'][available] == pri, True, False)
    elif any(available):
        apflog( "getNext(): No new stars available, going back to the previously observed list.",level="warn",echo=True)
        pri = max(star_table['apfpri'][available])
        sort_i = np.where(star_table['apfpri'][available] == pri, True, False)
    else:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None
    
    starstr = "getNext(): star table available: %s" % (sn[available][sort_i]) 
    apflog(starstr,echo=True)

    sort_j = cur_elevations[available][sort_i].argsort()[::-1]
    
    t_n = sn[available][sort_i][sort_j][0]

    elstr= "getNext(): star elevations %s" % (cur_elevations[available][sort_i][sort_j])
    apflog(elstr,echo=True)

    t_n = sn[available][sort_i][sort_j][0]

    apflog("getNext(): selected target %s" %( t_n) )

    idx, = np.where(sn == t_n)
    idx = idx[0]

    stars[idx].compute(apf_obs)
    
    res = dict()
        
    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = 0.0
    res['PM_DEC'] = 0.0
    res['VMAG']   = star_table['vmag'][idx]
    res['BV']     = star_table['b-v'][idx]
    res['COUNTS'] = star_table['counts'][idx]
    res['EXP_TIME'] = star_table['expt'][idx]
    res['NEXP'] = star_table['nshots'][idx]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['I2CNTS'] = i2cnts[idx]
    res['NAME']   = sn[idx]
    res['SCORE']  = star_table['apfpri'][idx]
    res['PRI']    = star_table['apfpri'][idx]
    res['SCRIPTOBS'] = makeScriptobsLine(res, do_flag[idx], dt, decker=confg['decker'], I2=confg['I2'], owner=owner)
    return res




if __name__ == '__main__':
    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 13.99, 5, 5.0, bstar=True, verbose=True)
    ot.write("%s\n" % (result["SCRIPTOBS"]))
    ot.close()
    starttime += 400
    for i in range(5):
        
        result = getNext(starttime, 13.99, 0.4, 5.0, bstar=False, verbose=True)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print ("Get None target")
        else:
            for k in result:
                print (k, result[k])
        ot = open(otfn,"a")
        ot.write("%s\n" % (result["SCRIPTOBS"]))
        ot.close()
        starttime += result["EXP_TIME"]
                
    print ("Done")
