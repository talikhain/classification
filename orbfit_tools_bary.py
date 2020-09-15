from __future__ import division
import pyOrbfit as orbfit
import ephem
import numpy as np

'''
Some tools for working with abg files instead of directly from observations.
   orbfit_abg(abgfile) -- reads an abgfile into the appropriate data structures, computes aei orbit
   predict_from_abg(abginfo, date, obscode=807) -- predicts position on a given date.
   orbfit_aei(aeifile) -- reads an aei file into simple Python data structures (but won't interact with orbfit)

   Example usage:
       abginfo = orbfit_abg('QO441.abg')
       pos_pred = predict_from_abg(abginfo, ephem.date('2015/01/01'))
       print ephem.hours(pos['ra']), ephem.degrees(pos['dec']), pos['err']

'''

GM = 4.*np.pi*np.pi/1.0000378 # solar gravitation, from orbfit
combinedMass = GM * 1.00134  #Alter GM to account for total SS mass, from orbfit
DAY = 1/365.25

def orbfit_abg(abgfile):
    with open(abgfile) as f:
        lines = [line for line in f if line[0]!='#']
        elem_abg = [float(e) for e in lines[0].split()]
        cov_arr = np.array([[float(de) for de in lines[i].split()] for i in range(1,7)],dtype=np.float64)
        last = lines[7].split()
        lat0 = float(last[0])
        lon0 = float(last[1])
        xBary = float(last[2])
        yBary = float(last[3])
        zBary = float(last[4])
        jd0 = float(last[5])
        abginfo = {'elem':elem_abg, 'cov':cov_arr, 'lat0':lat0, 'lon0':lon0, 
                   'xBary':xBary, 'yBary':yBary, 'zBary':zBary, 'jd0':jd0}
        jd0 = abginfo['jd0']     # zeropoint of time scale
        # Create space for various useful matrices
        cov_abg = orbfit.dmatrix(1,6,1,6)
        p_in = orbfit.PBASIS()
        orbit_xyz = orbfit.XVBASIS()
        orbit_aei = orbfit.ORBIT()
        p_in.a = abginfo['elem'][0]
        p_in.adot = abginfo['elem'][1]
        p_in.b = abginfo['elem'][2]
        p_in.bdot = abginfo['elem'][3]
        p_in.g = abginfo['elem'][4]
        p_in.gdot = abginfo['elem'][5]
        # Here follows kludginess to put the covariance matrix into C double pointers
        c = orbfit.doubleArray(36)
        cc = abginfo['cov'].reshape(36, order='F')
        for i in range(len(cc)):
            c[i] = cc[i]
        orbfit.unflatten_cov(c, 6, cov_abg)
        abginfo['pbasis'] = p_in
        abginfo['cov_abg']  = cov_abg
        cov_xyz = orbfit.dmatrix(1,6,1,6)
        cov_aei = orbfit.dmatrix(1,6,1,6)
        derivs  = orbfit.dmatrix(1,6,1,6)
    #   Transform the orbit basis and get the deriv. matrix 
        orbfit.cvar.jd0 = abginfo['jd0']
        orbfit.cvar.xBary = abginfo['xBary']
        orbfit.cvar.yBary = abginfo['yBary']
        orbfit.cvar.zBary = abginfo['zBary']
        orbfit.cvar.lon0 = abginfo['lon0']*np.pi/180
        orbfit.cvar.lat0 = abginfo['lat0']*np.pi/180 
        #CHANGE HERE
        orbfit.pbasis_to_bary(p_in, orbit_xyz, derivs) 
        orbfit.orbitElements(orbit_xyz, orbit_aei) # aei orbit elements 
        
        abginfo['orbit_aei']=orbit_aei
        abginfo['cov_aei']=cov_aei
        abginfo['elems_aei'] = {'a':orbit_aei.a, 'e':orbit_aei.e, 'i':orbit_aei.i, 'lan':orbit_aei.lan,
                               'aop':orbit_aei.aop, 'top':orbit_aei.T}
    return abginfo

def predict_from_abg(abginfo, date, obscode=807):

    p_in = abginfo['pbasis']
    cov_abg = abginfo['cov_abg']
    orbfit.cvar.jd0 = abginfo['jd0']
    orbfit.cvar.xBary = abginfo['xBary']
    orbfit.cvar.yBary = abginfo['yBary']
    orbfit.cvar.zBary = abginfo['zBary']
    orbfit.cvar.lon0 = abginfo['lon0']*np.pi/180
    orbfit.cvar.lat0 = abginfo['lat0']*np.pi/180
    sigxy = orbfit.dmatrix(1,2,1,2)
    derivs = orbfit.dmatrix(1,2,1,2)
    covecl = orbfit.dmatrix(1,2,1,2)
    coveq = orbfit.dmatrix(1,2,1,2)
    # Fill the OBSERVATION structure
    futobs = orbfit.OBSERVATION()
    futobs.obscode = obscode
    futobs.obstime = (ephem.julian_date(date)-orbfit.cvar.jd0)*orbfit.DAY
    futobs.xe = -999  # force evaluation of earth3D
    
    dx = orbfit.dvector(1,6)
    dy = orbfit.dvector(1,6)

    thetax, thetay = orbfit.kbo2d(p_in, futobs, dx, dy)

    # Predicted position, in abg basis:
    orbfit.predict_posn(p_in, cov_abg, futobs, sigxy)
    solar_elongation = orbfit.elongation(futobs)/orbfit.DTOR       # solar elongation in degrees
    opposition_angle = orbfit.opposition_angle(futobs)/orbfit.DTOR # opposition angle in degrees
    lat_ec, lon_ec = orbfit.proj_to_ec(futobs.thetax, futobs.thetay, orbfit.cvar.lat0, orbfit.cvar.lon0, derivs)  # project to ecliptic coords
    orbfit.covar_map(sigxy, derivs, covecl, 2, 2)    # map the covariance
    ra_eq, dec_eq = orbfit.ec_to_eq(lat_ec, lon_ec, derivs)    # transform to ICRS to compute ra, dec
    eq = ephem.Equatorial(ephem.Ecliptic(lat_ec, lon_ec))
    ra_eq2 = eq.ra*180/np.pi
    dec_eq2 = eq.dec*180/np.pi
    if ra_eq<0: ra_eq += 2*np.pi        # convert angle to (0,2pi) range
#        print ra_eq*180/np.pi, dec_eq*180/np.pi, ra_eq2, dec_eq2, ra_eq*180/np.pi-ra_eq2, dec_eq*180/np.pi-dec_eq2
    orbfit.covar_map(covecl, derivs, coveq, 2, 2)    # map the covariance
#        print ephem.hours(ra_eq), ephem.degrees(dec_eq)
    # Compute ICRS error ellipse
    c2d = orbfit.doubleArray(4)   # stoopid workaround for double pointers...
    orbfit.flatten_cov(coveq, 2, c2d)
    covar_eq = np.array([c2d[i] for i in range(4)]).reshape(2,2)
    xx = covar_eq[0][0]*np.cos(dec_eq)**2
    xy = covar_eq[0][1]*np.cos(dec_eq)
    yy = covar_eq[1][1]
    pos_angle = 0.5*np.arctan2(2.*xy,(xx-yy)) * 180./np.pi
    pos_angle = 90 - pos_angle     # convert to astronomy convention of measuring position angle North through East
    bovasqrd  = (xx + yy - np.sqrt((xx-yy)**2 + (2*xy)**2)) / (xx + yy + np.sqrt((xx-yy)**2 + (2*xy)**2))
    det = xx*yy-xy*xy
    a = (det/bovasqrd)**(1/4)/orbfit.ARCSEC   # semimajor, minor axes of error ellipse, in arcsec
    b = (det*bovasqrd)**(1/4)/orbfit.ARCSEC
    err_ellipse = dict(a=a, b=b, PA=pos_angle)   # store as a dictionary
    
    pos = dict(ra=ephem.hours(ra_eq), dec=ephem.degrees(dec_eq), err=err_ellipse, elong=solar_elongation, opp=opposition_angle)
    return pos

def djd(jd):
    '''Dublin Julian date from JD
    '''
    return jd-2415020

def mean_anomaly(elems, jd0):
    if np.isnan(elems['top']):
        mu = -999
    else:
        #CHANGE HERE combinedMass to GM
        mu = (jd0-elems['top'])*np.sqrt(combinedMass/elems['a']**3)*DAY*180/np.pi
        if mu<0: mu+=360
    return mu
    
def orbfit_aei(aeifile):
    with open(aeifile) as f:
        lines = [line for line in f if line[0]!='#' or 'epoch' in line]
        jd0 = float(lines[0].split('epoch')[1][:-2])
        elems = [float(e) for e in lines[1].split()]
        elem_aei = {'a':elems[0], 'e':elems[1], 'i':elems[2], 'lan':elems[3], 'aop':elems[4],'top':elems[5]}
        elem_aei_err = [float(e) for e in lines[2].split() if e != '+-']
        cov_arr = np.array([[float(de) for de in lines[i].split()] for i in range(3,9)],dtype=np.float64)
        name = aeifile.split('.aei')[0]
        aeiinfo = {'elem':elem_aei, 'cov':cov_arr, 'jd0':jd0, 'mu':mean_anomaly(elem_aei, jd0),
                  'name':name}
    return aeiinfo
 