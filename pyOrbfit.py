# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_pyOrbfit')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_pyOrbfit')
    _pyOrbfit = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pyOrbfit', [dirname(__file__)])
        except ImportError:
            import _pyOrbfit
            return _pyOrbfit
        try:
            _mod = imp.load_module('_pyOrbfit', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _pyOrbfit = swig_import_helper()
    del swig_import_helper
else:
    import _pyOrbfit
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class OBSERVATION_ARRAY(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBSERVATION_ARRAY, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBSERVATION_ARRAY, name)
    __repr__ = _swig_repr

    def __init__(self, nelements):
        this = _pyOrbfit.new_OBSERVATION_ARRAY(nelements)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_OBSERVATION_ARRAY
    __del__ = lambda self: None

    def __getitem__(self, index):
        return _pyOrbfit.OBSERVATION_ARRAY___getitem__(self, index)

    def __setitem__(self, index, value):
        return _pyOrbfit.OBSERVATION_ARRAY___setitem__(self, index, value)

    def cast(self):
        return _pyOrbfit.OBSERVATION_ARRAY_cast(self)
    if _newclass:
        frompointer = staticmethod(_pyOrbfit.OBSERVATION_ARRAY_frompointer)
    else:
        frompointer = _pyOrbfit.OBSERVATION_ARRAY_frompointer
OBSERVATION_ARRAY_swigregister = _pyOrbfit.OBSERVATION_ARRAY_swigregister
OBSERVATION_ARRAY_swigregister(OBSERVATION_ARRAY)

def OBSERVATION_ARRAY_frompointer(t):
    return _pyOrbfit.OBSERVATION_ARRAY_frompointer(t)
OBSERVATION_ARRAY_frompointer = _pyOrbfit.OBSERVATION_ARRAY_frompointer

class doubleArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, doubleArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, doubleArray, name)
    __repr__ = _swig_repr

    def __init__(self, nelements):
        this = _pyOrbfit.new_doubleArray(nelements)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_doubleArray
    __del__ = lambda self: None

    def __getitem__(self, index):
        return _pyOrbfit.doubleArray___getitem__(self, index)

    def __setitem__(self, index, value):
        return _pyOrbfit.doubleArray___setitem__(self, index, value)

    def cast(self):
        return _pyOrbfit.doubleArray_cast(self)
    if _newclass:
        frompointer = staticmethod(_pyOrbfit.doubleArray_frompointer)
    else:
        frompointer = _pyOrbfit.doubleArray_frompointer
doubleArray_swigregister = _pyOrbfit.doubleArray_swigregister
doubleArray_swigregister(doubleArray)

def doubleArray_frompointer(t):
    return _pyOrbfit.doubleArray_frompointer(t)
doubleArray_frompointer = _pyOrbfit.doubleArray_frompointer

MAXOBS = _pyOrbfit.MAXOBS
DEFAULT_EPHEM_FILE = _pyOrbfit.DEFAULT_EPHEM_FILE
EPHEM_ENVIRON = _pyOrbfit.EPHEM_ENVIRON
DEFAULT_OBSERVATORY_FILE = _pyOrbfit.DEFAULT_OBSERVATORY_FILE
OBS_ENVIRON = _pyOrbfit.OBS_ENVIRON
MAX_SITES = _pyOrbfit.MAX_SITES
OBSCODE_SSBARY = _pyOrbfit.OBSCODE_SSBARY
OBSCODE_GEOCENTER = _pyOrbfit.OBSCODE_GEOCENTER
OBSCODE_ORBITAL = _pyOrbfit.OBSCODE_ORBITAL
DEFAULT_DTHETA = _pyOrbfit.DEFAULT_DTHETA
PI = _pyOrbfit.PI
TPI = _pyOrbfit.TPI
DTOR = _pyOrbfit.DTOR
GM = _pyOrbfit.GM
SSMASS = _pyOrbfit.SSMASS
ARCSEC = _pyOrbfit.ARCSEC
DAY = _pyOrbfit.DAY
SPEED_OF_LIGHT = _pyOrbfit.SPEED_OF_LIGHT
ECL = _pyOrbfit.ECL
TSTEP = _pyOrbfit.TSTEP
EARTHMASS = _pyOrbfit.EARTHMASS
class PBASIS(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PBASIS, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PBASIS, name)
    __repr__ = _swig_repr
    __swig_setmethods__["a"] = _pyOrbfit.PBASIS_a_set
    __swig_getmethods__["a"] = _pyOrbfit.PBASIS_a_get
    if _newclass:
        a = _swig_property(_pyOrbfit.PBASIS_a_get, _pyOrbfit.PBASIS_a_set)
    __swig_setmethods__["adot"] = _pyOrbfit.PBASIS_adot_set
    __swig_getmethods__["adot"] = _pyOrbfit.PBASIS_adot_get
    if _newclass:
        adot = _swig_property(_pyOrbfit.PBASIS_adot_get, _pyOrbfit.PBASIS_adot_set)
    __swig_setmethods__["b"] = _pyOrbfit.PBASIS_b_set
    __swig_getmethods__["b"] = _pyOrbfit.PBASIS_b_get
    if _newclass:
        b = _swig_property(_pyOrbfit.PBASIS_b_get, _pyOrbfit.PBASIS_b_set)
    __swig_setmethods__["bdot"] = _pyOrbfit.PBASIS_bdot_set
    __swig_getmethods__["bdot"] = _pyOrbfit.PBASIS_bdot_get
    if _newclass:
        bdot = _swig_property(_pyOrbfit.PBASIS_bdot_get, _pyOrbfit.PBASIS_bdot_set)
    __swig_setmethods__["g"] = _pyOrbfit.PBASIS_g_set
    __swig_getmethods__["g"] = _pyOrbfit.PBASIS_g_get
    if _newclass:
        g = _swig_property(_pyOrbfit.PBASIS_g_get, _pyOrbfit.PBASIS_g_set)
    __swig_setmethods__["gdot"] = _pyOrbfit.PBASIS_gdot_set
    __swig_getmethods__["gdot"] = _pyOrbfit.PBASIS_gdot_get
    if _newclass:
        gdot = _swig_property(_pyOrbfit.PBASIS_gdot_get, _pyOrbfit.PBASIS_gdot_set)

    def __init__(self):
        this = _pyOrbfit.new_PBASIS()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_PBASIS
    __del__ = lambda self: None
PBASIS_swigregister = _pyOrbfit.PBASIS_swigregister
PBASIS_swigregister(PBASIS)

class XVBASIS(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, XVBASIS, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, XVBASIS, name)
    __repr__ = _swig_repr
    __swig_setmethods__["x"] = _pyOrbfit.XVBASIS_x_set
    __swig_getmethods__["x"] = _pyOrbfit.XVBASIS_x_get
    if _newclass:
        x = _swig_property(_pyOrbfit.XVBASIS_x_get, _pyOrbfit.XVBASIS_x_set)
    __swig_setmethods__["y"] = _pyOrbfit.XVBASIS_y_set
    __swig_getmethods__["y"] = _pyOrbfit.XVBASIS_y_get
    if _newclass:
        y = _swig_property(_pyOrbfit.XVBASIS_y_get, _pyOrbfit.XVBASIS_y_set)
    __swig_setmethods__["z"] = _pyOrbfit.XVBASIS_z_set
    __swig_getmethods__["z"] = _pyOrbfit.XVBASIS_z_get
    if _newclass:
        z = _swig_property(_pyOrbfit.XVBASIS_z_get, _pyOrbfit.XVBASIS_z_set)
    __swig_setmethods__["xdot"] = _pyOrbfit.XVBASIS_xdot_set
    __swig_getmethods__["xdot"] = _pyOrbfit.XVBASIS_xdot_get
    if _newclass:
        xdot = _swig_property(_pyOrbfit.XVBASIS_xdot_get, _pyOrbfit.XVBASIS_xdot_set)
    __swig_setmethods__["ydot"] = _pyOrbfit.XVBASIS_ydot_set
    __swig_getmethods__["ydot"] = _pyOrbfit.XVBASIS_ydot_get
    if _newclass:
        ydot = _swig_property(_pyOrbfit.XVBASIS_ydot_get, _pyOrbfit.XVBASIS_ydot_set)
    __swig_setmethods__["zdot"] = _pyOrbfit.XVBASIS_zdot_set
    __swig_getmethods__["zdot"] = _pyOrbfit.XVBASIS_zdot_get
    if _newclass:
        zdot = _swig_property(_pyOrbfit.XVBASIS_zdot_get, _pyOrbfit.XVBASIS_zdot_set)
    __swig_setmethods__["jd0"] = _pyOrbfit.XVBASIS_jd0_set
    __swig_getmethods__["jd0"] = _pyOrbfit.XVBASIS_jd0_get
    if _newclass:
        jd0 = _swig_property(_pyOrbfit.XVBASIS_jd0_get, _pyOrbfit.XVBASIS_jd0_set)

    def __init__(self):
        this = _pyOrbfit.new_XVBASIS()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_XVBASIS
    __del__ = lambda self: None
XVBASIS_swigregister = _pyOrbfit.XVBASIS_swigregister
XVBASIS_swigregister(XVBASIS)

class ORBIT(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ORBIT, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ORBIT, name)
    __repr__ = _swig_repr
    __swig_setmethods__["a"] = _pyOrbfit.ORBIT_a_set
    __swig_getmethods__["a"] = _pyOrbfit.ORBIT_a_get
    if _newclass:
        a = _swig_property(_pyOrbfit.ORBIT_a_get, _pyOrbfit.ORBIT_a_set)
    __swig_setmethods__["e"] = _pyOrbfit.ORBIT_e_set
    __swig_getmethods__["e"] = _pyOrbfit.ORBIT_e_get
    if _newclass:
        e = _swig_property(_pyOrbfit.ORBIT_e_get, _pyOrbfit.ORBIT_e_set)
    __swig_setmethods__["i"] = _pyOrbfit.ORBIT_i_set
    __swig_getmethods__["i"] = _pyOrbfit.ORBIT_i_get
    if _newclass:
        i = _swig_property(_pyOrbfit.ORBIT_i_get, _pyOrbfit.ORBIT_i_set)
    __swig_setmethods__["lan"] = _pyOrbfit.ORBIT_lan_set
    __swig_getmethods__["lan"] = _pyOrbfit.ORBIT_lan_get
    if _newclass:
        lan = _swig_property(_pyOrbfit.ORBIT_lan_get, _pyOrbfit.ORBIT_lan_set)
    __swig_setmethods__["aop"] = _pyOrbfit.ORBIT_aop_set
    __swig_getmethods__["aop"] = _pyOrbfit.ORBIT_aop_get
    if _newclass:
        aop = _swig_property(_pyOrbfit.ORBIT_aop_get, _pyOrbfit.ORBIT_aop_set)
    __swig_setmethods__["T"] = _pyOrbfit.ORBIT_T_set
    __swig_getmethods__["T"] = _pyOrbfit.ORBIT_T_get
    if _newclass:
        T = _swig_property(_pyOrbfit.ORBIT_T_get, _pyOrbfit.ORBIT_T_set)
    __swig_setmethods__["ma"] = _pyOrbfit.ORBIT_ma_set
    __swig_getmethods__["ma"] = _pyOrbfit.ORBIT_ma_get
    if _newclass:
        ma = _swig_property(_pyOrbfit.ORBIT_ma_get, _pyOrbfit.ORBIT_ma_set)

    def __init__(self):
        this = _pyOrbfit.new_ORBIT()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_ORBIT
    __del__ = lambda self: None
ORBIT_swigregister = _pyOrbfit.ORBIT_swigregister
ORBIT_swigregister(ORBIT)

class date_time(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, date_time, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, date_time, name)
    __repr__ = _swig_repr
    __swig_setmethods__["y"] = _pyOrbfit.date_time_y_set
    __swig_getmethods__["y"] = _pyOrbfit.date_time_y_get
    if _newclass:
        y = _swig_property(_pyOrbfit.date_time_y_get, _pyOrbfit.date_time_y_set)
    __swig_setmethods__["mo"] = _pyOrbfit.date_time_mo_set
    __swig_getmethods__["mo"] = _pyOrbfit.date_time_mo_get
    if _newclass:
        mo = _swig_property(_pyOrbfit.date_time_mo_get, _pyOrbfit.date_time_mo_set)
    __swig_setmethods__["d"] = _pyOrbfit.date_time_d_set
    __swig_getmethods__["d"] = _pyOrbfit.date_time_d_get
    if _newclass:
        d = _swig_property(_pyOrbfit.date_time_d_get, _pyOrbfit.date_time_d_set)
    __swig_setmethods__["h"] = _pyOrbfit.date_time_h_set
    __swig_getmethods__["h"] = _pyOrbfit.date_time_h_get
    if _newclass:
        h = _swig_property(_pyOrbfit.date_time_h_get, _pyOrbfit.date_time_h_set)
    __swig_setmethods__["mn"] = _pyOrbfit.date_time_mn_set
    __swig_getmethods__["mn"] = _pyOrbfit.date_time_mn_get
    if _newclass:
        mn = _swig_property(_pyOrbfit.date_time_mn_get, _pyOrbfit.date_time_mn_set)
    __swig_setmethods__["s"] = _pyOrbfit.date_time_s_set
    __swig_getmethods__["s"] = _pyOrbfit.date_time_s_get
    if _newclass:
        s = _swig_property(_pyOrbfit.date_time_s_get, _pyOrbfit.date_time_s_set)

    def __init__(self):
        this = _pyOrbfit.new_date_time()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_date_time
    __del__ = lambda self: None
date_time_swigregister = _pyOrbfit.date_time_swigregister
date_time_swigregister(date_time)

class OBSERVATION(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBSERVATION, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBSERVATION, name)
    __repr__ = _swig_repr
    __swig_setmethods__["thetax"] = _pyOrbfit.OBSERVATION_thetax_set
    __swig_getmethods__["thetax"] = _pyOrbfit.OBSERVATION_thetax_get
    if _newclass:
        thetax = _swig_property(_pyOrbfit.OBSERVATION_thetax_get, _pyOrbfit.OBSERVATION_thetax_set)
    __swig_setmethods__["dthetax"] = _pyOrbfit.OBSERVATION_dthetax_set
    __swig_getmethods__["dthetax"] = _pyOrbfit.OBSERVATION_dthetax_get
    if _newclass:
        dthetax = _swig_property(_pyOrbfit.OBSERVATION_dthetax_get, _pyOrbfit.OBSERVATION_dthetax_set)
    __swig_setmethods__["thetay"] = _pyOrbfit.OBSERVATION_thetay_set
    __swig_getmethods__["thetay"] = _pyOrbfit.OBSERVATION_thetay_get
    if _newclass:
        thetay = _swig_property(_pyOrbfit.OBSERVATION_thetay_get, _pyOrbfit.OBSERVATION_thetay_set)
    __swig_setmethods__["dthetay"] = _pyOrbfit.OBSERVATION_dthetay_set
    __swig_getmethods__["dthetay"] = _pyOrbfit.OBSERVATION_dthetay_get
    if _newclass:
        dthetay = _swig_property(_pyOrbfit.OBSERVATION_dthetay_get, _pyOrbfit.OBSERVATION_dthetay_set)
    __swig_setmethods__["obstime"] = _pyOrbfit.OBSERVATION_obstime_set
    __swig_getmethods__["obstime"] = _pyOrbfit.OBSERVATION_obstime_get
    if _newclass:
        obstime = _swig_property(_pyOrbfit.OBSERVATION_obstime_get, _pyOrbfit.OBSERVATION_obstime_set)
    __swig_setmethods__["obscode"] = _pyOrbfit.OBSERVATION_obscode_set
    __swig_getmethods__["obscode"] = _pyOrbfit.OBSERVATION_obscode_get
    if _newclass:
        obscode = _swig_property(_pyOrbfit.OBSERVATION_obscode_get, _pyOrbfit.OBSERVATION_obscode_set)
    __swig_setmethods__["xe"] = _pyOrbfit.OBSERVATION_xe_set
    __swig_getmethods__["xe"] = _pyOrbfit.OBSERVATION_xe_get
    if _newclass:
        xe = _swig_property(_pyOrbfit.OBSERVATION_xe_get, _pyOrbfit.OBSERVATION_xe_set)
    __swig_setmethods__["ye"] = _pyOrbfit.OBSERVATION_ye_set
    __swig_getmethods__["ye"] = _pyOrbfit.OBSERVATION_ye_get
    if _newclass:
        ye = _swig_property(_pyOrbfit.OBSERVATION_ye_get, _pyOrbfit.OBSERVATION_ye_set)
    __swig_setmethods__["ze"] = _pyOrbfit.OBSERVATION_ze_set
    __swig_getmethods__["ze"] = _pyOrbfit.OBSERVATION_ze_get
    if _newclass:
        ze = _swig_property(_pyOrbfit.OBSERVATION_ze_get, _pyOrbfit.OBSERVATION_ze_set)

    def __init__(self):
        this = _pyOrbfit.new_OBSERVATION()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyOrbfit.delete_OBSERVATION
    __del__ = lambda self: None
OBSERVATION_swigregister = _pyOrbfit.OBSERVATION_swigregister
OBSERVATION_swigregister(OBSERVATION)


def kbo3d(pin, t, xout, yout, dx, dy, dz):
    return _pyOrbfit.kbo3d(pin, t, xout, yout, dx, dy, dz)
kbo3d = _pyOrbfit.kbo3d

def earth3d(t, obscode, z):
    return _pyOrbfit.earth3d(t, obscode, z)
earth3d = _pyOrbfit.earth3d

def earth_ssbary(jd, obscode, z):
    return _pyOrbfit.earth_ssbary(jd, obscode, z)
earth_ssbary = _pyOrbfit.earth_ssbary

def body3d(t, body, z, vxyz):
    return _pyOrbfit.body3d(t, body, z, vxyz)
body3d = _pyOrbfit.body3d

def bodycenter_ssbary(jd, xyz, body, vxyz):
    return _pyOrbfit.bodycenter_ssbary(jd, xyz, body, vxyz)
bodycenter_ssbary = _pyOrbfit.bodycenter_ssbary

def geocenter_ssbary(jd, xyz):
    return _pyOrbfit.geocenter_ssbary(jd, xyz)
geocenter_ssbary = _pyOrbfit.geocenter_ssbary

def observatory_geocenter(jd, obscode, xobs, yobs, zobs):
    return _pyOrbfit.observatory_geocenter(jd, obscode, xobs, yobs, zobs)
observatory_geocenter = _pyOrbfit.observatory_geocenter

def read_observatories(fname):
    return _pyOrbfit.read_observatories(fname)
read_observatories = _pyOrbfit.read_observatories

def print_help(h):
    return _pyOrbfit.print_help(h)
print_help = _pyOrbfit.print_help

def zenith_horizon(obscode):
    return _pyOrbfit.zenith_horizon(obscode)
zenith_horizon = _pyOrbfit.zenith_horizon

def zenith_angle(o):
    return _pyOrbfit.zenith_angle(o)
zenith_angle = _pyOrbfit.zenith_angle

def is_visible(o):
    return _pyOrbfit.is_visible(o)
is_visible = _pyOrbfit.is_visible

def kbo2d(pin, obs, dx, dy):
    return _pyOrbfit.kbo2d(pin, obs, dx, dy)
kbo2d = _pyOrbfit.kbo2d

def kbo2d_linear(pin, obs, dx, dy):
    return _pyOrbfit.kbo2d_linear(pin, obs, dx, dy)
kbo2d_linear = _pyOrbfit.kbo2d_linear

def scan_observation(inbuff, obs):
    return _pyOrbfit.scan_observation(inbuff, obs)
scan_observation = _pyOrbfit.scan_observation

def fit_radec(fname, pbasis, orb):
    return _pyOrbfit.fit_radec(fname, pbasis, orb)
fit_radec = _pyOrbfit.fit_radec

def set_mpc_dtheta(d):
    return _pyOrbfit.set_mpc_dtheta(d)
set_mpc_dtheta = _pyOrbfit.set_mpc_dtheta

def set_ephem_file(fname):
    return _pyOrbfit.set_ephem_file(fname)
set_ephem_file = _pyOrbfit.set_ephem_file

def set_observatory_file(fname):
    return _pyOrbfit.set_observatory_file(fname)
set_observatory_file = _pyOrbfit.set_observatory_file

def read_options(iarg, argc, argv):
    return _pyOrbfit.read_options(iarg, argc, argv)
read_options = _pyOrbfit.read_options

def prelim_fit(obsarray, nobs, pout, covar):
    return _pyOrbfit.prelim_fit(obsarray, nobs, pout, covar)
prelim_fit = _pyOrbfit.prelim_fit

def predict_posn(pin, covar, obs, sigxy):
    return _pyOrbfit.predict_posn(pin, covar, obs, sigxy)
predict_posn = _pyOrbfit.predict_posn

def print_matrix(fptr, matrix, xdim, ydim):
    return _pyOrbfit.print_matrix(fptr, matrix, xdim, ydim)
print_matrix = _pyOrbfit.print_matrix

def fgets_nocomment(buffer, length, fptr, fpout):
    return _pyOrbfit.fgets_nocomment(buffer, length, fptr, fpout)
fgets_nocomment = _pyOrbfit.fgets_nocomment

def eq_to_ec(ra_eq, dec_eq, partials):
    return _pyOrbfit.eq_to_ec(ra_eq, dec_eq, partials)
eq_to_ec = _pyOrbfit.eq_to_ec

def xyz_eq_to_ec(x_eq, y_eq, z_eq, x_ec, y_ec, z_ec, partials):
    return _pyOrbfit.xyz_eq_to_ec(x_eq, y_eq, z_eq, x_ec, y_ec, z_ec, partials)
xyz_eq_to_ec = _pyOrbfit.xyz_eq_to_ec

def ec_to_eq(lat_ec, lon_ec, partials):
    return _pyOrbfit.ec_to_eq(lat_ec, lon_ec, partials)
ec_to_eq = _pyOrbfit.ec_to_eq

def xyz_ec_to_eq(x_ec, y_ec, z_ec, x_eq, y_eq, z_eq, partials):
    return _pyOrbfit.xyz_ec_to_eq(x_ec, y_ec, z_ec, x_eq, y_eq, z_eq, partials)
xyz_ec_to_eq = _pyOrbfit.xyz_ec_to_eq

def ec_to_proj(lat_ec, lon_ec, x_proj, y_proj, lat0, lon0, partials):
    return _pyOrbfit.ec_to_proj(lat_ec, lon_ec, x_proj, y_proj, lat0, lon0, partials)
ec_to_proj = _pyOrbfit.ec_to_proj

def proj_to_ec(x_proj, y_proj, lat0, lon0, partials):
    return _pyOrbfit.proj_to_ec(x_proj, y_proj, lat0, lon0, partials)
proj_to_ec = _pyOrbfit.proj_to_ec

def xyz_ec_to_proj(x_ec, y_ec, z_ec, x_p, y_p, z_p, lat0, lon0, partials):
    return _pyOrbfit.xyz_ec_to_proj(x_ec, y_ec, z_ec, x_p, y_p, z_p, lat0, lon0, partials)
xyz_ec_to_proj = _pyOrbfit.xyz_ec_to_proj

def xyz_proj_to_ec(x_p, y_p, z_p, x_ec, y_ec, z_ec, lat0, lon0, partials):
    return _pyOrbfit.xyz_proj_to_ec(x_p, y_p, z_p, x_ec, y_ec, z_ec, lat0, lon0, partials)
xyz_proj_to_ec = _pyOrbfit.xyz_proj_to_ec

def pbasis_to_bary(p, xv, partials):
    return _pyOrbfit.pbasis_to_bary(p, xv, partials)
pbasis_to_bary = _pyOrbfit.pbasis_to_bary

def matrix_multiply(m1, m2, mout, x1, y1, x2, y2):
    return _pyOrbfit.matrix_multiply(m1, m2, mout, x1, y1, x2, y2)
matrix_multiply = _pyOrbfit.matrix_multiply

def orbitElements(xv, orb):
    return _pyOrbfit.orbitElements(xv, orb)
orbitElements = _pyOrbfit.orbitElements

def elements_to_xv(o, jd, xv):
    return _pyOrbfit.elements_to_xv(o, jd, xv)
elements_to_xv = _pyOrbfit.elements_to_xv

def elements_to_pbasis(o, jd, obscode, p):
    return _pyOrbfit.elements_to_pbasis(o, jd, obscode, p)
elements_to_pbasis = _pyOrbfit.elements_to_pbasis

def covar_map(covar_in, derivs, covar_out, kin, kout):
    return _pyOrbfit.covar_map(covar_in, derivs, covar_out, kin, kout)
covar_map = _pyOrbfit.covar_map

def date_to_jd(date):
    return _pyOrbfit.date_to_jd(date)
date_to_jd = _pyOrbfit.date_to_jd

def aei_derivs(xv, daei_dxv):
    return _pyOrbfit.aei_derivs(xv, daei_dxv)
aei_derivs = _pyOrbfit.aei_derivs

def aei_derivs_helio(xv, daei_dxv):
    return _pyOrbfit.aei_derivs_helio(xv, daei_dxv)
aei_derivs_helio = _pyOrbfit.aei_derivs_helio

def elongation(obs):
    return _pyOrbfit.elongation(obs)
elongation = _pyOrbfit.elongation

def opposition_angle(obs):
    return _pyOrbfit.opposition_angle(obs)
opposition_angle = _pyOrbfit.opposition_angle

def fake_observation(p, obs):
    return _pyOrbfit.fake_observation(p, obs)
fake_observation = _pyOrbfit.fake_observation

def flatten_cov(cov, ndim, cov1d):
    return _pyOrbfit.flatten_cov(cov, ndim, cov1d)
flatten_cov = _pyOrbfit.flatten_cov

def unflatten_cov(cov1d, ndim, cov):
    return _pyOrbfit.unflatten_cov(cov1d, ndim, cov)
unflatten_cov = _pyOrbfit.unflatten_cov

def dmatrix(nrl, nrh, ncl, nch):
    return _pyOrbfit.dmatrix(nrl, nrh, ncl, nch)
dmatrix = _pyOrbfit.dmatrix

def dvector(nl, nh):
    return _pyOrbfit.dvector(nl, nh)
dvector = _pyOrbfit.dvector

def read_radec(*args):
    return _pyOrbfit.read_radec(*args)
read_radec = _pyOrbfit.read_radec

def fit_observations(obsarray, nobs, p, covar, logfile):
    return _pyOrbfit.fit_observations(obsarray, nobs, p, covar, logfile)
fit_observations = _pyOrbfit.fit_observations

def pbasis_to_helio(p, xv, partials):
    return _pyOrbfit.pbasis_to_helio(p, xv, partials)
pbasis_to_helio = _pyOrbfit.pbasis_to_helio

def add_to_obsarray(*args):
    return _pyOrbfit.add_to_obsarray(*args)
add_to_obsarray = _pyOrbfit.add_to_obsarray

def deghms(degr):
    return _pyOrbfit.deghms(degr)
deghms = _pyOrbfit.deghms

def orbitElements_helio(xv, orb):
    return _pyOrbfit.orbitElements_helio(xv, orb)
orbitElements_helio = _pyOrbfit.orbitElements_helio
# This file is compatible with both classic and new-style classes.

cvar = _pyOrbfit.cvar

