from .numcmp import *
import numpy as np
from mpmath import mp

def is_between(c, a, b, n = 2*np.pi):
    a = mp.fmod(a, n)
    b = mp.fmod(b, n)
    c = mp.fmod(c, n)
    
    if a < b:
        if a < c and c < b:
            return True;
        else:
            return False
    else: # b < a
        if b < c and c < a:
            return False  # if in [b, a] then not in [a, b]
        else:
            return True

def is_between_or_eq(c, a, b, n = 2*np.pi):
    a = mp.fmod(a, n)
    b = mp.fmod(b, n)
    c = mp.fmod(c, n)
    
    if a < b:
        if a <= c and c <= b:
            return True;
        else:
            return False
    else: # b < a
        if b <= c and c <= a:
            return False  # if in [b, a] then not in [a, b]
        else:
            return True

def diskp(z, r):
    """
    Returns whether a complex number z lies within the (open) disk or radius r
    """
    return lt(mp.fabs(z), r)

def closure_diskp(z, r):
    """
    Returns whether a complex number z lies within the closure of the open disk or radius r
    """
    return diskp(z, r) or eq(mp.fabs(z), r)

def wedgep(z, phi1, phi2):
    """
    Returns whether a complex number z lies within the wegde defined by (3.2)
    """
    return is_between(mp.arg(z), phi1, phi2)

def closure_wedgep(z, phi1, phi2):
    """
    Returns whether a complex number z lies within the wegde defined by (3.3)
    """
    return is_between_or_eq(mp.arg(z), phi1, phi2)


def in_G0(z, params):
    """
    Returns whether z lies in Region G0, eq. (3.4)
    """
    r0 = params['r0']
    return closure_diskp(z, r0)

def in_G1(z, params):
    """
    Returns whether z lies in Region G1, eq. (3.5)
    """
    r1 = params['r1']
    alpha = params['alpha']
    delta = params['delta']
    return (not diskp(z, r1)) and wedgep(z, -np.pi*alpha+delta, np.pi*alpha-delta) 

def in_G2(z, params):
    """
    Returns whether z lies in Region G2, eq. (3.6)
    """
    r1 = params['r1']
    alpha = params['alpha']
    deltat = params['deltat']
    return (not diskp(z, r1)) and wedgep(z, np.pi*alpha+deltat, -np.pi*alpha-deltat)

def in_G3(z, params):
    """
    Returns whether z lies in Region G3, eq. (3.7)
    """
    r1 = params['r1']
    alpha = params['alpha']
    delta = params['delta']
    deltat = params['deltat']
    #print((np.pi*alpha-delta)*180/pi, (np.pi*alpha+deltat)*180/pi)
    return (not diskp(z, r1)) and closure_wedgep(z, np.pi*alpha-delta, np.pi*alpha+deltat)

def in_G4(z, params):
    """
    Returns whether z lies in Region G4, eq. (3.8)
    """
    r1 = params['r1']
    alpha = params['alpha']
    delta = params['delta']
    deltat = params['deltat']
    return (not diskp(z, r1)) and closure_wedgep(z, -np.pi*alpha-deltat, -np.pi*alpha + delta)

def in_G5(z, params):
    """
    Returns whether z lies in Region G5, eq. (3.9)
    """
    r0 = params['r0']
    r1 = params['r1']
    alpha = params['alpha']
    return diskp(z, r1) and (closure_wedgep(z, -5*np.pi*alpha/6, 5*np.pi*alpha/6) and (not diskp(z, r0)))

def in_G6(z, params):
    """
    Returns whether z lies in Region G6, eq. (3.10)
    """
    r0 = params['r0']
    r1 = params['r1']
    alpha = params['alpha']
    return diskp(z, r1) and (wedgep(z, 5*np.pi*alpha/6, -5*np.pi*alpha/6) and (not diskp(z, r0)))
