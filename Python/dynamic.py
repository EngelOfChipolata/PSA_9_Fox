#-*- coding: utf-8 -*-
'''
Dynamic model for a 3 Degrees Of Freedom longitudinal aircraft  
'''
import math, numpy as np, scipy.optimize
import matplotlib.pyplot as plt

import utils as ut

''' naming of state and input components '''
s_y, s_h, s_va, s_a, s_th, s_q, s_size = range(0, 7)
i_dm, i_dth, i_wy, i_wz, i_size =  range(0, 5)

def get_mach(va, T, k=1.4, Rs=287.05): return va/math.sqrt(k*Rs*T)

def va_of_mach(m, h, k=1.4, Rs=287.05):
    p, rho, T = ut.isa(h)
    return m*math.sqrt(k*Rs*T)

def propulsion_model(X, U, P):
    p, rho, T = ut.isa(X[s_h])
    rho0 = 1.225
    mach = get_mach(X[s_va], T)
    return P.F0*math.pow(rho/rho0, 0.6)*(0.568+0.25*math.pow(1.2-mach, 3))*U[i_dth]
    
def get_aero_forces_and_moments(X, U, P):
    p, rho, T = ut.isa(X[s_h])
    pdyn = 0.5*rho*X[s_va]**2
    da = X[s_a]-P.a0
    CLw = P.CLa*da                                            
    at = X[s_a] - 0.25*da + U[i_dm] + P.CLq*P.lt/X[s_va]*X[s_q]
    CLt = P.St/P.S*P.CLat*at                                   
    CL = CLw + CLt
    L = pdyn*P.S*CL
    D = pdyn*P.S*(P.CD0 + P.ki*CL**2)
    Cm = P.Cm0 - P.ms*CLw + P.Cmq*P.lt/X[s_va]*X[s_q] + P.Cmd*U[i_dm]
    M = pdyn*P.S*P.cbar*Cm 
    return L, D, M

def dyn(X, t, U, P):
    '''  Dynamic model '''
    Xdot = np.zeros(s_size)
    gamma_a = X[s_th] - X[s_a]  # air path angle
    cg, sg = math.cos(gamma_a), math.sin(gamma_a)
    ca, sa = math.cos(X[s_a]), math.sin(X[s_a])
    L, D, M = get_aero_forces_and_moments(X, U, P)
    F = propulsion_model(X, U, P)
    Xdot[s_y] = X[s_va] * cg - U[i_wy]
    Xdot[s_h] = X[s_va] * sg - U[i_wz]
    Xdot[s_va] = (F*ca-D)/P.m-P.g*sg 
    Xdot[s_a] = X[s_q] - (L+F*sa)/P.m/X[s_va]  + P.g/X[s_va]*cg 
    Xdot[s_th] = X[s_q]
    Xdot[s_q] = M/P.Iyy
    return Xdot

def trim(P, args=None):
    va = args.get('va', 100.)
    gamma = args.get('gamma', 0.)
    h = args.get('h', 5000.)
    wy, wz = 0, 0
    def err_func((dm, dth, alpha)):
        theta = gamma + alpha
        U = np.array([dm, dth, wy, wz])
        X = np.array([0., h, va, alpha, theta, 0])
        L, D, M = get_aero_forces_and_moments(X, U, P)
        F = propulsion_model(X, U, P)
        cg, sg = math.cos(gamma), math.sin(gamma)
        ca, sa = math.cos(alpha), math.sin(alpha)
        return [(F*ca-D)/P.m-P.g*sg, -(L+F*sa)/P.m + P.g*cg, M]
    p0 = [ut.rad_of_deg(0.), 0.5, ut.rad_of_deg(1.)]
    sol = scipy.optimize.root(err_func, p0, method='hybr')
    dm, dth, alpha = sol.x
    X, U = [0, h, va, alpha, gamma+alpha, 0], [dm, dth, 0, 0]
    return X, U


class Param:
    def __init__(self):
        self.g = 9.81
        self.m_k = 0.5                
        self.ms = 0.3                 # static margin
        # aero
        self.a0 = ut.rad_of_deg(-2.)  # zero lift angle
        self.CD0 = 0.025              # zero drag coefficient
        self.Cm0 = -0.59              # zero moment coefficient
        self.CLq = 1.3

    def set_mass_and_static_margin(self, km, sm):
        self.m_k = km
        self.ms = sm
        self.compute_auxiliary()

    def compute_auxiliary(self):
        self.m = (1-self.m_k)*self.m_OWE + self.m_k*self.m_MTOW
        self.Iyy = 0.5*(1./12.*self.m*self.l_fus**2)

        self.lt = 0.5 * self.l_fus           # CG to tail distance

        self.CLa = math.pi*self._lambda/(1.+math.sqrt(1+(0.5*self._lambda)**2))
        self.ki = 1./(math.pi*self._lambda)

        self.CLat = math.pi*self._lambdat/(1.+math.sqrt(1+(0.5*self._lambdat)**2))
        
        Vt = self.lt*self.St/self.cbar/self.S # tail volume
        self.Cmd = -Vt*self.CLat
        self.Cmq = self.Cmd * self.CLq 


class Param_A320(Param):
    def __init__(self, m_k=0.5):
        Param.__init__(self)
        self.name = 'Airbus A-320'
        self.m_OWE  = 39733.
        self.m_MTOW = 73500.
        self.m_k = m_k

        self.l_fus    =  37.57          # length of fuselage
        self.cbar     =   4.19          # wing reference chord
        self.St       =  31.            # tail lifting surface
        self.S        = 122.44          # wing surface
        self._lambdat =   5.            # tail aspect ratio
        self._lambda  =   9.39          # wing aspect ratio
        
        self.F0       = 2.*111205       # engine max thrust
        self.eng_name = 'CFM 56-5A1' 

        self.compute_auxiliary()


class Param_737_800(Param):
    def __init__(self, ):
        Param.__init__(self)
        self.name = 'Boeing 737-800'
        self.m_OWE  = 41413.
        self.m_MTOW = 70534.
 
        self.l_fus    =  38.02          # length of fuselage
        self.cbar     =   4.17          # wing reference chord
        self.St       =  32.8           # tail lifting surface
        self.S        = 124.6           # wing surface
        self._lambdat =   6.28          # tail aspect ratio
        self._lambda  =   9.45          # wing aspect ratio
        
        self.F0       = 2.*106757       # engine max thrust
        self.eng_name = 'CFM 56-7B24' 

        self.compute_auxiliary()


class Param_A319(Param):
    def __init__(self):
        Param.__init__(self)
        self.name = 'Airbus A-319'
        self.m_OWE  = 39358.
        self.m_MTOW = 64000.

        self.l_fus    =  33.84          # length of fuselage
        self.cbar     =   4.19          # wing reference chord
        self.St       =  31.0           # tail lifting surface
        self.S        = 122.44          # wing surface
        self._lambdat =   5.            # tail aspect ratio
        self._lambda  =   9.39          # wing aspect ratio

        self.F0       = 2.*97860        # engine max thrust
        self.eng_name = 'CFM 56-5B5' 

        self.compute_auxiliary()


class Param_A321(Param):
    def __init__(self):
        Param.__init__(self)
        self.name = 'Airbus A-321'
        self.m_OWE  = 47000.
        self.m_MTOW = 89000.

        self.l_fus    =  44.51          # length of fuselage
        self.cbar     =   4.34          # wing reference chord
        self.St       =  31.0           # tail lifting surface
        self.S        = 126.00          # wing surface
        self._lambdat =   5.            # tail aspect ratio
        self._lambda  =   9.13          # wing aspect ratio

        self.F0       = 2.*133446       # engine max thrust
        self.eng_name = 'CFM 56-5B1' 

        self.compute_auxiliary()


class Param_737_700(Param):
    def __init__(self):
        Param.__init__(self)
        self.name = 'Boeing 737-700'
        self.m_OWE  = 37648.
        self.m_MTOW = 60326.

        self.l_fus    =  32.18          # length of fuselage
        self.cbar     =   4.17          # wing reference chord
        self.St       =  32.8           # tail lifting surface
        self.S        = 124.60          # wing surface
        self._lambdat =   6.28          # tail aspect ratio
        self._lambda  =   9.44          # wing aspect ratio

        self.F0       = 2.*61633        # engine max thrust
        self.eng_name = 'CFM 56-7B20' 

        self.compute_auxiliary()


class Param_737_300(Param):
    def __init__(self):
        Param.__init__(self)
        self.name = 'Boeing 737-300'
        self.m_OWE  = 31480.
        self.m_MTOW = 56473.

        self.l_fus    =  32.18          # length of fuselage
        self.cbar     =   3.73          # wing reference chord
        self.St       =  31.31          # tail lifting surface
        self.S        =  91.04          # wing surface
        self._lambdat =   5.15          # tail aspect ratio
        self._lambda  =   9.16          # wing aspect ratio

        self.F0       = 2.*88694        # engine max thrust
        self.eng_name = 'CFM 56-3B1' 

        self.compute_auxiliary()

all_ac_types = [Param_A320, Param_737_800, Param_A319, Param_A321, Param_737_700, Param_737_300]

def plot(time, X, U=None, figure=None, window_title="Trajectory"):
    figure = ut.prepare_fig(figure, window_title, (20.48, 10.24))
    plots = [("$y$",       "m",     X[:,s_y], None),
             ("$h$",       "m",     X[:,s_h], 2.),
             ("$v_a$",     "m/s",   X[:,s_va], 1.),
             ("$\\alpha$", "deg",   ut.deg_of_rad(X[:,s_a]), 2.),
             ("$\\theta$", "deg",   ut.deg_of_rad(X[:,s_th]), 2.),
             ("$q$",       "deg/s", ut.deg_of_rad(X[:,s_q]), 2.)]
    for i, (title, ylab, data, min_yspan) in enumerate(plots):
        ax = plt.subplot(3, 2, i+1)
        plt.plot(time, data)
        ut.decorate(ax, title=title, ylab=ylab, min_yspan=min_yspan)
    return figure
