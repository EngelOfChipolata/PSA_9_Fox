#!/usr/bin/env python
#-*- coding: utf-8 -*-
import math, numpy as np
import matplotlib.pyplot as plt
import pdb

import dynamic as dyn
import utils as ut

def plot_thrust(P, filename=None):
    figure = ut.prepare_fig(None, u'Poussée {}'.format(P.name))
    hs = np.linspace(3000, 11000, 5)
    for h in hs:
        Machs = np.linspace(0.5, 0.8, 30)
        thrusts = np.zeros(len(Machs))
        for i in range(0, len(Machs)):
            thrusts[i] = dyn.propulsion_model([0, h, dyn.va_of_mach(Machs[i], h), 0, 0, 0], [0, 1., 0, 0], P)
        plt.plot(Machs, thrusts)
    ut.decorate(plt.gca(), u'Poussée maximum {}'.format(P.name), r'Mach', '$N$', 
                ['{} m'.format(h) for h in hs])
    if filename<> None: plt.savefig(filename, dpi=160)
    return figure

def CL(P, alpha, dphr):
    St_over_S = P.St/P.S
    CL0 = (St_over_S*0.25*P.CLat - P.CLa)*P.a0
    CLa = P.CLa + St_over_S*P.CLat*(1-0.25)
    CLdphr = St_over_S*P.CLat
    CL = CL0 + CLa*alpha + CLdphr * dphr
    return CL

def plot_CL(P, filename=None):
    alphas = np.linspace(ut.rad_of_deg(-10), ut.rad_of_deg(20), 30)
    dms = np.linspace(ut.rad_of_deg(20), ut.rad_of_deg(-30), 3)
    figure = ut.prepare_fig(None, u'Coefficient de Portance {}'.format(P.name))
    for dm in dms:
        plt.plot(ut.deg_of_rad(alphas), CL(P, alphas, dm))
    ut.decorate(plt.gca(), u'Coefficient de Portance', r'$\alpha$ en degres', '$C_L$')
    plt.legend(['$\delta _{{PHR}} =  ${:.1f}'.format(ut.deg_of_rad(dm)) for dm in dms], loc='best')
    if filename<> None: plt.savefig(filename, dpi=160)

def Cm(P, alpha, ms1):
    Cma = -ms1*P.CLa
    return P.Cm0 + Cma*(alpha-P.a0)

def plot_Cm(P, filename=None):
    alphas = np.linspace(ut.rad_of_deg(-10), ut.rad_of_deg(20), 30)
    mss = np.array([-0.1, 0, 0.2, 1])
    figure = ut.prepare_fig(None, u'Coefficient du moment de tangage {}'.format(P.name))
    for ms1 in mss:
        plt.plot(ut.deg_of_rad(alphas), Cm(P, alphas, ms1))
    ut.decorate(plt.gca(), u'Coefficient du moment de tangage', r'$\alpha$ en degres', '$C_m$')
    plt.legend(['ms =  {}'.format(ms) for ms in mss], loc='best')
    if filename<> None: plt.savefig(filename, dpi=160)

def dphr_e(P, alpha):
    return 1./P.Cmd*(P.ms*P.CLa*(alpha-P.a0) - P.Cm0)
    
def plot_dphr_e(P, filename=None):
    alphas = np.linspace(ut.rad_of_deg(-10), ut.rad_of_deg(20), 30)
    mss = [-0.1, 0., 0.2, 1.]
    figure = ut.prepare_fig(None, u'Équilibre {}'.format(P.name))
    for ms in mss:
        P.set_mass_and_static_margin(0.5, ms)
        dmes = np.array([dphr_e(P, alpha) for alpha in alphas])
        plt.plot(ut.deg_of_rad(alphas), ut.deg_of_rad(dmes))
    ut.decorate(plt.gca(), r'$\delta_{PHR_e}$', r'$\alpha$ en degres', r'$\delta_{PHR_e}$ en degres',
                ['$ms =  ${: .1f}'.format(ms) for ms in mss])


def plot_CLe(P):
    alphas = np.linspace(ut.rad_of_deg(-10), ut.rad_of_deg(20), 30)
    figure = ut.prepare_fig(None, u'Coefficient de portance équilibrée {}'.format(P.name))
    sms = [0.2, 1]
    for sm in sms:
        P.set_mass_and_static_margin(0.5, sm)
        dphres = [dphr_e(P, alpha) for alpha in alphas]
        CLes = [CL(P, alpha, dphr) for alpha, dphr in zip(alphas, dphres)]
        plt.plot(ut.deg_of_rad(alphas), CLes)
    ut.decorate(plt.gca(), u'Coefficient de portance équilibrée', r'$\alpha$ en degres', r'$CL_e$',
                ['$ms =  ${: .1f}'.format(sm) for sm in sms])

def plot_polar(P):
    #
    # TODO
    #
    pass

aircraft = dyn.Param_737_300() # use assigned aircraft
# plot_thrust(aircraft)
# plot_CL(aircraft)
plot_Cm(aircraft)
#plot_dphr_e(aircraft)
#plot_CLe(aircraft)
#plot_polar(aircraft)
plt.show()
