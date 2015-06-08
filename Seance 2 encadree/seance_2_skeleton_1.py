#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
squelette de code pour la deuxième séance du projet de synthèse
'''
import math, numpy as np, scipy.integrate
import matplotlib.pyplot as plt
import pdb

import dynamic as dyn
import utils as ut
np.set_printoptions(precision=3, suppress=True, linewidth=200)


def get_trim(aircraft, h, Ma, sm, km):
    '''
    Calcul du trim pour un point de vol
    '''
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, h)
    Xe, Ue = dyn.trim(aircraft, {'va':va, 'h':h, 'gamma':0})
    return Xe, Ue

def get_all_trims(aircraft, hs, Mas, sms, kms):
    '''
    Calcul de trims pour une serie de points de vol
    '''
    trims = np.zeros((len(hs), len(Mas), len(sms), len(kms), 3))
    for i, h in enumerate(hs):
        for j, Ma in enumerate(Mas):
            for k, sm in enumerate(sms):
                for l,km in enumerate(kms):
                    Xe, Ue = get_trim(aircraft, h, Ma, sm, km)
                    trims[i, j, k, l] = (Xe[dyn. s_a], Ue[dyn.i_dm], Ue[dyn.i_dth])
    return trims

def plot_all_trims(aircraft, hs, Mas, sms, kms, trims, filename=None):
    '''
    Affichage graphique des valeurs de trim
    '''
    margins = (0.03, 0.05, 0.98, 0.95, 0.15, 0.46)
    fig = ut.prepare_fig(window_title='Trims {}'.format(aircraft.name), figsize=(20.48, 10.24), margins=margins)
    
    m=0
    for k, sm in enumerate(sms):
        for l,km in enumerate(kms):
            for i, h in enumerate(hs):
                for j, Ma in enumerate(Mas):
                    alpha, dphr, dth = trims[i, j, k, l]
                    fmt = 'alt {:5.0f} Ma {:.1f} sm {:.1f} km {:.1f} -> alpha {:5.2f} deg phr {:-5.1f} deg throttle {:.1f} %'
                    print fmt.format(h, Ma, sm, km, ut.deg_of_rad(alpha), ut.deg_of_rad(dphr), 100*dth)
            legend, params = ['Mach {}'.format(Ma) for Ma in Mas], r'\quad sm={} \quad km={}'.format(sm, km)
            ax = plt.subplot(4, 3, 3*m+1)
            plt.plot(hs, ut.deg_of_rad(trims[:, 0, k, l, 0]))
            plt.plot(hs, ut.deg_of_rad(trims[:, 1, k, l, 0]))
            ut.decorate(ax, r'$\alpha {}$'.format(params), r'altitude', '$deg$', legend=legend)
            ax = plt.subplot(4, 3, 3*m+2)
            plt.plot(hs, ut.deg_of_rad(trims[:, 0, k, l, 1]))
            plt.plot(hs, ut.deg_of_rad(trims[:, 1, k, l, 1]))
            ut.decorate(ax, r'$\delta_{{PHR}} {}$'.format(params), r'altitude', '$deg$', legend=legend)
            ax = plt.subplot(4, 3, 3*m+3)
            plt.plot(hs, trims[:, 0, k, l, 2]*100)
            plt.plot(hs, trims[:, 1, k, l, 2]*100)
            ut.decorate(ax, r'$\delta_{{th}} {}$'.format(params), r'altitude', '$\%$', legend=legend)
            m = m+1
    if filename<> None: plt.savefig(filename, dpi=160)


def plot_traj_trim(aircraft, h, Ma, sm, km):
    '''
    Affichage d'une trajectoire avec un point de trim comme condition initiale
    '''
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, h)
    Xe, Ue = dyn.trim(aircraft, {'va':va, 'h':h, 'gamma':0})
    time = np.arange(0., 100, 0.5)
    X = scipy.integrate.odeint(dyn.dyn, Xe, time, args=(Ue, aircraft))
    dyn.plot(time, X)


def get_linearized_model(aicraft, h, Ma, sm, km):
    '''
    Calcul numérique du modèle tangeant linéarisé pour un point de trim
    '''
    aircraft.set_mass_and_static_margin(km, sm)
    Xe, Ue = dyn.trim(aircraft, {'va':dyn.va_of_mach(Ma, h), 'h':h, 'gamma':0})
    A, B = ut.num_jacobian(Xe, Ue, aicraft, dyn.dyn)
    poles, vect_p = np.linalg.eig(A[dyn.s_va:, dyn.s_va:])
    return A, B, poles, vect_p


def plot_poles(aicraft, hs, Mas, sms, kms, filename=None):
    '''
    Affichage graphique de la position des poles pour une série de points de trim
    '''
    margins = (0.03, 0.05, 0.98, 0.95, 0.2, 0.38)
    fig = ut.prepare_fig(window_title='Poles {}'.format(aircraft.name), figsize=(20.48, 10.24), margins=margins)
    for i, h in enumerate(hs[::-1]):
        for j, Ma in enumerate(Mas):
            ax = plt.subplot(len(hs), len(Mas), i*len(Mas)+j+1)
            legend = []
            for k, sm in enumerate(sms):
                for l, km in enumerate(kms):
                    A, B, poles, vect_p = get_linearized_model(aicraft, h, Ma, sm, km)
                    print('{}'.format(poles))
                    plt.plot(poles.real, poles.imag, '*', markersize=10)
                    legend.append('ms={} km={}'.format(sm, km))
            ut.decorate(ax, r'$h={}m \quad  Ma={}$'.format(h, Ma), legend=legend,  xlim=[-3., 0.1])
    if filename<> None: plt.savefig(filename, dpi=160)                                    


aircraft = dyn.Param_A320()

hs, Mas = [3000, 11000], [0.5, 0.8]
sms, kms = [0.2, 1.], [0.1, 0.9]

trims = get_all_trims(aircraft, hs, Mas, sms, kms)

plot_all_trims(aircraft, hs, Mas, sms, kms, trims)

plot_traj_trim(aircraft, 5000, 0.5, 0.2, 0.5)

plot_poles(aircraft, hs, Mas, sms, kms)

plt.show()
