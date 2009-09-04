from __future__ import division
import TMM
reload(TMM)
from scipy import pi, sqrt, arange, vectorize, exp, c_, array, transpose, real, imag, rand, cos, sin, sinh, cosh, argmax, arange, eye, zeros, shape, poly1d
import pylab
import scipy
import pdb
import time
#import SAMII
import rwkascii
from rwkmisc import prettymat, colwise
from rwkos import FindFullPath
import shutil, os, sys, glob
import TMM.beam
reload(TMM.beam)
import TMM.rigid
reload(TMM.rigid)
import TMM.spring
reload(TMM.spring)
import TMM.velocitysource
reload(TMM.velocitysource)
import TMM.feedback
reload(TMM.feedback)
from TMM.beam import BeamElement
from TMM.rigid import RigidMass
from TMM.spring import TorsionalSpringDamper
from TMM.velocitysource import AngularVelocitySource
from TMM.velocitysource import AVSwThetaFB
from TMM.feedback import SAMIIAccelFB
from rwkdataproc import datastruct
import re
import rwkparse
import rwkbode
reload(rwkbode)
from rwkbode import bodeout

import copy

ms=4

#beam dimensions in inches
Li = 17.0
wi = 13.0/16
ti = 1.0/32#did not use calipers
#ti = 0.029

i_to_m = 25.4/1000.0
L = Li*i_to_m
w = wi*i_to_m
t = ti*i_to_m
A = w*t

#calculating rho from total mass
Lt = (19.0+5.0/8)*25.4/1000#total length (including the part in the clamp)
m = 65.75/1000.0#total mass (65.75 grams)
mu = m/Lt#mass per unit length

E = 200*10**9#N/m**2 or 200 GPa
rho = 7850.0#kg/m**3
I = 1.0/12*w*t**3
#mu = A*rho

EI = E*I

#from Frank
#EI = 0.4167
#mu = 0.2
#L = 0.508

class SLFRBeam(BeamElement):
    def __init__(self, beamparams={'EI':EI, \
                                   'mu':mu, \
                                   'L':L}, \
                 maxsize=ms, symname='Ubeam', \
                 symlabel='beam', symsub=True, usez=True):
        return BeamElement.__init__(self, beamparams, \
                                    maxsize=maxsize, \
                                    symlabel=symlabel, \
                                    symname=symname, \
                                    symsub=symsub, \
                                    usez=usez)


#accel dimiensions in inches
hi = 5.0/8
di = 7.0/16

#convert to meters
h = hi*i_to_m
d = di*i_to_m
Ra = d/2.0#radius of the cylinder

#mass in kilograms
ma = 7.7/1000#i.e. 8 grams

#I about centroid
Ic = 1.0/12*ma*(3*Ra**2+h**2)
#parallel axis theorem
Ia = Ic + ma*(h/2.0)**2

#note that the accel is attached at 16.5 inches from the base
class accel_mass(RigidMass):
    def __init__(self, accel_params={'m':ma, 'L':d, \
                                     'r':Ra,'I':Ia}, \
                 maxsize=ms, symname='Uaccel', symlabel='accel', \
                 symsub=True, usez=True):
        return RigidMass.__init__(self, accel_params,\
                                  maxsize=maxsize, symlabel=symlabel, \
                                  symname=symname, symsub=symsub, \
                                  usez=usez)


class AVS_only_model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, K, tau):
        self.avs=AngularVelocitySource({'K':K,'tau':tau}, \
                                       maxsize=ms, \
                                       symname='Uact', \
                                       symlabel='act', \
                                       unknownparams=['K','tau'])
        bodeout1 = {'input':'v', 'output':'th', 'type':'abs', \
                    'ind':self.avs, 'post':'', 'dof':1, \
                    'gain':180.0/pi*1024.0/360.0}
        self.bodeout1 = bodeout(**bodeout1)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                              [self.avs],
                              bodeouts=[self.bodeout1])
        

class SLFR_TMM_OL_model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, xf, actvect=[], include_spring=True):
        t = xf[-1]
        A = w*t
        I = 1.0/12*w*t**3
        mu = A*rho
        EI = E*I

        self.spring = TorsionalSpringDamper({'k':xf[0],'c':xf[1]}, \
                                       maxsize=ms,symname='Usp', \
                                       symlabel='sp', \
                                       unknownparams=['k','c'])
        self.beam = SLFRBeam(beamparams={'EI':EI, 'L':L, 'mu':mu})
        if actvect:#an empty actvect means the actuator is unknown
            self.avs = AngularVelocitySource({'K':actvect[0], \
                                              'tau':actvect[1]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act')
            ind=2
        else:
            self.avs = AngularVelocitySource({'K':xf[2],'tau':xf[3]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act', \
                                             unknownparams=['K','tau'])
            ind=4
        self.accel_tip = accel_mass()#a rigid mass at the tip of the beam
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.accel_tip, 'post':'accel', 'dof':0, \
                  'gain':xf[-2],'gainknown':False}
        if include_spring:
            b2_ind = self.spring
            my_list = [self.avs, self.spring, self.beam, self.accel_tip]
        else:
            b2_ind = self.avs
            my_list = [self.avs, self.beam, self.accel_tip]
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])


def calc_beam_props(t=1.0/32*i_to_m, w=13.0/16*i_to_m, L=17.0*i_to_m, \
                    rho=7850.0, E=200*10**9):
    A = w*t
    I = 1.0/12*w*t**3
    mu = A*rho
    EI = E*I
    return {'EI':EI, 'mu':mu, 'L':L}

    
class Beam_Only_Model(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SLFRBeam(beamparams=beamparams)
        my_list = [self.beam]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list)


class system_w_beam(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SLFRBeam(beamparams=beamparams)
    
class Beam_w_Accel_Mass(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        self.beam = SLFRBeam(beamparams=beamparams)
        self.accel = accel_mass()
        my_list = [self.beam, self.accel]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)
    
class Beam_Base_Spring_Accel_Mass(system_w_beam):
    def __init__(self, beamparams=None, k=1000.0):
        system_w_beam.__init__(self, beamparams)
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':0}, \
                                                  maxsize=ms)
        self.accel = accel_mass()
        my_list = [self.clamp_spring, self.beam, self.accel]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)
        
class Two_Piece_Beam_w_Accel_Mass(TMM.TMMSystem.ClampedFreeTMMSystem):
    def __init__(self, beamparams=None):
        if beamparams is None:
            beamparams = calc_beam_props()
        L1 = (16.0+11.0/32)*25.4/1000
        L2 = (5.0/16)*25.4/1000
        d1 = copy.copy(beamparams)
        d1['L'] = L1
        d2 = copy.copy(beamparams)
        d2['L'] = L2
        self.beam = SLFRBeam(beamparams=d1)
        self.accel = accel_mass()
        self.tip_beam = SLFRBeam(beamparams=d2)
        my_list = [self.beam, self.accel, self.tip_beam]
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self, \
                                                           my_list)


class SLFR_TMM_OL_model_v2(system_w_beam):
    def __init__(self, xf, actvect=[], include_spring=True, \
                 beamparams=None):
        k = xf[-1]#default should be 45.0
        if beamparams is None:
            beamparams = calc_beam_props()
        beamparams['mu'] = beamparams['mu']*xf[-3]
        system_w_beam.__init__(self, beamparams)
        self.clamp_spring = TorsionalSpringDamper({'k':k,'c':0}, \
                                                  maxsize=ms)
        self.spring = TorsionalSpringDamper({'k':xf[0],'c':xf[1]}, \
                                       maxsize=ms,symname='Usp', \
                                       symlabel='sp', \
                                       unknownparams=['k','c'])
        if actvect:#an empty actvect means the actuator is unknown
            self.avs = AngularVelocitySource({'K':actvect[0], \
                                              'tau':actvect[1]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act')
            ind=2
        else:
            self.avs = AngularVelocitySource({'K':xf[2],'tau':xf[3]}, \
                                             maxsize=ms, \
                                             symname='Uact', \
                                             symlabel='act', \
                                             unknownparams=['K','tau'])
            ind=4
        self.accel_tip = accel_mass()#a rigid mass at the tip of the beam
        bodeout1={'input':'v', 'output':'atip', 'type':'abs', \
                  'ind':self.accel_tip, 'post':'accel', 'dof':0, \
                  'gain':xf[-2],'gainknown':False}
        if include_spring:
            b2_ind = self.spring
            my_list = [self.avs, self.spring, self.clamp_spring, \
                       self.beam, self.accel_tip]
        else:
            b2_ind = self.avs
            my_list = [self.avs, self.beam, self.accel_tip]
        bodeout2={'input':'v', 'output':'th', 'type':'abs', \
                  'ind':b2_ind, 'post':'', 'dof':1, \
                  'gain':180.0/pi*1024.0/360.0}
        self.bodeout1=bodeout(**bodeout1)
        self.bodeout2=bodeout(**bodeout2)
        return TMM.TMMSystem.ClampedFreeTMMSystem.__init__(self,
                                                           my_list, \
                              bodeouts=[self.bodeout1, self.bodeout2])
