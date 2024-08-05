# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # XX: Simulating Wavefront Propagation through initial part of a 
# Soft X-Ray Undulator Radiation Beamline containing Variable Line Spacing (VLS) Grating
# Based on SRW example #12
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

__author__ = ['Rafael Celestre']
__contact__ = 'rafael.celestre@synchrotron-soleil.fr'
__license__ = 'GPL-3.0'
__copyright__ = 'Synchrotron SOLEIL, Saint Aubin, France'
__created__ = '01/AUG/2024'
__changed__ = '03/AUG/2024'

import sys

from copy import deepcopy
from oasys_srw.srwlib import *
from oasys_srw.uti_plot import *
from numpy import sin, cos, pi

DEGREE = pi/180

def main():
    ####################################################################################
    # initial wavefron definition
    ####################################################################################

    ini_wft_npix = [100, 100]
    ini_wft_range = [5e-3, 5e-3]
    beam_energy = 719.9
    pupil_position = 18.151
    sampling_factor = 0.2

    mesh = SRWLRadMesh(_eStart=beam_energy,
                            _eFin  =beam_energy,
                            _ne    =1,
                            _xStart= -ini_wft_range[0]/2,
                            _xFin  = ini_wft_range[0]/2,
                            _nx    = ini_wft_npix[0],
                            _yStart= -ini_wft_range[1]/2,
                            _yFin  = ini_wft_range[1]/2,
                            _ny    = ini_wft_npix[1],
                            _zStart=pupil_position)

    wfr = SRWLWfr()
    wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
    wfr.unitElFld = 1
    wfr.mesh = mesh

    ####################################################################################
    # electron beam
    ####################################################################################
    print('> Generating the electron beam ... ', end='')

    e_beam = SRWLPartBeam()
    e_beam.Iavg = 0.5

    e_beam.partStatMom1.x = 0.0
    e_beam.partStatMom1.y = 0.0
    e_beam.partStatMom1.z = 0.0 
    e_beam.partStatMom1.xp = 0.0
    e_beam.partStatMom1.yp = 0.0
    e_beam.partStatMom1.gamma = 5381.615754828152

    e_beam.arStatMom2[0] = 5.2873549684158404e-08       # <(x-<x>)^2>
    e_beam.arStatMom2[1] = 1.6376749051980198e-09       # <(x-<x>)(x'-<x'>)>
    e_beam.arStatMom2[2] = 9.363301879518787e-10        # <(x'-<x'>)^2>
    e_beam.arStatMom2[3] = 1.6661552574257426e-10       # <(y-<y>)^2>
    e_beam.arStatMom2[4] = 4.177394430693069e-11        # <(y-<y>)(y'-<y'>)>
    e_beam.arStatMom2[5] = 2.321201599884592e-11        # <(y'-<y'>)^2>
    e_beam.arStatMom2[10] = 0.001025**2  # <(E-<E>)^2>/<E>^2

    wfr.partBeam = e_beam
    print('completed')

    ####################################################################################
    # undulator
    ####################################################################################

    print('> Generating the magnetic structure ... ', end='')

    eTraj = 0

    e_beam.partStatMom1.z = -0.5*0.064*(28 + 4)

    und = SRWLMagFldU()
    und.set_sin(_per=0.064,
                _len=0.064*28, 
                _bx=0, 
                _by=0.1768940926246925, 
                _phx=0, _phy=0, _sx=-1, _sy=-1)
    magFldCnt = SRWLMagFldC(_arMagFld=[und],
                            _arXc=array('d', [0.0]),
                            _arYc=array('d', [0.0]),
                            _arZc=array('d', [0.0]))
    print('completed')

    ####################################################################################
    # generation of initial wavefront
    ####################################################################################

    arPrecSR = [0]*7
    arPrecSR[0] = 1     # SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    arPrecSR[1] = 0.01  # relative precision
    arPrecSR[2] = 0     # longitudinal position to start integration (effective if < zEndInteg)
    arPrecSR[3] = 0     # longitudinal position to finish integration (effective if > zStartInteg)
    arPrecSR[4] = 30000 # Number of points for trajectory calculation
    arPrecSR[5] = 1     # Use "terminating terms"  or not (1 or 0 respectively)
    arPrecSR[6] = sampling_factor # sampling factor for adjusting nx, ny (effective if > 0)

    print('>  Undulator initial electric field calculation ... ', end='')
    srwl.CalcElecFieldSR(wfr, eTraj, magFldCnt, arPrecSR)
    print('completed')

    wavefront_info(wfr)
    srw_quick_plot(wfr, 'before')

    ####################################################################################
    # beamline definition
    ####################################################################################

    pwfr = deepcopy(wfr)
    oe_array = []
    pp_array = []

    # pupil
    pupil = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=ini_wft_range[0], _Dy=ini_wft_range[0])
    oe_array.append(pupil)
    pp_array.append([0, 0, 1.0, 1, 0, 1., 5., 1., 5., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # M1A - bounce right mirror
    theta = 2.5*DEGREE
    m1a = SRWLOptMirPl(_size_tang=200e-3, _size_sag=18e-3, _npt=100, _nps=100,
                       _nvx=cos(theta),  _nvy=0, _nvz=-sin(theta), 
                       _tvx=-sin(theta), _tvy=0)

    oe_array.append(m1a)
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # drift
    oe_array.append(SRWLOptD(0.470422))
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # M1B - bounce left mirror

    m1b = SRWLOptMirTor(_rt=126.7433,_rs=1.802572,_size_tang=200e-3, _size_sag=10e-3, _npt=1001, _nps=1001,
                       _nvx=-cos(theta), _nvy=0, _nvz=-sin(theta), 
                       _tvx=-sin(theta), _tvy=0)

    oe_array.append(m1b)
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # drift
    oe_array.append(SRWLOptD(3.2095982206327456))
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # mono entrance slit
    slit = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=ini_wft_range[0], _Dy=ini_wft_range[0])
    oe_array.append(slit)
    pp_array.append([0, 0, 1.0, 0, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # drift - grating illumination
    oe_array.append(SRWLOptD(0.6))
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # grating 450 l/mm - bounce up
    grazG = 2.302402*DEGREE

    substrate = SRWLOptMirPl(_size_tang=80e-3, _size_sag=5e-3, _npt=1001, _nps=1001,
                             _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), 
                             _tvx=0, _tvy=-sin(grazG))

    vls = [450.00015185438964, -4.9135366794975125, -18.437531686217524, 4.9580003419891]

    g450 = SRWLOptG(_mirSub=substrate, _m=-1, _grDen=vls[0], _grDen1=vls[1], _grDen2=vls[2], _grDen3=vls[3], _cff=0.2, _ang_graz=grazG)
    deflection = g450.ang2cff(beam_energy,grazG)[1]*2

    oe_array.append(g450)
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, sin(deflection), cos(deflection), 1, 0])

    oe_array.append(SRWLOptD(0.311199))
    pp_array.append([0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    ####################################################################################
    # wavefront propagation
    ###################################################################################

    optBL = SRWLOptC(oe_array, pp_array)
    print('> Simulating electric field propagation ... ', end='')
    srwl.PropagElecField(pwfr, optBL)
    print('completed')
    wavefront_info(pwfr)
    srw_quick_plot(pwfr)


def srw_quick_plot(wfr, keyword):
    "Wrapper of SRW quick plotter"
    mesh0 = deepcopy(wfr.mesh)
    arI = array('f', [0]*mesh0.nx*mesh0.ny)
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
    arIx = array('f', [0]*mesh0.nx)
    srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
    arIy = array('f', [0]*mesh0.ny)
    srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)

    arP = array('d', [0]*mesh0.nx*mesh0.ny)
    srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, mesh0.eStart, 0, 0)
    arPx = array('d', [0]*mesh0.nx)
    srwl.CalcIntFromElecField(arPx, wfr, 0, 4, 1, mesh0.eStart, 0, 0)
    arPy = array('d', [0]*mesh0.ny)
    srwl.CalcIntFromElecField(arPy, wfr, 0, 4, 2, mesh0.eStart, 0, 0)

    plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
    plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
    uti_plot2d1d(arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', f'Intensity {keyword} Propagation'])
    uti_plot2d1d(arP, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', f'Phase {keyword} Propagation'])

    uti_plot_show()


def wavefront_info(wft):
    "Prints detailed information about the wavefront properties."

    Dx = wft.mesh.xFin - wft.mesh.xStart
    dx = Dx/(wft.mesh.nx-1)
    Dy = wft.mesh.yFin - wft.mesh.yStart
    dy = Dy/(wft.mesh.ny-1)

    print('\nWavefront information:')
    print(f'Nx = {wft.mesh.nx}, Ny = {wft.mesh.ny}')
    print(f'dx = {dx * 1E6:.4f} um, dy = {dy * 1E6:.4f} um')
    print(f'range x = {Dx * 1E3:.4f} mm, range y = {Dy * 1E3:.4f} mm')
    print(f'Rx = {wft.Rx:.6f} m, Ry = {wft.Ry:.6f} m')


if __name__=="__main__":
    main()