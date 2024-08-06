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
import os

DEGREE = pi/180

def main():
    startTime = time.time()
    ####################################################################################
    # initial wavefron definition
    ####################################################################################

    ini_wft_npix = [100, 100]
    ini_wft_range = [5e-3, 5e-3]
    beam_energy = 719.9
    pupil_position = 18.151
    sampling_factor = 0.2

    nMacroElec = 100

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

    meshPartCoh = deepcopy(wfr.mesh)

    ####################################################################################
    # electron beam
    ####################################################################################
    if srwl_uti_proc_is_master(): print('> Generating the electron beam ... ', end='')

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
    if srwl_uti_proc_is_master(): print('completed')

    ####################################################################################
    # undulator
    ####################################################################################

    if srwl_uti_proc_is_master(): print('> Generating the magnetic structure ... ', end='')

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
    if srwl_uti_proc_is_master(): print('completed')

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

    sampFactNxNyForProp = arPrecSR[6] 

    if srwl_uti_proc_is_master(): print('>  Undulator initial electric field calculation ... ', end='') 
    srwl.CalcElecFieldSR(wfr, eTraj, magFldCnt, arPrecSR)
    if srwl_uti_proc_is_master(): print('completed')

    if srwl_uti_proc_is_master(): wavefront_info(wfr)
    # srw_quick_plot(wfr, 'before')

    ####################################################################################
    # beamline definition
    ####################################################################################

    # pupil
    pupil = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=ini_wft_range[0], _Dy=ini_wft_range[0])
    pp_pupil = [0, 0, 1.0, 1, 0, 1., 5., 1., 5., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # M1A - bounce right mirror
    theta = 2.5*DEGREE
    m1a = SRWLOptMirPl(_size_tang=200e-3, _size_sag=18e-3, _npt=100, _nps=100,
                       _nvx=cos(theta),  _nvy=0, _nvz=-sin(theta), 
                       _tvx=-sin(theta), _tvy=0)

    pp_m1a = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # drift
    d2m1b = SRWLOptD(0.470422)
    pp_d2m1b = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # M1B - bounce left mirror
    m1b = SRWLOptMirTor(_rt=126.7433,_rs=1.802572,_size_tang=200e-3, _size_sag=10e-3, _npt=1001, _nps=1001,
                       _nvx=-cos(theta), _nvy=0, _nvz=-sin(theta), 
                       _tvx=-sin(theta), _tvy=0)

    pp_m1b = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # drift
    d2S1 = SRWLOptD(3.2095982206327456)
    pp_d2S1 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # mono entrance slit
    S1 = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=ini_wft_range[0], _Dy=ini_wft_range[0])
    pp_S1 = [0, 0, 1.0, 0, 0, 2/3, 3/2., 1, 4., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # drift - grating illumination
    d2G = SRWLOptD(0.6)
    pp_d2G = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # grating 450 l/mm - bounce up
    grazG = 2.302402*DEGREE
    outG = 0.460361*DEGREE
    deflection = (grazG+outG)*(1+0.840e-4)

    substrate = SRWLOptMirPl(_size_tang=80e-3, _size_sag=5e-3, _npt=5001, _nps=5001,
                             _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), 
                             _tvx=0, _tvy=-sin(grazG))

    vls = [450.00015185438964, -4.9135366794975125, -18.437531686217524, 4.9580003419891]

    G450 = SRWLOptG(_mirSub=substrate, _m=-1, _grDen=vls[0], _grDen1=vls[1]*0, _grDen2=vls[2]*0, _grDen3=vls[3]*0, _cff=0.2, _ang_graz=grazG)

    pp_G450 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, sin(2*deflection), cos(2*deflection), 1, 0]

    # drift
    d2m2 = SRWLOptD(0.311199)
    pp_d2m2 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    m2 = SRWLOptMirPl(_size_tang=120e-3, _size_sag=15e-3, _npt=500, _nps=500,
                      _nvx=0, _nvy=-cos(0.5*deflection), _nvz=-sin(0.5*deflection), 
                      _tvx=0, _tvy=-sin(0.5*deflection))

    pp_m2 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # drift
    d2m3 = SRWLOptD(0.239162)
    pp_d2m3 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # M3
    m3 = SRWLOptMirTor(_rt=83.0,_rs=0.1462124,_size_tang=120e-3, _size_sag=10e-3, _npt=1001, _nps=1001,
                            _nvx=-cos(1.2*DEGREE), _nvy=0, _nvz=-sin(1.2*DEGREE), 
                            _tvx= sin(1.2*DEGREE), _tvy=0)

    pp_m3 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # drift
    d2S2 = SRWLOptD(3.5)
    pp_d2S2 = [0, 0, 1.0, 1, 0, 1., 1., 1., 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # mono exit slit
    S2 = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=ini_wft_range[0], _Dy=ini_wft_range[0])
    pp_S2 = [0, 0, 1.0, 0, 0, 1., 1., 1, 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ####################################################################################
    # wavefront propagation
    ###################################################################################

    OE = [pupil,    m1a,    d2m1b,    m1b,    d2S1,    S1,    d2G,    G450,    d2m2,    m2,    d2m3,    m3,    d2S2,    S2]
    PP = [pp_pupil, pp_m1a, pp_d2m1b, pp_m1b, pp_d2S1, pp_S1, pp_d2G, pp_G450, pp_d2m2, pp_m2, pp_d2m3, pp_m3, pp_d2S2, pp_S2]

    branch = "mono"

    oe_branch = []
    pp_branch = []
        
    optBL = SRWLOptC(OE, PP)

    if srwl_uti_proc_is_master() is True and nMacroElec==1:
        print('> Simulating electric field propagation ... ', end='')
        pwfr = deepcopy(wfr)
        srwl.PropagElecField(pwfr, optBL)
        print('completed')
        wavefront_info(pwfr)
        srw_quick_plot(pwfr, me=1)

    if nMacroElec>1:
        calculation = 0
        strDataFolderName = "./results/"
        strIntPrtlChrnc = branch + '_' + str(nMacroElec/1000)+'k_ME.dat'
        print('- Simulating Partially-Coherent Wavefront Propagation... ') if(srwl_uti_proc_is_master()) else 0
        nMacroElecAvgPerProc = 10   # number of macro-electrons / wavefront to average on worker processes
        nMacroElecSavePer = 100     # intermediate data saving periodicity (in macro-electrons)
        srCalcMeth = 1              # SR calculation method
        srCalcPrec = 0.01           # SR calculation rel. accuracy
        radStokesProp = srwl_wfr_emit_prop_multi_e(e_beam, magFldCnt, meshPartCoh, srCalcMeth, srCalcPrec, nMacroElec,
                                                nMacroElecAvgPerProc, nMacroElecSavePer, os.path.join(os.getcwd(),
                                                strDataFolderName, strIntPrtlChrnc), sampFactNxNyForProp, optBL, _char=calculation)
        print('>> multi electron electron calculations: done') if(srwl_uti_proc_is_master()) else 0

    deltaT = time.time() - startTime
    hours, minutes = divmod(deltaT, 3600)
    minutes, seconds = divmod(minutes, 60)
    print("\n>>>> Elapsed time: " + str(int(hours)) + "h " + str(int(minutes)) + "min " + str(seconds) + "s ") if(srwl_uti_proc_is_master()) else 0




def srw_quick_plot(wfr, phase=False, me=0, backend="srw"):
    
    mesh0 = deepcopy(wfr.mesh)
    arI = array('f', [0]*mesh0.nx*mesh0.ny)
    srwl.CalcIntFromElecField(arI, wfr, 6, me, 3, mesh0.eStart, 0, 0)
    arIx = array('f', [0]*mesh0.nx)
    srwl.CalcIntFromElecField(arIx, wfr, 6, me, 1, mesh0.eStart, 0, 0)
    arIy = array('f', [0]*mesh0.ny)
    srwl.CalcIntFromElecField(arIy, wfr, 6, me, 2, mesh0.eStart, 0, 0)

    if phase:
        arP = array('d', [0]*mesh0.nx*mesh0.ny)
        srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, mesh0.eStart, 0, 0)
        arPx = array('d', [0]*mesh0.nx)
        srwl.CalcIntFromElecField(arPx, wfr, 0, 4, 1, mesh0.eStart, 0, 0)
        arPy = array('d', [0]*mesh0.ny)
        srwl.CalcIntFromElecField(arPy, wfr, 0, 4, 2, mesh0.eStart, 0, 0)

    plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
    plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
    uti_plot2d1d(arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity'])

    if phase:
        uti_plot2d1d(arP, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Phase'])

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