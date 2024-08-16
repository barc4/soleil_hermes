# PyOptiX benchmarking: preparing our beamlines for the SOLEIL upgrade

![alt text](https://github.com/barc4/soleil_hermes/blob/main/poster.png)

## PyOptiX

You can get PyOptiX [here](https://github.com/ddennetiere/pyoptix) and its C core (here)[https://github.com/ddennetiere/pyoptix].


# HERMES beamline
repository for [Hermes](https://doi.org/10.1107/S1600577515007778) simulations using [PyOptiX](https://github.com/ddennetiere/pyoptix), [SHADOW3](https://github.com/oasys-kit/shadow3) and [4](https://github.com/oasys-kit/shadow4/) as well as [SRW](https://github.com/ochubar/SRW)


## The beamline

HERMES is a phase III beamline dedicated to X-ray microscopy. It combines two types of microscopy. The first one is a Photon-Photon microscopy: Scanning Transmission X-Ray Microscopy (STXM). The second type is Photo-Electron microscopy: X-Ray Photoemitted ElectronMicroscopy (X-PEEM). This original approach offers to users two very complementary techniques, each one being adapted to specific sample environments.  

The STXM essentially enables to probe the sample volume’s properties, with depths of analysis around few hundred nanometers. Besides, X-PEEM  is a surface sensitive technique adapted to an ultra-vacuum environment. It mainly enables to scan the surface’s first few nanometers.
Both techniques carry out X-Ray spectroscopy as a contrast method, in this case a chemical contrast.     

Alongside, other contrast means can be added by making use of the spectroscopic techniques’ specificities. Thus, using circular and linear polarization of light gives an access to magnetic domains, ferromagnetic and antiferromagnetic (XMCD and XMLD)…
Finally, both methods can be used not only to get images of samples, but also to realize some measurement in local spectroscopy (XAS, XANES, XPS, ARPES…) at a nanoscopic scale.

## Technical data

HERMES beamline is installed on a medium straight section (IM9). The beamline covers the energy range going from 70 eV to 2,5 keV. The beamline starts right below the Si2p-edge, covers the K-edge of light elements (C, N, O...) , L-edge of transition metals, M-edge of rare earth and eventually the K-edge of Si, S and P. Furthermore, the beamline enables to work with the variable linear polarization (horizontal to vertical), in  addition to circular and elliptical polarizations.
To enable this, the beamline has two undulators, Apple-II type : HU-64 (1.7m, 25 periods) covering the 70-600 eV range; HU-42 (1.8m, 42 periods) for the high energy range, going from 0,5 to 2,5 keV.
