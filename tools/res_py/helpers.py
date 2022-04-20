#
# helper functions
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license GPLv2
#
# @desc for reso algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984)
# @desc for vertical scattering modification: [eck20] G. Eckold, personal communication, 2020.
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

# requires numpy version >= 1.10
import numpy as np


#--------------------------------------------------------------------------
# scattering triangle
# see: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/tascalc.py

ksq2E = 2.072124836832


def k2lam(k):
    return 2.*np.pi / k


def get_E(ki, kf):
        return ksq2E * (ki**2. - kf**2.)


def get_scattering_angle(ki, kf, Q):
    c = (ki**2. + kf**2. - Q**2.) / (2.*ki*kf)
    return np.arccos(c)


def get_angle_ki_Q(ki, kf, Q):
    c = (ki**2. + Q**2. - kf**2.) / (2.*ki*Q)
    return np.arccos(c)


def get_angle_kf_Q(ki, kf, Q):
    c = (ki**2. - Q**2. - kf**2.) / (2.*kf*Q)
    return np.arccos(c)


def get_mono_angle(k, d):
    s = np.pi/(d*k)
    angle = np.arcsin(s)
    return angle
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# helpers

#
# z rotation matrix
#
def rotation_matrix_3d_z(angle):
    s = np.sin(angle)
    c = np.cos(angle)

    return np.array([
        [c, -s, 0],
        [s,  c, 0],
        [0,  0, 1]])


def mirror_matrix(iSize, iComp):
    mat = np.identity(iSize)
    mat[iComp, iComp] = -1.

    return mat;


#
# thin lens equation: 1/f = 1/lenB + 1/lenA
#
def focal_len(lenBefore, lenAfter):
    f_inv = 1./lenBefore + 1./lenAfter
    return 1. / f_inv


#
# optimal mono/ana curvature,
# see e.g.
#    - (Shirane 2002) p. 66
#    - or nicos/nicos-core.git/tree/nicos/devices/tas/mono.py in nicos
#    - or Monochromator_curved.comp in McStas
#
def foc_curv(lenBefore, lenAfter, tt, bVert):
    f = focal_len(lenBefore, lenAfter)
    s = np.abs(np.sin(0.5*tt))

    if bVert:
        curv = 2. * f*s
    else:
        curv = 2. * f/s

    return curv
#--------------------------------------------------------------------------
