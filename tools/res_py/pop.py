#
# implementation of the popovici algo
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019, apr-2022
# @license GPLv2
#
# @desc for algorithm: [pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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
import numpy.linalg as la
import reso
import helpers


IDX_SRC_Y      = 0
IDX_SRC_Z      = 1
IDX_MONO_X     = 2
IDX_MONO_Y     = 3
IDX_MONO_Z     = 4
IDX_SAMPLE_X   = 5
IDX_SAMPLE_Y   = 6
IDX_SAMPLE_Z   = 7
NUM_MONO_POS   = 8

IDX_SRC_H      = 0
IDX_SRC_V      = 1
IDX_MONO_H     = 2
IDX_MONO_V     = 3
IDX_SAMPLE_H   = 4
IDX_SAMPLE_V   = 5
IDX_ANA_H      = 6
IDX_ANA_V      = 7
NUM_MONO_COMPS = 4
NUM_COMPS      = 8

IDX_KI_H       = 0
IDX_KI_V       = 1
IDX_KF_H       = 2
IDX_KF_V       = 3
NUM_KI         = 2
NUM_KIKF       = 4



#
# get matrices for ki axis
#
def get_mono_trafos(dist_src_mono, dist_mono_sample,
	thetam, thetas, inv_curvh, inv_curvv):

    s_th_m = np.sin(thetam)
    c_th_m = np.cos(thetam)
    s_th_s = np.sin(thetas)
    c_th_s = np.cos(thetas)


    # D matrix, [pop75], Appendix 2
    D = np.zeros([NUM_MONO_COMPS, NUM_MONO_POS])

    D[IDX_SRC_H, IDX_SRC_Y] = -1. / dist_src_mono
    D[IDX_SRC_H, IDX_MONO_X] = -c_th_m / dist_src_mono
    D[IDX_SRC_H, IDX_MONO_Y] = s_th_m / dist_src_mono

    D[IDX_SRC_V, IDX_SRC_Z] = -1. / dist_src_mono
    D[IDX_SRC_V, IDX_MONO_Z] = 1. / dist_src_mono

    D[IDX_MONO_H, IDX_MONO_X] = c_th_m / dist_mono_sample
    D[IDX_MONO_H, IDX_MONO_Y] = s_th_m / dist_mono_sample
    D[IDX_MONO_H, IDX_SAMPLE_X] = s_th_s / dist_mono_sample
    D[IDX_MONO_H, IDX_SAMPLE_Y] = c_th_s / dist_mono_sample

    D[IDX_MONO_V, IDX_MONO_Z] = -1. / dist_mono_sample
    D[IDX_MONO_V, IDX_SAMPLE_Z] = 1. / dist_mono_sample


    # T matrix, [pop75], Appendix 2
    T = np.zeros([NUM_KI, NUM_MONO_POS])

    T[IDX_KI_H, IDX_SRC_Y] = 0.5 * D[IDX_SRC_H, IDX_SRC_Y]
    T[IDX_KI_H, IDX_MONO_X] = 0.5 * (D[IDX_SRC_H, IDX_MONO_X] + D[IDX_SRC_H, IDX_MONO_X])
    T[IDX_KI_H, IDX_MONO_Y] = 0.5 * (D[IDX_SRC_H, IDX_MONO_Y] + D[IDX_SRC_H, IDX_MONO_Y]) - inv_curvh
    T[IDX_KI_H, IDX_SAMPLE_X] = 0.5 * D[IDX_MONO_H, IDX_SAMPLE_X]
    T[IDX_KI_H, IDX_SAMPLE_Y] = 0.5 * D[IDX_MONO_H, IDX_SAMPLE_Y]

    T[IDX_KI_V, IDX_SRC_Z] = 0.5 * D[IDX_SRC_V, IDX_SRC_Z] / s_th_m
    T[IDX_KI_V, IDX_MONO_Z] = 0.5 * (D[IDX_SRC_V, IDX_MONO_Z] - D[IDX_SRC_V, IDX_MONO_Z]) / s_th_m - inv_curvv
    T[IDX_KI_V, IDX_SAMPLE_Z] = -0.5 * D[IDX_MONO_V, IDX_SAMPLE_Z] / s_th_m


    return [D, T]



#
# unite trafo matrices
#
def combine_mono_ana_trafos(Dm, Tm, Da, Ta):
    N = Dm.shape[0]
    M = Dm.shape[1]
    D = np.zeros([2*N, 2*M - 3])

    D[0:N, 0:M] = Dm
    # add ana matrix (without sample columns) in lower right block
    D[N:2*N, M:2*M-3] = Da[:, 0:IDX_SAMPLE_X]
    # unite mono and ana sample columns
    D[N:2*N, IDX_SAMPLE_X:IDX_SAMPLE_Z+1] = Da[:, IDX_SAMPLE_X:IDX_SAMPLE_Z+1]


    N = Tm.shape[0]
    M = Tm.shape[1]
    T = np.zeros([2*N, 2*M - 3])

    T[0:N, 0:M] = Tm
    # add ana matrix (without sample columns) in lower right block
    T[N:2*N, M:2*M-3] = Ta[:, 0:IDX_SAMPLE_X]
    # unite mono and ana sample columns
    T[N:2*N, IDX_SAMPLE_X:IDX_SAMPLE_Z+1] = Ta[:, IDX_SAMPLE_X:IDX_SAMPLE_Z+1]

    return [D, T]



#
# trafo from the [dki_x, dki_y, dki_z, dkf_x, dkf_y, dkf_z]
#       into the [dQ_x, dQ_y, dQ_z, dE] system
# see [mit84], equ. A.15
#
def get_QE_trafo(ki, kf, ki_Q, kf_Q):
    B = np.zeros([6, 6])

    B[0:2, 0:2] = +helpers.rotation_matrix_2d(ki_Q)
    B[0:2, 3:5] = -helpers.rotation_matrix_2d(kf_Q)

    # dQ = dki - dkf
    B[2, 2] = +1.
    B[2, 5] = -1.

    #  E ~ ki^2 - kf^2
    # dE ~ 2ki*dki - 2kf*dkf
    B[3, 0] = +2.*ki * helpers.ksq2E
    B[3, 3] = -2.*kf * helpers.ksq2E

    return B



#
# resolution algorithm
#
def calc(param):
    twotheta = param["twotheta"] * param["sample_sense"]
    thetas = twotheta * 0.5
    thetam = param["thetam"] * param["mono_sense"]
    thetaa = param["thetaa"] * param["ana_sense"]
    ki_Q = param["angle_ki_Q"] * param["sample_sense"]
    kf_Q = param["angle_kf_Q"] * param["sample_sense"]

    ki = param["ki"]
    kf = param["kf"]
    E = param["E"]
    Q = param["Q"]

    lam = helpers.k2lam(ki)


    # --------------------------------------------------------------------
    # mono/ana focus
    mono_curvh = param["mono_curvh"]
    mono_curvv = param["mono_curvv"]
    ana_curvh = param["ana_curvh"]
    ana_curvv = param["ana_curvv"]

    if param["mono_is_optimally_curved_h"]:
        mono_curvh = foc_curv(param["dist_src_mono"], \
            param["dist_mono_sample"], np.abs(2.*thetam), False)
    if param["mono_is_optimally_curved_v"]:
        mono_curvv = foc_curv(param["dist_src_mono"], \
            param["dist_mono_sample"], np.abs(2.*thetam), True)
    if param["ana_is_optimally_curved_h"]:
        ana_curvh = foc_curv(param["dist_sample_ana"], \
            param["dist_ana_det"], np.abs(2.*thetaa), False)
    if param["ana_is_optimally_curved_v"]:
        ana_curvv = foc_curv(param["dist_sample_ana"], \
            param["dist_ana_det"], np.abs(2.*thetaa), True)

    inv_mono_curvh = 0.
    inv_mono_curvv = 0.
    inv_ana_curvh = 0.
    inv_ana_curvv = 0.

    if param["mono_is_curved_h"]:
        inv_mono_curvh = 1./mono_curvh
    if param["mono_is_curved_v"]:
        inv_mono_curvv = 1./mono_curvv
    if param["ana_is_curved_h"]:
        inv_ana_curvh = 1./ana_curvh
    if param["ana_is_curved_v"]:
        inv_ana_curvv = 1./ana_curvv
    # --------------------------------------------------------------------


    # dict with results
    res = {}

    res["Q_avg"] = np.array([ Q, 0., 0., E ])


    # -------------------------------------------------------------------------
    # - if the instruments works in kf=const mode and the scans are counted for
    #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
    # - if the instrument works in ki=const mode the kf^3 factor is needed.

    # TODO
    #tupScFact = get_scatter_factors(param.flags, param.thetam, param.ki, param.thetaa, param.kf);
    tupScFact = [1., 1., 1.]

    dmono_refl = param["dmono_refl"] * tupScFact[0]
    dana_effic = param["dana_effic"] * tupScFact[1]
    dxsec = tupScFact[2]
    #if param.mono_refl_curve:
    #    dmono_refl *= (*param.mono_refl_curve)(param.ki)
    #if param.ana_effic_curve:
    #    dana_effic *= (*param.ana_effic_curve)(param.kf)
    # -------------------------------------------------------------------------


    # -------------------------------------------------------------------------
    # collimators and crystal mosaics
    coll_h_pre_mono = param["coll_h_pre_mono"]
    coll_v_pre_mono = param["coll_v_pre_mono"]

    if param["use_guide"]:
        coll_h_pre_mono = lam*param["guide_div_h"]
        coll_v_pre_mono = lam*param["guide_div_v"]

    # G matrix, [pop75], Appendix 1
    G = np.zeros([NUM_COMPS, NUM_COMPS])
    G[IDX_SRC_H, IDX_SRC_H] = 1. / coll_h_pre_mono**2.
    G[IDX_SRC_V, IDX_SRC_V] = 1. / coll_v_pre_mono**2.
    G[IDX_MONO_H, IDX_MONO_H] = 1. / param["coll_h_pre_sample"]**2.
    G[IDX_MONO_V, IDX_MONO_V] = 1. / param["coll_v_pre_sample"]**2.
    G[IDX_SAMPLE_H, IDX_SAMPLE_H] = 1. / param["coll_h_post_sample"]**2.
    G[IDX_SAMPLE_V, IDX_SAMPLE_V] = 1. / param["coll_v_post_sample"]**2.
    G[IDX_ANA_H, IDX_ANA_H] = 1. / param["coll_h_post_ana"]**2.
    G[IDX_ANA_V, IDX_ANA_V] = 1. / param["coll_v_post_ana"]**2.

    # F matrix, [pop75], Appendix 1
    F = np.zeros([NUM_KIKF, NUM_KIKF])
    F[IDX_KI_H, IDX_KI_H] = 1. / param["mono_mosaic"]**2.
    F[IDX_KI_V, IDX_KI_V] = 1. / param["mono_mosaic_v"]**2.
    F[IDX_KF_H, IDX_KF_H] = 1. / param["ana_mosaic"]**2.
    F[IDX_KF_V, IDX_KF_V] = 1. / param["ana_mosaic_v"]**2.
    # -------------------------------------------------------------------------


    # -------------------------------------------------------------------------
    # ki and kf trafo matrices
    [Dm, Tm] = get_mono_trafos(param["dist_src_mono"], param["dist_mono_sample"], \
        thetam, thetas, inv_mono_curvh, inv_mono_curvv)
    [Da, Ta] = get_mono_trafos(param["dist_ana_det"], param["dist_sample_ana"], \
        thetaa, thetas, inv_ana_curvh, inv_ana_curvv)

    [D, T] = combine_mono_ana_trafos(Dm, Tm, Da, Ta)

    B = get_QE_trafo(ki, kf, ki_Q, kf_Q)
    # -------------------------------------------------------------------------


    #
    # TODO
    #


    res["ok"] = False
    return res
