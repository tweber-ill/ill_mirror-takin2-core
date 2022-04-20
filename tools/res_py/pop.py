#
# implementation of the popovici algo
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license GPLv2
#
# @desc for algorithm: [pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
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
import numpy.linalg as la
import reso
import helpers


#
# get matrices for ki axis 
#
def get_mono_vals(src_w, src_h, mono_w, mono_h,
	dist_src_mono, dist_mono_sample,
	ki, thetam, thetas,
	coll_h_pre_mono, coll_h_pre_sample,
	coll_v_pre_mono, coll_v_pre_sample,
	mono_mosaic, mono_mosaic_v,
	inv_mono_curvh, inv_mono_curvv,
	refl):

    s_th_m = np.sin(thetam)
    c_th_m = np.cos(thetam)
    s_th_s = np.sin(thetas)
    c_th_s = np.cos(thetas)


    IDX_SRC_Y    = 0
    IDX_MONO_X   = 1
    IDX_MONO_Y   = 2
    IDX_SAMPLE_X = 3
    IDX_SAMPLE_Y = 4
    IDX_SRC_Z    = 5
    IDX_MONO_Z   = 6
    IDX_SAMPLE_Z = 7


    # D matrix, [pop75], Appendix 2
    D = np.zeros([4, 8])

    # POP_DIV_PREMONO_H
    D[0, IDX_SRC_Y] = -1. / dist_src_mono
    D[0, IDX_MONO_X] = -c_th_m / dist_src_mono
    D[0, IDX_MONO_Y] = s_th_m / dist_src_mono

    # POP_DIV_PRESAMPLE_H
    D[1, IDX_MONO_X] = c_th_m / dist_mono_sample
    D[1, IDX_MONO_Y] = s_th_m / dist_mono_sample
    D[1, IDX_SAMPLE_X] = s_th_s / dist_mono_sample
    D[1, IDX_SAMPLE_Y] = c_th_s / dist_mono_sample

    # POP_DIV_PREMONO_V
    D[2, IDX_SRC_Z] = -1. / dist_src_mono
    D[2, IDX_MONO_Z] = 1. / dist_src_mono

    # POP_DIV_PRESAMPLE_V
    D[3, IDX_MONO_Z] = -1. / dist_mono_sample
    D[3, IDX_SAMPLE_Z] = 1. / dist_mono_sample


    # T matrix, [pop75], Appendix 2
    T = np.zeros([2, 8])

    # horizontal
    T[0, IDX_SRC_Y] = 0.5 * D[0, IDX_SRC_Y]
    T[0, IDX_MONO_X] = 0.5 * (D[0, IDX_MONO_X] + D[1, IDX_MONO_X])
    T[0, IDX_MONO_Y] = 0.5 * (D[0, IDX_MONO_Y] + D[1, IDX_MONO_Y]) - inv_mono_curvh
    T[0, IDX_SAMPLE_X] = 0.5 * D[1, IDX_SAMPLE_X]
    T[0, IDX_SAMPLE_Y] = 0.5 * D[1, IDX_SAMPLE_Y]

    # vertical
    T[1, IDX_SRC_Z] = 0.5 * D[2, IDX_SRC_Z] / s_th_m
    T[1, IDX_MONO_Z] = 0.5 * (D[2, IDX_MONO_Z] - D[3, IDX_MONO_Z]) / s_th_m - inv_mono_curvv
    T[1, IDX_SAMPLE_Z] = -0.5 * D[3, IDX_SAMPLE_Z] / s_th_m



#
# resolution algorithm
#
def calc_pop(param):
    twotheta = param["twotheta"] * param["sample_sense"]
    thetam = param["thetam"] * param["mono_sense"]
    thetaa = param["thetaa"] * param["ana_sense"]
    ki_Q = param["angle_ki_Q"] * param["sample_sense"]
    kf_Q = param["angle_kf_Q"] * param["sample_sense"]

    ki = param["ki"]
    kf = param["kf"]
    E = param["E"]
    Q = param["Q"]


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


    lam = helpers.k2lam(ki)

    coll_h_pre_mono = param["coll_h_pre_mono"]
    coll_v_pre_mono = param["coll_v_pre_mono"]

    if param["use_guide"]:
        coll_h_pre_mono = lam*param["guide_div_h"]
        coll_v_pre_mono = lam*param["guide_div_v"]



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



    #
    # TODO
    #


    return res
