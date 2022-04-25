#
# tests the resolution calculation
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license GPLv2
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
import reso
import helpers
import pop
import eck

np.set_printoptions(floatmode = "fixed",  precision = 4)


verbose = True

cm2A = 1e8
min2rad = 1./ 60. / 180.*np.pi
rad2deg = 180. / np.pi

d_mono = 3.355
d_ana = 3.355

ki = 1.4
kf = 1.5
Q = 1.5
E = helpers.get_E(ki, kf)

sc_senses = [ 1., -1., 1.]

params = {
    # scattering triangle
    "ki" : ki, "kf" : kf, "E" : E, "Q" : Q,

    # angles
    "twotheta" : helpers.get_scattering_angle(ki, kf, Q),
    "thetam" : helpers.get_mono_angle(ki, d_mono),
    "thetaa" : helpers.get_mono_angle(kf, d_ana),
    "angle_ki_Q" : helpers.get_angle_ki_Q(ki, kf, Q),
    "angle_kf_Q" : helpers.get_angle_kf_Q(ki, kf, Q),

     # scattering senses
    "mono_sense" : sc_senses[0],
    "sample_sense" : sc_senses[1],
    "ana_sense" : sc_senses[2],

    # distances
    "dist_src_mono" : 100. * cm2A,
    "dist_mono_sample" : 100. * cm2A,
    "dist_sample_ana" : 100. * cm2A,
    "dist_ana_det" : 100. * cm2A,

    # component sizes
    "src_w" : 10. * cm2A,
    "src_h" : 10. * cm2A,
    "mono_d" : 0.2 * cm2A,
    "mono_w" : 10. * cm2A,
    "mono_h" : 10. * cm2A,
    "sample_d" : 1. * cm2A,
    "sample_w" : 1. * cm2A,
    "sample_h" : 1. * cm2A,
    "det_w" : 10. * cm2A,
    "det_h" : 10. * cm2A,
    "ana_d" : 0.2 * cm2A,
    "ana_w" : 10. * cm2A,
    "ana_h" : 10. * cm2A,

    # focusing
    "mono_curvh" : 0.,
    "mono_curvv" : 0.,
    "ana_curvh" : 0.,
    "ana_curvv" : 0.,
    "mono_is_optimally_curved_h" : False,
    "mono_is_optimally_curved_v" : False,
    "ana_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_v" : False,
    "mono_is_curved_h" : False,
    "mono_is_curved_v" : False,
    "ana_is_curved_h" : False,
    "ana_is_curved_v" : False,

    # collimation
    "coll_h_pre_mono" : 9999. *min2rad,
    "coll_v_pre_mono" : 9999. *min2rad,
    "coll_h_pre_sample" : 30. *min2rad,
    "coll_v_pre_sample" : 9999. *min2rad,
    "coll_h_post_sample" : 30. *min2rad,
    "coll_v_post_sample" : 9999. *min2rad,
    "coll_h_post_ana" : 9999. *min2rad,
    "coll_v_post_ana" : 9999. *min2rad,

    # guide
    "use_guide" : False,
    "guide_div_h" : 9999. *min2rad,
    "guide_div_v" : 9999. *min2rad,

    # mosaics
    "mono_mosaic" : 60. *min2rad,
    "mono_mosaic_v" : 60. *min2rad,
    "sample_mosaic" : 5. *min2rad,
    "sample_mosaic_v" : 5. *min2rad,
    "ana_mosaic" : 60. *min2rad,
    "ana_mosaic_v" : 60. *min2rad,

    # crystal reflectivities
    # TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * cm2A,
    "pos_y" : 0. * cm2A,
    "pos_z" : 0. * cm2A,

    "kf_vert" : False,
}


# calculate resolution ellipsoid
#res = eck.calc(params)
res = pop.calc(params)

if not res["ok"]:
    print("RESOLUTION CALCULATION FAILED!")
    exit(-1)

if verbose:
    print("2theta = %g, thetam = %g, thetaa = %g, ki_Q = %g, kf_Q = %g\n" %
        (params["twotheta"]*rad2deg, params["thetam"]*rad2deg, params["thetaa"]*rad2deg,
        params["angle_ki_Q"]*rad2deg, params["angle_kf_Q"]*rad2deg))
    print("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
    print("Resolution matrix:\n%s" % res["reso"])
    print("Resolution vector: %s" % res["reso_v"])
    print("Resolution scalar: %g" % res["reso_s"])


# describe and plot ellipses
ellipses = reso.calc_ellipses(res["reso"], verbose)
reso.plot_ellipses(ellipses, verbose)
