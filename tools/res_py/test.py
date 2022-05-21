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
import helpers
import reso
import pop
import eck

np.set_printoptions(floatmode = "fixed",  precision = 4)


# settings
ki = 1.4
kf = 1.4
Q = 1.777

reso_method = "eck"    # "eck" or "pop"
verbose = True


#
# resolution parameters
# NOTE: not all parameters are used by all calculation backends
#
params = {
    # options
    "verbose" : verbose,

    # scattering triangle
    "ki" : ki,
    "kf" : kf,
    "E" : helpers.get_E(ki, kf),
    "Q" : Q,

    # d spacings
    "mono_xtal_d" : 3.355,
    "ana_xtal_d" : 3.355,

     # scattering senses
    "mono_sense" : -1.,
    "sample_sense" : 1.,
    "ana_sense" : -1.,
    "mirror_Qperp" : False,

    # distances
    "dist_src_mono" : 10. * helpers.cm2A,
    "dist_mono_sample" : 200. * helpers.cm2A,
    "dist_sample_ana" : 115. * helpers.cm2A,
    "dist_ana_det" : 85. * helpers.cm2A,

    # shapes
    "src_shape" : "rectangular",     # "rectangular" or "circular"
    "sample_shape" : "cylindrical",  # "cuboid" or "cylindrical"
    "det_shape" : "rectangular",     # "rectangular" or "circular"

    # component sizes
    "src_w" : 6. * helpers.cm2A,
    "src_h" : 12. * helpers.cm2A,
    "mono_d" : 0.15 * helpers.cm2A,
    "mono_w" : 12. * helpers.cm2A,
    "mono_h" : 8. * helpers.cm2A,
    "sample_d" : 1. * helpers.cm2A,
    "sample_w" : 1. * helpers.cm2A,
    "sample_h" : 1. * helpers.cm2A,
    "ana_d" : 0.3 * helpers.cm2A,
    "ana_w" : 12. * helpers.cm2A,
    "ana_h" : 8. * helpers.cm2A,
    "det_w" : 1.5 * helpers.cm2A,
    "det_h" : 5. * helpers.cm2A,

    # horizontal collimation
    "coll_h_pre_mono" : 30. * helpers.min2rad,
    "coll_h_pre_sample" : 30. * helpers.min2rad,
    "coll_h_post_sample" : 30. * helpers.min2rad,
    "coll_h_post_ana" : 9999. * helpers.min2rad,

    # vertical collimation
    "coll_v_pre_mono" : 9999. * helpers.min2rad,
    "coll_v_pre_sample" : 9999. * helpers.min2rad,
    "coll_v_post_sample" : 9999. * helpers.min2rad,
    "coll_v_post_ana" : 9999. * helpers.min2rad,

    # horizontal focusing
    "mono_curvh" : 0.,
    "ana_curvh" : 0.,
    "mono_is_curved_h" : False,
    "ana_is_curved_h" : False,
    "mono_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_h" : False,

    # vertical focusing
    "mono_curvv" : 0.,
    "ana_curvv" : 0.,
    "mono_is_curved_v" : True,
    "ana_is_curved_v" : False,
    "mono_is_optimally_curved_v" : True,
    "ana_is_optimally_curved_v" : False,

    # guide before monochromator
    "use_guide" : True,
    "guide_div_h" : 15. *helpers.min2rad,
    "guide_div_v" : 15. *helpers.min2rad,

    # horizontal mosaics
    "mono_mosaic" : 45. *helpers.min2rad,
    "sample_mosaic" : 5. *helpers.min2rad,
    "ana_mosaic" : 45. *helpers.min2rad,

    # vertical mosaics
    "mono_mosaic_v" : 45. *helpers.min2rad,
    "sample_mosaic_v" : 5. *helpers.min2rad,
    "ana_mosaic_v" : 45. *helpers.min2rad,

    # calculate R0 factor (not needed if only the ellipses are to be plotted)
    "calc_R0" : True,

    # crystal reflectivities
    # TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * helpers.cm2A,
    "pos_y" : 0. * helpers.cm2A,
    "pos_z" : 0. * helpers.cm2A,

    # vertical scattering in kf, keep "False" for normal TAS
    "kf_vert" : False,
}


# calculate resolution ellipsoid using the given backend
if reso_method == "eck":
    res = eck.calc(params)
elif reso_method == "pop":
    res = pop.calc(params)
else:
    raise "ResPy: Invalid resolution calculation method selected."


if not res["ok"]:
    print("RESOLUTION CALCULATION FAILED!")
    exit(-1)

if verbose:
    print("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
    print("Resolution matrix:\n%s" % res["reso"])
    print("Resolution vector: %s" % res["reso_v"])
    print("Resolution scalar: %g" % res["reso_s"])


# describe and plot ellipses
ellipses = reso.calc_ellipses(res["reso"], verbose)
reso.plot_ellipses(ellipses, verbose)
