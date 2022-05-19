/**
 * r0 resolution function normalisation/intensity factor
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2019
 * @license GPLv2
 *
 * @desc see:
 *	[ch73] N. J. Chesser and J. D. Axe, Acta Cryst. A 29, 160 (1973)
 *	[mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984)
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#ifndef __TAKIN_R0_H__
#define __TAKIN_R0_H__

#include "defs.h"
#include "tlibs/phys/neutrons.h"


/**
 * R0 factor from formula (2) in [ch73]
 */
extern t_real_reso chess_R0(
	tl::t_wavenumber_si<t_real_reso> ki, tl::t_wavenumber_si<t_real_reso> kf,
        tl::t_angle_si<t_real_reso> theta_m, tl::t_angle_si<t_real_reso> theta_a,
	tl::t_angle_si<t_real_reso> twotheta_s,
        tl::t_angle_si<t_real_reso> mos_m, tl::t_angle_si<t_real_reso> mos_a,
	tl::t_angle_si<t_real_reso> coll_pre_mono_v, tl::t_angle_si<t_real_reso> coll_post_ana_v,
        t_real_reso refl_m, t_real_reso refl_a);


/**
 * general R0 normalisation factor from [mit84], equ. A.57
 */
template<class t_real = double>
t_real mitch_R0(t_real dmono_refl, t_real dana_effic,
	t_real dKiVol, t_real dKfVol, t_real dResVol,
	bool bNormToResVol=0)
{
	t_real dR0 = dmono_refl * dKiVol * dana_effic * dKfVol;

	// not needed for MC simulations, because the gaussian generated
	// with std::normal_distribution is already normalised
	// see: tools/test/tst_norm.cpp
	if(bNormToResVol)
		dR0 /= (dResVol * tl::get_pi<t_real>() * t_real{3});

	return dR0;
}

#endif
