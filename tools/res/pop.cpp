/**
 * popovici calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the file rc_popma.m of the
 *		rescal5 package by Zinkin, McMorrow, Tennant, Farhi, and Wildes (ca. 1995-2007):
 *		http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 * @desc see: - [pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
 *            - [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
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

#include "pop.h"
#include "r0.h"
#include "helper.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/math.h"

#include <string>
#include <iostream>


typedef t_real_reso t_real;
typedef ublas::matrix<t_real> t_mat;
typedef ublas::vector<t_real> t_vec;

using angle = tl::t_angle_si<t_real>;
using wavenumber = tl::t_wavenumber_si<t_real>;
using energy = tl::t_energy_si<t_real>;
using length = tl::t_length_si<t_real>;
using inv_length = tl::t_length_inverse_si<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();
static const auto meV = tl::get_one_meV<t_real>();
static const auto cm = tl::get_one_centimeter<t_real>();
static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();

enum PopPosIdx : std::size_t
{
	POP_SRC_Y = 0, POP_SRC_Z,                 // w, h
	POP_MONO_X, POP_MONO_Y, POP_MONO_Z,       // d, w, h
	POP_SAMPLE_X, POP_SAMPLE_Y, POP_SAMPLE_Z, // perp Q, para Q, h
	POP_ANA_X, POP_ANA_Y, POP_ANA_Z,          // d, w, h
	POP_DET_Y, POP_DET_Z,                     // w, h

	POP_NUM_POS
};

enum PopMosaicIdx : std::size_t
{
	POP_MONO_MOSAIC_Y = 0, POP_MONO_MOSAIC_Z,
	POP_ANA_MOSAIC_Y, POP_ANA_MOSAIC_Z,

	POP_NUM_MOSAIC
};

enum PopDivIdx : std::size_t
{
	POP_DIV_PREMONO_H = 0, POP_DIV_PRESAMPLE_H,
	POP_DIV_PREMONO_V, POP_DIV_PRESAMPLE_V,
	POP_DIV_POSTSAMPLE_H, POP_DIV_POSTANA_H,
	POP_DIV_POSTSAMPLE_V, POP_DIV_POSTANA_V,

	POP_NUM_DIV
};

enum PopKiKfIdx : std::size_t
{
	POP_KI_X = 0, POP_KI_Y, POP_KI_Z,
	POP_KF_X, POP_KF_Y, POP_KF_Z,

	POP_NUM_KIKF
};


ResoResults calc_pop(const PopParams& pop)
{
	ResoResults res;

	res.Q_avg.resize(4);
	res.Q_avg[0] = pop.Q * angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = pop.E / meV;


	length lam = tl::k2lam(pop.ki);
	angle twotheta = pop.twotheta;
	angle thetaa = pop.thetaa * pop.dana_sense;
	angle thetam = pop.thetam * pop.dmono_sense;
	angle ki_Q = pop.angle_ki_Q;
	angle kf_Q = pop.angle_kf_Q;
	//kf_Q = ki_Q + twotheta;

	twotheta *= pop.dsample_sense;
	ki_Q *= pop.dsample_sense;
	kf_Q *= pop.dsample_sense;

	// B matrix, [pop75], Appendix 1 -> U matrix in CN
	t_mat B_trafo_QE = get_trafo_dkidkf_dQdE(ki_Q, kf_Q, pop.ki, pop.kf);
	B_trafo_QE.resize(4, POP_NUM_KIKF, true);

	angle coll_h_pre_mono = pop.coll_h_pre_mono;
	angle coll_v_pre_mono = pop.coll_v_pre_mono;

	if(pop.bGuide)
	{
		coll_h_pre_mono = lam*(pop.guide_div_h/angs);
		coll_v_pre_mono = lam*(pop.guide_div_v/angs);
	}


	// collimator covariance matrix G, [pop75], Appendix 1
	t_mat G_collis = tl::zero_matrix(POP_NUM_DIV, POP_NUM_DIV);

	G_collis(POP_DIV_PREMONO_H, POP_DIV_PREMONO_H) = 
		t_real(1) / (coll_h_pre_mono*coll_h_pre_mono /rads/rads);
	G_collis(POP_DIV_PRESAMPLE_H, POP_DIV_PRESAMPLE_H) = 
		t_real(1) / (pop.coll_h_pre_sample*pop.coll_h_pre_sample /rads/rads);

	G_collis(POP_DIV_POSTSAMPLE_H, POP_DIV_POSTSAMPLE_H) = 
		t_real(1) / (pop.coll_h_post_sample*pop.coll_h_post_sample /rads/rads);
	G_collis(POP_DIV_POSTANA_H, POP_DIV_POSTANA_H) = 
		t_real(1) / (pop.coll_h_post_ana*pop.coll_h_post_ana /rads/rads);

	G_collis(POP_DIV_PREMONO_V, POP_DIV_PREMONO_V) = 
		t_real(1) / (coll_v_pre_mono*coll_v_pre_mono /rads/rads);
	G_collis(POP_DIV_PRESAMPLE_V, POP_DIV_PRESAMPLE_V) = 
		t_real(1) / (pop.coll_v_pre_sample*pop.coll_v_pre_sample /rads/rads);

	G_collis(POP_DIV_POSTSAMPLE_V, POP_DIV_POSTSAMPLE_V) = 
		t_real(1) / (pop.coll_v_post_sample*pop.coll_v_post_sample /rads/rads);
	G_collis(POP_DIV_POSTANA_V, POP_DIV_POSTANA_V) = 
		t_real(1) / (pop.coll_v_post_ana*pop.coll_v_post_ana /rads/rads);


	const angle mono_mosaic_z = pop.mono_mosaic;
	const angle ana_mosaic_z = pop.ana_mosaic;
	const angle sample_mosaic_z = pop.sample_mosaic;

	// crystal mosaic covariance matrix F, [pop75], Appendix 1
	t_mat F_mosaics = tl::zero_matrix(POP_NUM_MOSAIC, POP_NUM_MOSAIC);
	F_mosaics(POP_MONO_MOSAIC_Y, POP_MONO_MOSAIC_Y) =
		t_real(1)/(pop.mono_mosaic*pop.mono_mosaic /rads/rads);
	F_mosaics(POP_MONO_MOSAIC_Z, POP_MONO_MOSAIC_Z) =
		t_real(1)/(mono_mosaic_z*mono_mosaic_z /rads/rads);
	F_mosaics(POP_ANA_MOSAIC_Y, POP_ANA_MOSAIC_Y) =
		t_real(1)/(pop.ana_mosaic*pop.ana_mosaic /rads/rads);
	F_mosaics(POP_ANA_MOSAIC_Z, POP_ANA_MOSAIC_Z) =
		t_real(1)/(ana_mosaic_z*ana_mosaic_z /rads/rads);


	// A matrix, [pop75], Appendix 1
	t_mat A_div_kikf_trafo = ublas::zero_matrix<t_real>(POP_NUM_KIKF, POP_NUM_DIV);
	A_div_kikf_trafo(POP_KI_X, POP_DIV_PREMONO_H) = t_real(0.5) * pop.ki*angs * 
		units::cos(thetam)/units::sin(thetam);
	A_div_kikf_trafo(POP_KI_X, POP_DIV_PRESAMPLE_H) = t_real(-0.5) * pop.ki*angs * 
		units::cos(thetam)/units::sin(thetam);
	A_div_kikf_trafo(POP_KI_Y, POP_DIV_PRESAMPLE_H) = pop.ki * angs;
	A_div_kikf_trafo(POP_KI_Z, POP_DIV_PRESAMPLE_V) = pop.ki * angs;

	A_div_kikf_trafo(POP_KF_X, POP_DIV_POSTSAMPLE_H) = t_real(0.5) * pop.kf*angs * 
		units::cos(thetaa)/units::sin(thetaa);
	A_div_kikf_trafo(POP_KF_X, POP_DIV_POSTANA_H) = t_real(-0.5) * pop.kf*angs * 
		units::cos(thetaa)/units::sin(thetaa);
	A_div_kikf_trafo(POP_KF_Y, POP_DIV_POSTSAMPLE_H) = pop.kf * angs;
	A_div_kikf_trafo(POP_KF_Z, POP_DIV_POSTSAMPLE_V) = pop.kf * angs;


	// covariance matrix of component geometries, S, [pop75], Appendices 2 and 3
	t_real dMultSrc = pop.bSrcRect ? 1./12. : 1./16.;
	t_real dMultSample = pop.bSampleCub ? 1./12. : 1./16.;
	t_real dMultDet = pop.bDetRect ? 1./12. : 1./16.;

	t_mat SI_geo = tl::zero_matrix(POP_NUM_POS, POP_NUM_POS);
	SI_geo(POP_SRC_Y, POP_SRC_Y) = dMultSrc * pop.src_w*pop.src_w /cm/cm;
	SI_geo(POP_SRC_Z, POP_SRC_Z) = dMultSrc * pop.src_h*pop.src_h /cm/cm;

	SI_geo(POP_MONO_X, POP_MONO_X) = t_real(1./12.) * pop.mono_thick*pop.mono_thick /cm/cm;
	SI_geo(POP_MONO_Y, POP_MONO_Y) = t_real(1./12.) * pop.mono_w*pop.mono_w /cm/cm;
	SI_geo(POP_MONO_Z, POP_MONO_Z) = t_real(1./12.) * pop.mono_h*pop.mono_h /cm/cm;

	SI_geo(POP_SAMPLE_X, POP_SAMPLE_X) = dMultSample * pop.sample_w_perpq*pop.sample_w_perpq /cm/cm;
	SI_geo(POP_SAMPLE_Y, POP_SAMPLE_Y) = dMultSample * pop.sample_w_q*pop.sample_w_q /cm/cm;
	SI_geo(POP_SAMPLE_Z, POP_SAMPLE_Z) = t_real(1./12.) * pop.sample_h*pop.sample_h /cm/cm;

	SI_geo(POP_ANA_X, POP_ANA_X) = t_real(1./12.) * pop.ana_thick*pop.ana_thick /cm/cm;
	SI_geo(POP_ANA_Y, POP_ANA_Y) = t_real(1./12.) * pop.ana_w*pop.ana_w /cm/cm;
	SI_geo(POP_ANA_Z, POP_ANA_Z) = t_real(1./12.) * pop.ana_h*pop.ana_h /cm/cm;

	SI_geo(POP_DET_Y, POP_DET_Y) = dMultDet * pop.det_w*pop.det_w /cm/cm;
	SI_geo(POP_DET_Z, POP_DET_Z) = dMultDet * pop.det_h*pop.det_h /cm/cm;

	SI_geo *= sig2fwhm*sig2fwhm;

	t_mat S_geo;
	if(!tl::inverse(SI_geo, S_geo))
	{
		res.bOk = false;
		res.strErr = "S matrix cannot be inverted.";
		return res;
	}


	// --------------------------------------------------------------------
	// mono/ana focus
	length mono_curvh = pop.mono_curvh, mono_curvv = pop.mono_curvv;
	length ana_curvh = pop.ana_curvh, ana_curvv = pop.ana_curvv;

	if(pop.bMonoIsOptimallyCurvedH)
		mono_curvh = tl::foc_curv(pop.dist_src_mono, pop.dist_mono_sample, units::abs(t_real(2)*thetam), false);
	if(pop.bAnaIsOptimallyCurvedH)
		ana_curvh = tl::foc_curv(pop.dist_sample_ana, pop.dist_ana_det, units::abs(t_real(2)*thetaa), false);
	if(pop.bMonoIsOptimallyCurvedV)
		mono_curvv = tl::foc_curv(pop.dist_src_mono, pop.dist_mono_sample, units::abs(t_real(2)*thetam), true);
	if(pop.bAnaIsOptimallyCurvedV)
		ana_curvv = tl::foc_curv(pop.dist_sample_ana, pop.dist_ana_det, units::abs(t_real(2)*thetaa), true);

	mono_curvh *= pop.dmono_sense;
	ana_curvh *= pop.dana_sense;
	mono_curvv *= pop.dmono_sense;
	ana_curvv *= pop.dana_sense;

	inv_length inv_mono_curvh = pop.bMonoIsCurvedH ? t_real(1)/mono_curvh : t_real(0)/cm;
	inv_length inv_ana_curvh = pop.bAnaIsCurvedH ? t_real(1)/ana_curvh : t_real(0)/cm;
	inv_length inv_mono_curvv = pop.bMonoIsCurvedV ? t_real(1)/mono_curvv : t_real(0)/cm;
	inv_length inv_ana_curvv = pop.bAnaIsCurvedV ? t_real(1)/ana_curvv : t_real(0)/cm;

	const auto tupScFact = get_scatter_factors(pop.flags, pop.thetam, pop.ki, pop.thetaa, pop.kf);

	t_real dmono_refl = pop.dmono_refl * std::get<0>(tupScFact);
	t_real dana_effic = pop.dana_effic * std::get<1>(tupScFact);
	if(pop.mono_refl_curve) dmono_refl *= (*pop.mono_refl_curve)(pop.ki);
	if(pop.ana_effic_curve) dana_effic *= (*pop.ana_effic_curve)(pop.kf);
	t_real dxsec = std::get<2>(tupScFact);
	// --------------------------------------------------------------------


	// T matrix to transform the mosaic cov. matrix, [pop75], Appendix 2
	t_mat T_mosaic_trafo = ublas::zero_matrix<t_real>(POP_NUM_MOSAIC, POP_NUM_POS);
	T_mosaic_trafo(POP_MONO_MOSAIC_Y, POP_SRC_Y) = t_real(-0.5) / (pop.dist_src_mono / cm);
	T_mosaic_trafo(POP_MONO_MOSAIC_Y, POP_MONO_X) = t_real(0.5) * units::cos(thetam) * (
		t_real(1)/(pop.dist_mono_sample/cm) -
		t_real(1)/(pop.dist_src_mono/cm));
	T_mosaic_trafo(POP_MONO_MOSAIC_Y, POP_MONO_Y) = t_real(0.5) * units::sin(thetam) * (
		t_real(1)/(pop.dist_src_mono/cm) +
		t_real(1)/(pop.dist_mono_sample/cm) -
		t_real(2)*inv_mono_curvh*cm / units::sin(thetam));
	T_mosaic_trafo(POP_MONO_MOSAIC_Y, POP_SAMPLE_X) = t_real(0.5) *
		units::sin(t_real(0.5)*twotheta) / (pop.dist_mono_sample/cm);
	T_mosaic_trafo(POP_MONO_MOSAIC_Y, POP_SAMPLE_Y) = t_real(0.5) *
		units::cos(t_real(0.5)*twotheta) / (pop.dist_mono_sample/cm);

	T_mosaic_trafo(POP_ANA_MOSAIC_Y, POP_DET_Y) = t_real(0.5) / (pop.dist_ana_det / cm);
	T_mosaic_trafo(POP_ANA_MOSAIC_Y, POP_ANA_X) = t_real(-0.5) * units::cos(thetaa) * (
		t_real(1)/(pop.dist_sample_ana/cm) -
		t_real(1)/(pop.dist_ana_det/cm));
	T_mosaic_trafo(POP_ANA_MOSAIC_Y, POP_ANA_Y) = t_real(0.5) * units::sin(thetaa) * (
		t_real(1)/(pop.dist_ana_det/cm) +
		t_real(1)/(pop.dist_sample_ana/cm) -
		t_real(2)*inv_ana_curvh*cm / units::sin(thetaa));
	T_mosaic_trafo(POP_ANA_MOSAIC_Y, POP_SAMPLE_X) = t_real(0.5) *
		units::sin(t_real(0.5)*twotheta) / (pop.dist_sample_ana/cm);
	T_mosaic_trafo(POP_ANA_MOSAIC_Y, POP_SAMPLE_Y) = t_real(-0.5) *
		units::cos(t_real(0.5)*twotheta) / (pop.dist_sample_ana/cm);

	T_mosaic_trafo(POP_MONO_MOSAIC_Z, POP_SRC_Z) = t_real(-0.5) / (
		pop.dist_src_mono/cm * units::sin(thetam));
	T_mosaic_trafo(POP_MONO_MOSAIC_Z, POP_MONO_Z) = t_real(0.5) * (
		t_real(1)/(pop.dist_src_mono/cm * units::sin(thetam)) +
		t_real(1)/(pop.dist_mono_sample/cm * units::sin(thetam)) -
		t_real(2)*inv_mono_curvv*cm);
	T_mosaic_trafo(POP_MONO_MOSAIC_Z, POP_SAMPLE_Z) = t_real(-0.5) / (
		pop.dist_mono_sample/cm * units::sin(thetam));

	T_mosaic_trafo(POP_ANA_MOSAIC_Z, POP_DET_Z) = t_real(-0.5) / (
		pop.dist_ana_det/cm * units::sin(thetaa));
	T_mosaic_trafo(POP_ANA_MOSAIC_Z, POP_ANA_Z) = t_real(0.5) * (
		t_real(1)/(pop.dist_ana_det/cm * units::sin(thetaa)) +
		t_real(1)/(pop.dist_sample_ana/cm * units::sin(thetaa)) -
		t_real(2)*inv_ana_curvv*cm);
	T_mosaic_trafo(POP_ANA_MOSAIC_Z, POP_SAMPLE_Z) = t_real(-0.5) / (
		pop.dist_sample_ana/cm * units::sin(thetaa));


	// D matrix to transform spatial to divergence variables, [pop75], Appendix 2
	t_mat D_geo_div_trafo = ublas::zero_matrix<t_real>(POP_NUM_DIV, POP_NUM_POS);
	D_geo_div_trafo(POP_DIV_PREMONO_H, POP_SRC_Y) = t_real(-1) / (pop.dist_src_mono/cm);
	D_geo_div_trafo(POP_DIV_PREMONO_H, POP_MONO_X) = -cos(thetam) / (pop.dist_src_mono/cm);
	D_geo_div_trafo(POP_DIV_PREMONO_H, POP_MONO_Y) = sin(thetam) / (pop.dist_src_mono/cm);

	D_geo_div_trafo(POP_DIV_POSTANA_H, POP_DET_Y) = t_real(1) / (pop.dist_ana_det/cm);
	D_geo_div_trafo(POP_DIV_POSTANA_H, POP_ANA_X) = cos(thetaa) / (pop.dist_ana_det/cm);
	D_geo_div_trafo(POP_DIV_POSTANA_H, POP_ANA_Y) = sin(thetaa) / (pop.dist_ana_det/cm);

	D_geo_div_trafo(POP_DIV_PRESAMPLE_H, POP_MONO_X) = cos(thetam) / (pop.dist_mono_sample/cm);
	D_geo_div_trafo(POP_DIV_PRESAMPLE_H, POP_MONO_Y) = sin(thetam) / (pop.dist_mono_sample/cm);
	D_geo_div_trafo(POP_DIV_PRESAMPLE_H, POP_SAMPLE_X) = sin(t_real(0.5)*twotheta) / (pop.dist_mono_sample/cm);
	D_geo_div_trafo(POP_DIV_PRESAMPLE_H, POP_SAMPLE_Y) = cos(t_real(0.5)*twotheta) / (pop.dist_mono_sample/cm);

	D_geo_div_trafo(POP_DIV_POSTSAMPLE_H, POP_ANA_X) = -cos(thetaa) / (pop.dist_sample_ana/cm);
	D_geo_div_trafo(POP_DIV_POSTSAMPLE_H, POP_ANA_Y) = sin(thetaa) / (pop.dist_sample_ana/cm);
	D_geo_div_trafo(POP_DIV_POSTSAMPLE_H, POP_SAMPLE_X) = sin(t_real(0.5)*twotheta) / (pop.dist_sample_ana/cm);
	D_geo_div_trafo(POP_DIV_POSTSAMPLE_H, POP_SAMPLE_Y) = -cos(t_real(0.5)*twotheta) / (pop.dist_sample_ana/cm);

	D_geo_div_trafo(POP_DIV_PREMONO_V, POP_SRC_Z) = t_real(-1) / (pop.dist_src_mono/cm);
	D_geo_div_trafo(POP_DIV_PREMONO_V, POP_MONO_Z) = t_real(1) / (pop.dist_src_mono/cm);

	D_geo_div_trafo(POP_DIV_PRESAMPLE_V, POP_MONO_Z) = t_real(-1) / (pop.dist_mono_sample/cm);
	D_geo_div_trafo(POP_DIV_PRESAMPLE_V, POP_SAMPLE_Z) = t_real(1) / (pop.dist_mono_sample/cm);

	D_geo_div_trafo(POP_DIV_POSTSAMPLE_V, POP_SAMPLE_Z) = t_real(-1) / (pop.dist_sample_ana/cm);
	D_geo_div_trafo(POP_DIV_POSTSAMPLE_V, POP_ANA_Z) = t_real(1) / (pop.dist_sample_ana/cm);

	D_geo_div_trafo(POP_DIV_POSTANA_V, POP_ANA_Z) = t_real(-1) / (pop.dist_ana_det/cm);
	D_geo_div_trafo(POP_DIV_POSTANA_V, POP_DET_Z) = t_real(1) / (pop.dist_ana_det/cm);


	// [pop75], equ. 20
	// [T] = 1/cm, [F] = 1/rad^2, [pop75], equ. 15
	t_mat K = S_geo + tl::transform(F_mosaics, T_mosaic_trafo, true);
	t_mat Ki;
	if(!tl::inverse(K, Ki))
	{
		res.bOk = false;
		res.strErr = "Matrix K cannot be inverted.";
		return res;
	}

	// [pop75], equ. 17
	t_mat Hi = tl::transform_inv(Ki, D_geo_div_trafo, true);
	t_mat H;
	if(!tl::inverse(Hi, H))
	{
		res.bOk = false;
		res.strErr = "Matrix H^(-1) cannot be inverted.";
		return res;
	}

	// [pop75], equ. 18
	t_mat H_G = H + G_collis;
	t_mat H_Gi;
	if(!tl::inverse(H_G, H_Gi))
	{
		res.bOk = false;
		res.strErr = "Matrix H+G cannot be inverted.";
		return res;
	}

	// [pop75], equ. 20
	t_mat BA = ublas::prod(B_trafo_QE, A_div_kikf_trafo);
	t_mat cov = tl::transform_inv(H_Gi, BA, true);
	cov(1,1) += pop.Q*pop.Q*angs*angs * pop.sample_mosaic*pop.sample_mosaic /rads/rads;
	cov(2,2) += pop.Q*pop.Q*angs*angs * sample_mosaic_z*sample_mosaic_z /rads/rads;

	if(!tl::inverse(cov, res.reso))
	{
		res.bOk = false;
		res.strErr = "Covariance matrix cannot be inverted.";
		return res;
	}


	// -------------------------------------------------------------------------


	res.reso *= sig2fwhm*sig2fwhm;
	res.reso_v = ublas::zero_vector<t_real>(4);
	res.reso_s = 0.;

	if(pop.dsample_sense < 0.)
	{
		// mirror Q_perp
		t_mat matMirror = tl::mirror_matrix<t_mat>(res.reso.size1(), 1);
		res.reso = tl::transform(res.reso, matMirror, true);
		res.reso_v[1] = -res.reso_v[1];
	}


	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	res.dR0 = 0.;
	if(pop.flags & CALC_R0)
	{
		const t_real pi = tl::get_pi<t_real>();

		// resolution volume, [pop75], equ. 13a & 16
		// [D] = 1/cm, [SI] = cm^2
		t_mat DSiDt = tl::transform_inv(SI_geo, D_geo_div_trafo, true);
		t_mat DSiDti;
		if(!tl::inverse(DSiDt, DSiDti))
		{
			res.bOk = false;
			res.strErr = "R0 factor cannot be calculated.";
			return res;
		}
		DSiDti += G_collis;

		t_real dDetS = tl::determinant(S_geo);
		t_real dDetF = tl::determinant(F_mosaics);
		t_real dDetK = tl::determinant(K);
		t_real dDetDSiDti = tl::determinant(DSiDti);

		// [pop75], equs. 13a & 16
		res.dR0 = dmono_refl*dana_effic * t_real((2.*pi)*(2.*pi)*(2.*pi)*(2.*pi));
		res.dR0 *= std::sqrt(dDetS*dDetF/(dDetK * dDetDSiDti));
		res.dR0 /= t_real(8.*pi*8.*pi) * units::sin(thetam)*units::sin(thetaa);
		res.dR0 *= dxsec;

		// rest of the prefactors, equ. 1 in [pop75], together with the mono and and ana reflectivities
		// (defining the resolution volume) these give the same correction as in [mit84] equ. A.57
		// NOTE: these factors are not needed, because the normalisation of the 4d gaussian distribution
		// is already taken care of in the MC step by the employed std::normal_distribution function
		//res.dR0 *= std::sqrt(std::abs(tl::determinant(res.reso))) / (2.*pi*2.*pi);
		// except for the (unimportant) prefactors this is the same as dividing by the resolution volume
		//res.dR0 /= res.dResVol * pi * t_real(3.);
	}

	// Bragg widths
	const std::vector<t_real> vecFwhms = calc_bragg_fwhms(res.reso);
	std::copy(vecFwhms.begin(), vecFwhms.end(), res.dBraggFWHMs);


	if(tl::is_nan_or_inf(res.dR0) || tl::is_nan_or_inf(res.reso))
	{
		res.strErr = "Invalid result.";
		res.bOk = false;
		return res;
	}

	res.bOk = true;
	return res;
}
