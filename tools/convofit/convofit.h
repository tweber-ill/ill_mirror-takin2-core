/**
 * Convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2
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

#ifndef __TAKIN_CONVOFIT_H__
#define __TAKIN_CONVOFIT_H__

#include "scan.h"
#include "model.h"
#include "tlibs/gfx/gnuplot.h"

#include <boost/signals2.hpp>
namespace sig = boost::signals2;


/**
 * set resolution parameters from scan file
 */
static inline void set_tasreso_params_from_scan(TASReso& reso, const Scan& sc)
{
	reso.SetLattice(sc.sample.a, sc.sample.b, sc.sample.c,
		sc.sample.alpha, sc.sample.beta, sc.sample.gamma,
		tl::make_vec({sc.plane.vec1[0], sc.plane.vec1[1], sc.plane.vec1[2]}),
		tl::make_vec({sc.plane.vec2[0], sc.plane.vec2[1], sc.plane.vec2[2]}));
	reso.SetKiFix(sc.bKiFixed);
	reso.SetKFix(sc.dKFix);
}


/**
 * set s(q,w) model parameters from scan file
 */
static inline void set_model_params_from_scan(SqwFuncModel& mod, const Scan& sc)
{
	mod.SetScanOrigin(sc.vecScanOrigin[0], sc.vecScanOrigin[1],
		sc.vecScanOrigin[2], sc.vecScanOrigin[3]);
	mod.SetScanDir(sc.vecScanDir[0], sc.vecScanDir[1],
		sc.vecScanDir[2], sc.vecScanDir[3]);

	std::pair<decltype(sc.vecX)::const_iterator, decltype(sc.vecX)::const_iterator> xminmax
		= std::minmax_element(sc.vecX.begin(), sc.vecX.end());
	mod.SetPrincipalScanAxisMinMax(t_real_mod(*xminmax.first), t_real_mod(*xminmax.second));

	mod.SetOtherParams(sc.dTemp, sc.dField);
}


// --------------------------------------------------------------------
// global command-line overrides
extern bool g_bVerbose;
extern bool g_bSkipFit;
extern bool g_bUseValuesFromModel;
extern unsigned int g_iNumNeutrons;
extern std::string g_strSetParams;
extern std::string g_strOutFileSuffix;
extern unsigned int g_iPlotPoints;
extern unsigned int g_iPlotSkipBegin;
extern unsigned int g_iPlotSkipEnd;
// --------------------------------------------------------------------


class Convofit
{
	public:
		using t_sigInitPlotter = sig::signal<void*(const std::string&)>;
		using t_sigDeinitPlotter = sig::signal<void(void*&)>;
		using t_sigPlot = sig::signal<void(void*, const char*, const char*,
			const tl::PlotObj<t_real_reso>&, const tl::PlotObj<t_real_reso>&,
			bool)>;

	protected:
		// user-supplied plotter
		void *m_pPlt = nullptr;

		t_sigInitPlotter m_sigInitPlotter;
		t_sigDeinitPlotter m_sigDeinitPlotter;
		t_sigPlot m_sigPlot;

	public:
		Convofit(bool bUseDefaultPlotter=1);
		~Convofit();

		bool run_job(const std::string& _strJob);

		void addsig_initplotter(const typename t_sigInitPlotter::slot_type& conn)
		{ m_sigInitPlotter.connect(conn); }
		void addsig_deinitplotter(const typename t_sigDeinitPlotter::slot_type& conn)
		{ m_sigDeinitPlotter.connect(conn); }
		void addsig_plot(const typename t_sigPlot::slot_type& conn)
		{ m_sigPlot.connect(conn); }
};

#endif
