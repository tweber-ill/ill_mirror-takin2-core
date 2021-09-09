/**
 * 3d unit cell drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2016
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

#ifndef __TAZ_REAL_3D__
#define __TAZ_REAL_3D__

#include <QDialog>
#include <QStatusBar>

#include <memory>

#include "libs/plotgl.h"
#include "tlibs/math/linalg.h"
#include "libs/spacegroups/latticehelper.h"
#include "tlibs/phys/bz.h"
#include "libs/globals.h"


class Real3DDlg : public QDialog
{Q_OBJECT
protected:
	QSettings *m_pSettings = nullptr;
	QStatusBar *m_pStatus = nullptr;
	std::unique_ptr<PlotGl> m_pPlot;

public:
	Real3DDlg(QWidget* pParent, QSettings* = 0);
	virtual ~Real3DDlg();

	void CalcPeaks(const tl::Brillouin3D<t_real_glob>& ws,
		const xtl::LatticeCommon<t_real_glob>& realcommon);

protected:
	virtual void hideEvent(QHideEvent*) override;
	virtual void showEvent(QShowEvent*) override;
	virtual void closeEvent(QCloseEvent*) override;

	virtual void keyPressEvent(QKeyEvent*) override;
};

#endif
