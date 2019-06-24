/*
 * Scattering Triangle tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2014
 * @license GPLv2
 */

#ifndef __TASOPTIONS_H__
#define __TASOPTIONS_H__

#include "tlibs/math/linalg.h"
#include "libs/globals.h"
#include <QPointF>

namespace ublas = boost::numeric::ublas;


static inline ublas::vector<t_real_glob> qpoint_to_vec(const QPointF& pt)
{
	ublas::vector<t_real_glob> vec(2);
	vec[0] = t_real_glob(pt.x());
	vec[1] = t_real_glob(pt.y());

	return vec;
}

static inline QPointF vec_to_qpoint(const ublas::vector<t_real_glob>& vec)
{
	if(vec.size() < 2)
		return QPointF(0., 0.);

	return QPointF(vec[0], vec[1]);
}

struct TriangleOptions
{
	bool bChangedTheta = 0;
	bool bChangedTwoTheta = 0;
	bool bChangedAnaTwoTheta = 0;
	bool bChangedMonoTwoTheta = 0;

	bool bChangedMonoD = 0;
	bool bChangedAnaD = 0;

	bool bChangedAngleKiVec0 = 0;


	t_real_glob dTheta;
	t_real_glob dTwoTheta;
	t_real_glob dAnaTwoTheta;
	t_real_glob dMonoTwoTheta;

	t_real_glob dMonoD;
	t_real_glob dAnaD;

	t_real_glob dAngleKiVec0;


	bool operator==(const TriangleOptions& op) const
	{
		if(!tl::float_equal(dTheta, op.dTheta, g_dEps)) return false;
		if(!tl::float_equal(dTwoTheta, op.dTwoTheta, g_dEps)) return false;
		if(!tl::float_equal(dAnaTwoTheta, op.dAnaTwoTheta, g_dEps)) return false;
		if(!tl::float_equal(dMonoTwoTheta, op.dMonoTwoTheta, g_dEps)) return false;
		if(!tl::float_equal(dMonoD, op.dMonoD, g_dEps)) return false;
		if(!tl::float_equal(dAnaD, op.dAnaD, g_dEps)) return false;
		if(!tl::float_equal(dAngleKiVec0, op.dAngleKiVec0, g_dEps)) return false;
		return true;
	}

	bool IsAnythingChanged() const
	{
		return bChangedTheta || bChangedTwoTheta || bChangedAnaTwoTheta ||
			bChangedMonoTwoTheta || bChangedMonoD || bChangedAnaD ||
			bChangedAngleKiVec0;
	}

	void clear()
	{
		bChangedTheta = bChangedTwoTheta = bChangedAnaTwoTheta =
			bChangedMonoTwoTheta = bChangedMonoD = bChangedAnaD =
			bChangedAngleKiVec0 = 0;

		dTheta = dTwoTheta = dAnaTwoTheta = dMonoTwoTheta =
			dMonoD = dAnaD = dAngleKiVec0 = t_real_glob(0);
	}
};

struct CrystalOptions
{
	bool bChangedLattice = 0;
	bool bChangedLatticeAngles = 0;
	bool bChangedSpacegroup = 0;
	bool bChangedPlane1 = 0;
	bool bChangedPlane2 = 0;
	bool bChangedSampleName = 0;

	t_real_glob dLattice[3];
	t_real_glob dLatticeAngles[3];
	std::string strSpacegroup;
	t_real_glob dPlane1[3];
	t_real_glob dPlane2[3];

	std::string strSampleName;


	bool IsAnythingChanged() const
	{
		return bChangedLattice || bChangedLatticeAngles ||
			bChangedSpacegroup || bChangedPlane1 || bChangedPlane2 ||
			bChangedSampleName;
	}

	void clear()	// struct is no POD, cannot use memset
	{
		bChangedLattice = bChangedLatticeAngles =
			bChangedSpacegroup = bChangedPlane1 =
			bChangedPlane2 = bChangedSampleName = 0;

		strSpacegroup = strSampleName = "";

		for(char i=0; i<3; ++i)
		{
			dLattice[i] = dLatticeAngles[i] = dPlane1[i] = dPlane2[i] = t_real_glob(0);
		}
	}
};


#endif
