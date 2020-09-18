/**
 * Loads tabulated spacegroups
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
 * @license GPLv2
 */

#ifndef __SG_TAB_H__
#define __SG_TAB_H__

#include "crystalsys.h"
#include "sghelper.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"

#include <map>
#include <mutex>


namespace xtl {

namespace ublas = boost::numeric::ublas;


template<class t_real = double>
class SpaceGroup
{
public:
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;

protected:
	unsigned int m_iNr = 0;
	std::string m_strName;
	std::string m_strLaue, m_strPoint;
	CrystalSystem m_crystalsys = CRYS_NOT_SET;
	std::string m_strCrystalSysName;

	std::vector<t_mat> m_vecTrafos, m_vecCentringTrafos;
	std::vector<unsigned int> m_vecInvTrafos, m_vecPrimTrafos,
		m_vecCenterTrafos, m_vecTrans;

public:
	SpaceGroup() = default;
	~SpaceGroup() = default;

	SpaceGroup(const SpaceGroup& sg)
		: m_iNr(sg.m_iNr), m_strName(sg.m_strName), m_strLaue(sg.m_strLaue),
		m_strPoint(sg.m_strPoint), m_crystalsys(sg.m_crystalsys),
		m_strCrystalSysName(sg.m_strCrystalSysName), m_vecTrafos(sg.m_vecTrafos),
		m_vecCentringTrafos(sg.m_vecCentringTrafos),
		m_vecInvTrafos(sg.m_vecInvTrafos), m_vecPrimTrafos(sg.m_vecPrimTrafos),
		m_vecCenterTrafos(sg.m_vecCenterTrafos), m_vecTrans(sg.m_vecTrans)
	{}

	SpaceGroup(SpaceGroup&& sg)
		: m_iNr(sg.m_iNr), m_strName(std::move(sg.m_strName)), m_strLaue(std::move(sg.m_strLaue)),
		m_strPoint(std::move(sg.m_strPoint)), m_crystalsys(std::move(sg.m_crystalsys)),
		m_strCrystalSysName(std::move(sg.m_strCrystalSysName)), m_vecTrafos(std::move(sg.m_vecTrafos)),
		m_vecCentringTrafos(sg.m_vecCentringTrafos),
		m_vecInvTrafos(std::move(sg.m_vecInvTrafos)), m_vecPrimTrafos(std::move(sg.m_vecPrimTrafos)),
		m_vecCenterTrafos(std::move(sg.m_vecCenterTrafos)), m_vecTrans(std::move(sg.m_vecTrans))
	{}


	/**
	 * reflection systematically allowed? if not, return index of trafo which forbids it
	 */
	bool HasReflection(int h, int k, int l, std::size_t* pTrafoIdx=nullptr) const
	{
		std::pair<bool, std::size_t> pair =
			is_reflection_allowed<std::vector, t_mat, t_vec>
				(h,k,l, m_vecTrafos);

		if(pTrafoIdx)
			*pTrafoIdx = pair.second;
		return pair.first;
	}


	/*bool HasReflection2(int h, int k, int l) const
	{
		return is_reflection_allowed_2<std::vector, t_mat, t_vec>
		(h,k,l, m_vecTrafos);
	}*/


	/**
	 * is reflection allowed based only on centring?
	 */
	bool HasGenReflection(int h, int k, int l, std::size_t* pTrafoIdx=nullptr) const
	{
		bool bAllowed = 1;

		if(m_vecCentringTrafos.size()) // calculate from space group
		{
			std::pair<bool, std::size_t> pair =
				is_reflection_allowed<std::vector, t_mat, t_vec>
					(h,k,l, m_vecCentringTrafos);
			bAllowed = pair.first;
			if(pTrafoIdx)
				*pTrafoIdx = pair.second;
		}
		else	// use equivalent pre-calculated rules alternatively
		{
			bAllowed = is_centering_reflection_allowed<int>(GetName(), h,k,l);
			if(pTrafoIdx)
				*pTrafoIdx = 0;
		}

		return bAllowed;
	}


	void SetNr(unsigned int iNr)
	{
		m_iNr = iNr;
		SetCrystalSystem(get_crystal_system_from_space_group(int(m_iNr)));
	}
	unsigned int GetNr() const { return m_iNr; }


	void SetName(const std::string& str)
	{
		m_strName = str;
		m_strPoint = get_pointgroup(m_strName);
	}
	const std::string& GetName() const { return m_strName; }


	void SetCrystalSystem(const CrystalSystem& crys)
	{
		if(m_crystalsys == CRYS_NOT_SET)
		{
			m_crystalsys = crys;
			m_strCrystalSysName = get_crystal_system_name(m_crystalsys);
		}
	}
	CrystalSystem GetCrystalSystem() const { return m_crystalsys; }
	const std::string& GetCrystalSystemName() const { return m_strCrystalSysName; }


	void SetLaueGroup(const std::string& str)
	{
		m_strLaue = str;
		SetCrystalSystem(get_crystal_system_from_laue_group(m_strLaue.c_str()));
	}
	const std::string& GetLaueGroup() const { return m_strLaue; }
	const std::string& GetPointGroup() const { return m_strPoint; }


	void SetTrafos(std::vector<t_mat>&& vecTrafos) { m_vecTrafos = std::move(vecTrafos); }
	void SetTrafos(const std::vector<t_mat>& vecTrafos) { m_vecTrafos = vecTrafos; }
	const std::vector<t_mat>& GetTrafos() const { return m_vecTrafos; }

	void SetInvTrafos(std::vector<unsigned int>&& vecTrafos) { m_vecInvTrafos = std::move(vecTrafos); }
	void SetInvTrafos(const std::vector<unsigned int>& vecTrafos) { m_vecInvTrafos = vecTrafos; }
	void SetPrimTrafos(std::vector<unsigned int>&& vecTrafos) { m_vecPrimTrafos = std::move(vecTrafos); }
	void SetPrimTrafos(const std::vector<unsigned int>& vecTrafos) { m_vecPrimTrafos = vecTrafos; }
	void SetCenterTrafos(std::vector<unsigned int>&& vecTrafos)
	{
		m_vecCentringTrafos.clear();
		m_vecCenterTrafos = std::move(vecTrafos);
		for(unsigned int iIdx : m_vecCenterTrafos)
			m_vecCentringTrafos.push_back(m_vecTrafos[iIdx]);
	}
	void SetCenterTrafos(const std::vector<unsigned int>& vecTrafos)
	{
		m_vecCentringTrafos.clear();
		m_vecCenterTrafos = vecTrafos;
		for(unsigned int iIdx : m_vecCenterTrafos)
			m_vecCentringTrafos.push_back(m_vecTrafos[iIdx]);
	}
	void SetTransTrafos(std::vector<unsigned int>&& vecTrafos) { m_vecTrans = std::move(vecTrafos); }
	void SetTransTrafos(const std::vector<unsigned int>& vecTrafos) { m_vecTrans = vecTrafos; }


	const std::vector<unsigned int>& GetInvTrafos() const { return m_vecInvTrafos; }
	const std::vector<unsigned int>& GetPrimTrafos() const { return m_vecPrimTrafos; }
	const std::vector<unsigned int>& GetCenterTrafos() const { return m_vecCenterTrafos; }
	const std::vector<unsigned int>& GetTransTrafos() const { return m_vecTrans; }
};


template<class t_real = double>
class SpaceGroups
{
	public:
		typedef std::vector<const SpaceGroup<t_real>*> t_vecSpaceGroups;
		typedef std::map<std::string, SpaceGroup<t_real>> t_mapSpaceGroups;

	private:
		static std::shared_ptr<SpaceGroups<t_real>> s_inst;
		static std::mutex s_mutex;

#ifdef _SGR_NO_SINGLETON
	public:
#endif
		SpaceGroups(const std::string& strFile, const std::string& strXmlRoot="");

	protected:
		t_mapSpaceGroups g_mapSpaceGroups;
		t_vecSpaceGroups g_vecSpaceGroups;
		std::string s_strSrc, s_strUrl;
		bool m_bOk = 0;

	protected:
		bool LoadSpaceGroups(const std::string& strFile, bool bMandatory=1, const std::string& strXmlRoot="");

	public:
		~SpaceGroups();
		static std::shared_ptr<const SpaceGroups<t_real>> GetInstance(const char *pcFile = nullptr);

		const t_mapSpaceGroups* get_space_groups() const;
		const t_vecSpaceGroups* get_space_groups_vec() const;
		const std::string& get_sgsource(bool bUrl=0) const;

		const SpaceGroup<t_real>* Find(const std::string& pcSG) const;

		bool IsOk() const { return m_bOk; }
};

}
#endif
