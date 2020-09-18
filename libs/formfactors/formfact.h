/**
 * Form factors, elements, and scattering length tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2
 */

#ifndef __XTL_FFACT_H__
#define __XTL_FFACT_H__

#include <string>
#include <vector>
#include <complex>
#include <mutex>
#include <boost/optional.hpp>

#include "tlibs/helper/array.h"
#include "tlibs/phys/atoms.h"
#include "tlibs/phys/mag.h"


namespace xtl {


template<typename T=double>
class PeriodicElement
{
	template<typename> friend class PeriodicSystem;

	public:
		typedef T value_type;

	protected:
		std::string strAtom;
		int iNr = -1, iPeriod = -1, iGroup = -1;

		T dMass = T(-1);
		T dRadCov = T(-1), dRadVdW = T(-1);
		T dEIon = T(-1), dEAffin = T(-1);
		T dTMelt = T(-1), dTBoil = T(-1);
		std::string strOrbitals, strBlock;

	public:
		const std::string& GetAtomIdent() const { return strAtom; }

		int GetNr() const { return iNr; }
		int GetPeriod() const { return iPeriod; }
		int GetGroup() const { return iGroup; }

		T GetMass() const { return dMass; }
		T GetRadiusCov() const { return dRadCov; }
		T GetRadiusVdW() const { return dRadVdW; }
		T GetEIon() const { return dEIon; }
		T GetEAffin() const { return dEAffin; }
		T GetTMelt() const { return dTMelt; }
		T GetTBoil() const { return dTBoil; }

		const std::string& GetOrbitals() const { return strOrbitals; }
		const std::string& GetBlock() const { return strBlock; }
};

template<typename T/*=double*/>
class PeriodicSystem
{
	public:
		typedef PeriodicElement<T> elem_type;
		typedef typename elem_type::value_type value_type;

	private:
		static std::shared_ptr<PeriodicSystem> s_inst;
		static std::mutex s_mutex;

#ifdef _FF_NO_SINGLETON
	public:
#endif
		PeriodicSystem(const std::string& strFile, const std::string& strXmlRoot="");

	protected:
		std::vector<elem_type> s_vecAtoms;
		std::string s_strSrc, s_strSrcUrl;

	public:
		~PeriodicSystem();
		static std::shared_ptr<const PeriodicSystem> GetInstance(const char* pcFile=nullptr);

		std::size_t GetNumAtoms() const { return s_vecAtoms.size(); }
		const elem_type& GetAtom(std::size_t i) const
		{ return s_vecAtoms[i]; }

		const elem_type* Find(const std::string& strElem) const;

		const std::string& GetSource() const { return s_strSrc; }
		const std::string& GetSourceUrl() const { return s_strSrcUrl; }
};


// ----------------------------------------------------------------------------



template<typename T=double>
class Formfact
{
	template<typename> friend class FormfactList;

	public:
		typedef T value_type;

	protected:
		std::string strAtom;

		std::vector<T> a;
		std::vector<T> b;
		T c;

	public:
		const std::string& GetAtomIdent() const { return strAtom; }

		T GetFormfact(T G) const
		{
			return tl::formfact<T, std::vector>(G, a, b, c);
		}
};

template<typename T/*=double*/>
class FormfactList
{
	public:
		typedef Formfact<T> elem_type;
		typedef typename elem_type::value_type value_type;

	private:
		static std::shared_ptr<FormfactList> s_inst;
		static std::mutex s_mutex;

#ifdef _FF_NO_SINGLETON
	public:
#endif
		FormfactList(const std::string& strFile, const std::string& strXmlRoot="");

	protected:
		std::vector<elem_type> s_vecAtoms, s_vecIons;
		std::string s_strSrc, s_strSrcUrl;

	public:
		~FormfactList();
		static std::shared_ptr<const FormfactList> GetInstance(const char* pcFile=nullptr);

		std::size_t GetNumAtoms() const { return s_vecAtoms.size(); }
		const elem_type& GetAtom(std::size_t iFormfact) const
		{ return s_vecAtoms[iFormfact]; }

		std::size_t GetNumIons() const { return s_vecIons.size(); }
		const elem_type& GetIon(std::size_t iFormfact) const
		{ return s_vecIons[iFormfact]; }

		const elem_type* Find(const std::string& strElem) const;

		const std::string& GetSource() const { return s_strSrc; }
		const std::string& GetSourceUrl() const { return s_strSrcUrl; }
};


// ----------------------------------------------------------------------------


template<typename T=double>
class MagFormfact
{
	template<typename> friend class MagFormfactList;

	public:
		typedef T value_type;

	protected:
		std::string strAtom;
		std::vector<T> A0, a0;
		std::vector<T> A2, a2;
		std::vector<T> A4, a4;

	public:
		const std::string& GetAtomIdent() const { return strAtom; }

		T GetFormfact(T Q, T g=2) const
		{
			T F;
			std::tie(std::ignore,std::ignore,F) =
				tl::mag_formfact_d<T, std::vector>
					(Q, g, A0,a0, A2,a2);
			return F;
		}

		T GetFormfact(T Q, T L, T S, T J) const
		{
			T F;
			std::tie(std::ignore,std::ignore,F) =
				tl::mag_formfact_f<T, std::vector>
					(Q, L,S,J, A0,a0, A2,a2);
			return F;
		}
};


template<typename T/*=double*/>
class MagFormfactList
{
	public:
		typedef MagFormfact<T> elem_type;
		typedef typename elem_type::value_type value_type;

	private:
		static std::shared_ptr<MagFormfactList> s_inst;
		static std::mutex s_mutex;

#ifdef _FF_NO_SINGLETON
	public:
#endif
		MagFormfactList(const std::string& strFile, const std::string& strXmlRoot="");

	protected:
		std::vector<elem_type> s_vecAtoms;
		std::string s_strSrc, s_strSrcUrl;

	public:
		~MagFormfactList();
		static std::shared_ptr<const MagFormfactList> GetInstance(const char* pcFile=nullptr);

		std::size_t GetNumAtoms() const { return s_vecAtoms.size(); }
		const elem_type& GetAtom(std::size_t iFormfact) const
		{ return s_vecAtoms[iFormfact]; }

		const elem_type* Find(const std::string& strElem) const;

		const std::string& GetSource() const { return s_strSrc; }
		const std::string& GetSourceUrl() const { return s_strSrcUrl; }
};


// ----------------------------------------------------------------------------


template<typename T = std::complex<double>>
class Scatlen
{
	template<typename> friend class ScatlenList;

	public:
		typedef T value_type;	// complex type
		typedef typename value_type::value_type real_type;

	protected:
		std::string strAtom;
		value_type coh;
		value_type incoh;

		value_type xsec_coh;
		value_type xsec_incoh;
		value_type xsec_scat;
		value_type xsec_abs;

		std::vector<const Scatlen<value_type>*> m_vecIsotopes;	// for mixtures: vector of isotopes
		boost::optional<real_type> abund;	// natural abundance of isotope
		boost::optional<real_type> hl;		// in a

	public:
		const std::string& GetAtomIdent() const { return strAtom; }

		const value_type& GetCoherent() const { return coh; }
		const value_type& GetIncoherent() const { return incoh; }

		const value_type& GetXSecCoherent() const { return xsec_coh; }
		const value_type& GetXSecIncoherent() const { return xsec_incoh; }
		const value_type& GetXSecScatter() const { return xsec_scat; }
		const value_type& GetXSecAbsorption() const { return xsec_abs; }

		const boost::optional<real_type>& GetAbundance() const { return abund; }
		const boost::optional<real_type>& GetHalflife() const { return hl; }

		const std::vector<const Scatlen<value_type>*> GetIsotopes() const { return m_vecIsotopes; }
};


template<typename T/*=double*/>
class ScatlenList
{
	public:
		typedef Scatlen<std::complex<T>> elem_type;
		typedef typename elem_type::value_type value_type;
		typedef typename elem_type::real_type real_type;

	private:
		static std::shared_ptr<ScatlenList> s_inst;
		static std::mutex s_mutex;

#ifdef _FF_NO_SINGLETON
	public:
#endif
		ScatlenList(const std::string& strFile, const std::string& strXmlRoot="");

	protected:
		std::vector<elem_type> s_vecElems, s_vecIsotopes;
		std::string s_strSrc, s_strSrcUrl;

	public:
		~ScatlenList();
		static std::shared_ptr<const ScatlenList> GetInstance(const char* pcFile=nullptr);

		std::size_t GetNumElems() const { return s_vecElems.size(); }
		const elem_type& GetElem(std::size_t i) const
		{ return s_vecElems[i]; }

		std::size_t GetNumIsotopes() const { return s_vecIsotopes.size(); }
		const elem_type& GetIsotope(std::size_t i) const
		{ return s_vecIsotopes[i]; }

		const elem_type* Find(const std::string& strElem) const;

		const std::string& GetSource() const { return s_strSrc; }
		const std::string& GetSourceUrl() const { return s_strSrcUrl; }
};


}
#endif
