/**
 * Creates needed tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2
 */

#include <iostream>
#include <sstream>
#include <set>
#include <limits>
#include <unordered_map>

#include "tlibs/string/string.h"
#include "tlibs/file/prop.h"
#include "tlibs/log/log.h"
#include "tlibs/math/linalg.h"
#include "libs/spacegroups/sghelper.h"
#ifndef NO_CLP
	#include "libs/spacegroups/spacegroup_clp.h"
#endif

#ifndef USE_BOOST_REX
	#include <regex>
	namespace rex = ::std;
#else
	#include <boost/tr1/regex.hpp>
	namespace rex = ::boost;
#endif

#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>

namespace prop = boost::property_tree;
namespace algo = boost::algorithm;
namespace ublas = boost::numeric::ublas;

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_mat = ublas::matrix<t_real>;


unsigned g_iPrec = std::numeric_limits<t_real>::max_digits10-1;

// ============================================================================

#include "gentab_web.cpp"
#ifndef NO_CLP
	#include "gentab_clp.cpp"
#endif

// ============================================================================


static std::unordered_map<std::string, int> g_mapElems;

/**
 * periodic table of elements
 */
bool gen_elements()
{
	g_mapElems.clear();

	using t_propval = tl::Prop<std::string>::t_propval;

	tl::Prop<std::string> propIn, propOut;
	propIn.SetSeparator('/');
	propOut.SetSeparator('.');

	if(!propIn.Load("tmp/elements.xml", tl::PropType::XML))
	{
		tl::log_err("Cannot load periodic table of elements \"tmp/elements.xml\".");
		return false;
	}


	// iterate over all elements
	std::vector<t_propval> vecElems = propIn.GetFullChildNodes("/list");
	std::size_t iElem = 0;
	for(const t_propval& elem : vecElems)
	{
		if(elem.first != "atom") continue;

		try
		{
			tl::Prop<std::string> propelem(elem.second, '/');

			std::string strName = propelem.Query<std::string>("<xmlattr>/id", "");
			if(strName == "" || strName == "Xx") continue;

			t_real dMass = t_real(-1);
			t_real dRadCov=t_real(-1), dRadVdW=t_real(-1);
			t_real dEIon=t_real(-1), dEAffin(-1);
			t_real dTMelt=t_real(-1), dTBoil=t_real(-1);
			int iNr=-1, iPeriod=-1, iGroup=-1;
			std::string strConfig, strBlock;

			// iterate over all properties
			for(auto iterVal=propelem.GetProp().begin(); iterVal!=propelem.GetProp().end(); ++iterVal)
			{
				tl::Prop<std::string> propVal(iterVal->second, '/');
				std::string strKey = propVal.Query<std::string>("<xmlattr>/dictRef", "");
				std::string strVal = propVal.Query<std::string>("/", "");
				//std::cout << strKey << " = " << strVal << std::endl;

				if(strKey.find("atomicNumber") != std::string::npos)
					iNr = tl::str_to_var<int>(strVal);
				else if(strKey.find("electronicConfiguration") != std::string::npos)
					strConfig = strVal;
				else if(strKey.find("periodTableBlock") != std::string::npos)
					strBlock = strVal;
				else if(strKey.find("period") != std::string::npos)
					iPeriod = tl::str_to_var<int>(strVal);
				else if(strKey.find("group") != std::string::npos)
					iGroup = tl::str_to_var<int>(strVal);
				else if(strKey.find("exactMass") != std::string::npos)
					dMass = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("radiusCovalent") != std::string::npos)
					dRadCov = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("radiusVDW") != std::string::npos)
					dRadVdW = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("ionization") != std::string::npos)
					dEIon = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("electronAffinity") != std::string::npos)
					dEAffin = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("melting") != std::string::npos)
					dTMelt = tl::str_to_var<t_real>(strVal);
				else if(strKey.find("boiling") != std::string::npos)
					dTBoil = tl::str_to_var<t_real>(strVal);
			}

			std::ostringstream ostr;
			ostr << "pte.elem_" << iElem;
			std::string strElem = ostr.str();

			propOut.Add(strElem + ".name", strName);
			propOut.Add(strElem + ".num", tl::var_to_str(iNr, g_iPrec));
			propOut.Add(strElem + ".period", tl::var_to_str(iPeriod, g_iPrec));
			propOut.Add(strElem + ".group", tl::var_to_str(iGroup, g_iPrec));
			propOut.Add(strElem + ".orbitals", strConfig);
			propOut.Add(strElem + ".block", strBlock);
			propOut.Add(strElem + ".m", tl::var_to_str(dMass, g_iPrec));
			propOut.Add(strElem + ".r_cov", tl::var_to_str(dRadCov, g_iPrec));
			propOut.Add(strElem + ".r_vdW", tl::var_to_str(dRadVdW, g_iPrec));
			propOut.Add(strElem + ".E_ion", tl::var_to_str(dEIon, g_iPrec));
			propOut.Add(strElem + ".E_affin", tl::var_to_str(dEAffin, g_iPrec));
			propOut.Add(strElem + ".T_melt", tl::var_to_str(dTMelt, g_iPrec));
			propOut.Add(strElem + ".T_boil", tl::var_to_str(dTBoil, g_iPrec));

			g_mapElems[strName] = iNr;
		}
		catch(const std::exception& ex)
		{
			tl::log_err("Element ", iElem, ": ", ex.what());
		}

		++iElem;
	}


	propOut.Add("pte.num_elems", iElem);

	propOut.Add("pte.source", "Periodic table of the elements obtained from the "
		"<a href=\"http://dx.doi.org/10.1021/ci050400b\">Blue Obelisk Data Repository</a>.");
	propOut.Add("pte.source_url", "https://github.com/egonw/bodr/blob/master/bodr/elements/elements.xml");

	if(!propOut.Save("res/data/elements.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/elements.xml.gz\".");
		return false;
	}

	return true;
}



// ============================================================================


/**
 * space groups (alternative)
 */
bool gen_spacegroups()
{
	using t_propval = tl::Prop<std::string>::t_propval;

	tl::Prop<std::string> propIn, propOut;
	propIn.SetSeparator('/');
	propOut.SetSeparator('.');

	if(!propIn.Load("tmp/space-groups.xml", tl::PropType::XML))
	{
		tl::log_err("Cannot load space group table \"tmp/space-groups.xml\".");
		return false;
	}


	// iterate over all groups
	std::vector<t_propval> vecGroups = propIn.GetFullChildNodes("/list");
	std::size_t iGroup = 0;
	std::set<std::string> setNames;
	for(const t_propval& grp : vecGroups)
	{
		if(grp.first != "group") continue;

		try
		{
			tl::Prop<std::string> propgrp(grp.second, '/');

			std::string strId = propgrp.Query<std::string>("<xmlattr>/id", "");
			if(strId == "") continue;
			int iNr = tl::str_to_var<int>(strId);

			std::string strName = propgrp.Query<std::string>("<xmlattr>/HM", "");
			if(strName == "") continue;
			tl::find_all_and_replace<std::string>(strName, ":1", "");
			tl::find_all_and_replace<std::string>(strName, ":2", "");
			tl::find_all_and_replace<std::string>(strName, ":3", "");
			xtl::convert_hm_symbol(strName);

			// find an unique name
			std::size_t iNameCtr = 1;
			while(1)
			{
				std::ostringstream ostrNewName;
				ostrNewName << strName;
				if(iNameCtr >= 2)
					ostrNewName << " [" << iNameCtr << "]";
				if(setNames.find(ostrNewName.str()) == setNames.end())
				{
					strName = ostrNewName.str();
					setNames.insert(strName);
					break;
				}

				++iNameCtr;
			}

			std::ostringstream ostr;
			ostr << "sgroups.group_" << iGroup;
			std::string strGroup = ostr.str();

			std::vector<t_propval> vecTrafos = propgrp.GetFullChildNodes("/");
			std::size_t iTrafo = 0;
			for(const t_propval& trafo : vecTrafos)
			{
				if(trafo.first != "transform") continue;

				tl::Prop<std::string> proptrafo(trafo.second, '/');
				std::string strTrafo = proptrafo.Query<std::string>("/", "");

				std::ostringstream ostrTrafo;
				ostrTrafo << strGroup << ".trafo_" << iTrafo;
				std::string strTrafoKey = ostrTrafo.str();

				t_mat matTrafo = xtl::get_desc_trafo<t_real>(strTrafo);
				std::string strAttribs = "; ";
				if(tl::is_identity_matrix(matTrafo) ||
					tl::is_inverting_matrix<t_mat>(tl::submatrix(matTrafo,3,3)))
					strAttribs += "i";
				if(tl::is_identity_matrix(matTrafo) ||
					tl::is_centering_matrix<t_mat>(matTrafo))
					strAttribs += "c";
				/*if(tl::is_identity_matrix(matTrafo) ||
					(tl::has_translation_components<t_mat>(matTrafo) &&
						!tl::is_centering_matrix<t_mat>(matTrafo) &&
						!tl::is_inverting_matrix<t_mat>(tl::submatrix(matTrafo,3,3))))
					strAttribs += "s";*/
				if(tl::is_identity_matrix(matTrafo) ||
					tl::has_translation_components<t_mat>(matTrafo))
					strAttribs += "t";

				propOut.Add(strTrafoKey, tl::var_to_str(matTrafo, g_iPrec) + strAttribs);
				++iTrafo;
			}

			propOut.Add(strGroup + ".number", tl::var_to_str(iNr, g_iPrec));
			propOut.Add(strGroup + ".name", strName);
			propOut.Add(strGroup + ".num_trafos", tl::var_to_str(iTrafo, g_iPrec));
		}
		catch(const std::exception& ex)
		{
			tl::log_err("Space group ", iGroup, ": ", ex.what());
		}

		++iGroup;
	}


	propOut.Add("sgroups.num_groups", iGroup);

	propOut.Add("sgroups.source", "Space groups obtained from the "
		"<a href=\"http://dx.doi.org/10.1021/ci050400b\">Blue Obelisk Data Repository</a>.");
	propOut.Add("sgroups.source_url", "https://github.com/egonw/bodr/blob/master/bodr/crystal/space-groups.xml");

	if(!propOut.Save("res/data/sgroups.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/sgroups.xml.gz\".");
		return false;
	}

	return true;
}



// ============================================================================


struct ffact
{
	std::string strName;

	t_cplx cCohb, cIncb;

	t_real dScatXs, dAbsXs;
	t_real dCohXs, dIncXs;

	bool bIncCohXs = 0;
	bool bIncIncXs = 0;
	bool bIncScatXs = 0;

	std::string strAbund;
	t_real dAbOrHL;
	bool bAb;
};

bool gen_scatlens_npy()
{
	tl::Prop<std::string> propIn, propOut;
	propIn.SetSeparator('/');
	propOut.SetSeparator('.');

	if(!propIn.Load("tmp/scattering_lengths.json", tl::PropType::JSON))
	{
		tl::log_err("Cannot load scattering length table \"tmp/scattering_lengths.json\".");
		return false;
	}

	std::vector<std::string> vecNuclei = propIn.GetChildNodes("/");
	std::vector<ffact> vecFfacts;

	for(const std::string& strNucl : vecNuclei)
	{
		ffact ff;

		ff.strName = strNucl;
		ff.cCohb = propIn.Query<t_real>("/" + strNucl + "/Coh b");
		ff.cIncb = propIn.Query<t_real>("/" + strNucl + "/Inc b");
		ff.dAbsXs = propIn.Query<t_real>("/" + strNucl + "/Abs xs") * t_real(100); // in fm^2
		ff.dCohXs = propIn.Query<t_real>("/" + strNucl + "/Coh xs") * t_real(100);
		ff.dIncXs = propIn.Query<t_real>("/" + strNucl + "/Inc xs") * t_real(100);
		ff.dScatXs = propIn.Query<t_real>("/" + strNucl + "/Scatt xs");
		ff.strAbund = propIn.Query<std::string>("/" + strNucl + "/conc");

		// include total scattering cross-section only if it is not the sum
		// of the coherent and incoherent cross-sections
		if(!tl::float_equal(ff.dScatXs, ff.dCohXs + ff.dIncXs, 1.))
		{
			ff.bIncScatXs = 1;
			//tl::log_warn("Mismatch in scattering cross-section for ", ff.strName, ": ",
			//	ff.dCohXs + ff.dIncXs, " != ", ff.dScatXs, ".");
		}

		if(!tl::float_equal(ff.dCohXs, (ff.cCohb*std::conj(ff.cCohb)).real()*t_real(4)*tl::get_pi<t_real>(), 1.))
		{
			ff.bIncCohXs = 1;
			//tl::log_warn("Mismatch in coherent cross-section for ", ff.strName, ": ",
			//	(ff.cCohb*std::conj(ff.cCohb)).real()*t_real(4)*tl::get_pi<t_real>(),
			//	" != ", ff.dCohXs, ".");
		}

		if(!tl::float_equal(ff.dIncXs, (ff.cIncb*std::conj(ff.cIncb)).real()*t_real(4)*tl::get_pi<t_real>(), 1.))
		{
			ff.bIncIncXs = 1;
			//tl::log_warn("Mismatch in incoherent cross-section for ", ff.strName, ": ",
			//	(ff.cIncb*std::conj(ff.cIncb)).real()*t_real(4)*tl::get_pi<t_real>(),
			//	" != ", ff.dIncXs, ".");
		}


		// complex?
		auto vecValsCohb = propIn.GetChildValues<t_real>("/" + strNucl + "/Coh b");
		auto vecValsIncb = propIn.GetChildValues<t_real>("/" + strNucl + "/Inc b");

		if(vecValsCohb.size() >= 2)
		{
			ff.cCohb.real(vecValsCohb[0]);
			ff.cCohb.imag(vecValsCohb[1]);
		}
		if(vecValsIncb.size() >= 2)
		{
			ff.cIncb.real(vecValsIncb[0]);
			ff.cIncb.imag(vecValsIncb[1]);
		}

		ff.dAbOrHL = t_real(0);
		ff.bAb = get_abundance_or_hl(ff.strAbund, ff.dAbOrHL);

		vecFfacts.emplace_back(std::move(ff));
	}


	// sort elements if elements map is not empty
	if(g_mapElems.size())
	{
		std::stable_sort(vecFfacts.begin(), vecFfacts.end(),
			[](const ffact& ff1, const ffact& ff2) -> bool
			{
				std::string strName1 = tl::remove_chars(ff1.strName, std::string("+-0123456789"));
				std::string strName2 = tl::remove_chars(ff2.strName, std::string("+-0123456789"));

				auto iter1 = g_mapElems.find(strName1);
				auto iter2 = g_mapElems.find(strName2);

				if(iter1 == g_mapElems.end())
				{
					tl::log_err("Element ", strName1, " not in table!");
					return 0;
				}
				if(iter2 == g_mapElems.end())
				{
					tl::log_err("Element ", strName2, " not in table!");
					return 0;
				}

				return iter1->second < iter2->second;
			});
	}


	// write database
	std::size_t iNucl = 0;
	for(const ffact& ff : vecFfacts)
	{
		std::ostringstream ostr;
		ostr << "scatlens.atom_" << iNucl;
		std::string strAtom = ostr.str();

		propOut.Add(strAtom + ".name", ff.strName);

		// scattering lengths
		propOut.Add(strAtom + ".coh", tl::var_to_str(ff.cCohb, g_iPrec));
		propOut.Add(strAtom + ".incoh", tl::var_to_str(ff.cIncb, g_iPrec));

		// cross-sections
		if(ff.bIncCohXs)
			propOut.Add(strAtom + ".xsec_coh", tl::var_to_str(ff.dCohXs, g_iPrec));
		if(ff.bIncIncXs)
			propOut.Add(strAtom + ".xsec_incoh", tl::var_to_str(ff.dIncXs, g_iPrec));
		if(ff.bIncScatXs)
			propOut.Add(strAtom + ".xsec_scat", tl::var_to_str(ff.dScatXs, g_iPrec));
		propOut.Add(strAtom + ".xsec_absorp", tl::var_to_str(ff.dAbsXs, g_iPrec));

		// abundances
		if(ff.bAb)
			propOut.Add(strAtom + ".abund", tl::var_to_str(ff.dAbOrHL, g_iPrec));
		else
			propOut.Add(strAtom + ".hl", tl::var_to_str(ff.dAbOrHL, g_iPrec));

		++iNucl;
	}

	propOut.Add("scatlens.num_atoms", tl::var_to_str(vecNuclei.size()));

	propOut.Add("scatlens.source", "Scattering lengths and cross-sections extracted from NeutronPy (by D. Fobes)"
		" (which itself is based on <a href=\"http://dx.doi.org/10.1080/10448639208218770\">this paper</a>).");
	propOut.Add("scatlens.source_url", "https://github.com/neutronpy/neutronpy/blob/master/neutronpy/database/scattering_lengths.json");

	if(!propOut.Save("res/data/scatlens.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/scatlens.xml.gz\".");
		return false;
	}

	return true;
}



// ============================================================================


struct mag_ffact
{
	std::string strName;
	std::string strJ0_A, strJ0_a;
	std::string strJ2_A, strJ2_a;
	std::string strJ4_A, strJ4_a;
};

bool gen_magformfacts_npy()
{
	tl::Prop<std::string> propIn, propOut;
	propIn.SetSeparator('/');
	propOut.SetSeparator('.');

	if(!propIn.Load("tmp/magnetic_form_factors.json", tl::PropType::JSON))
	{
		tl::log_err("Cannot load scattering length table \"tmp/magnetic_form_factors.json\".");
		return false;
	}

	std::vector<std::string> vecNuclei = propIn.GetChildNodes("/");
	std::vector<mag_ffact> vecFfacts;

	for(const std::string& strNucl : vecNuclei)
	{
		auto vecJ0 = propIn.GetChildValues<t_real>("/" + strNucl + "/j0");
		auto vecJ2 = propIn.GetChildValues<t_real>("/" + strNucl + "/j2");
		auto vecJ4 = propIn.GetChildValues<t_real>("/" + strNucl + "/j4");

		mag_ffact ffact;
		ffact.strName = strNucl;
		std::string& strJ0A = ffact.strJ0_A;
		std::string& strJ2A = ffact.strJ2_A;
		std::string& strJ4A = ffact.strJ4_A;
		std::string& strJ0a = ffact.strJ0_a;
		std::string& strJ2a = ffact.strJ2_a;
		std::string& strJ4a = ffact.strJ4_a;

		for(std::size_t iJ=0; iJ<vecJ0.size(); ++iJ)
		{
			t_real dVal = vecJ0[iJ];
			bool bEven = tl::is_even(iJ);
			if(bEven)
			{
				if(strJ0A != "") strJ0A += "; ";
				strJ0A += tl::var_to_str(dVal, g_iPrec);
			}
			else
			{
				if(strJ0a != "") strJ0a += "; ";
				strJ0a += tl::var_to_str(dVal, g_iPrec);
			}
		}
		for(std::size_t iJ=0; iJ<vecJ2.size(); ++iJ)
		{
			t_real dVal = vecJ2[iJ];
			bool bEven = tl::is_even(iJ);
			if(bEven)
			{
				if(strJ2A != "") strJ2A += "; ";
				strJ2A += tl::var_to_str(dVal, g_iPrec);
			}
			else
			{
				if(strJ2a != "") strJ2a += "; ";
				strJ2a += tl::var_to_str(dVal, g_iPrec);
			}
		}
		for(std::size_t iJ=0; iJ<vecJ4.size(); ++iJ)
		{
			t_real dVal = vecJ4[iJ];
			bool bEven = tl::is_even(iJ);
			if(bEven)
			{
				if(strJ4A != "") strJ4A += "; ";
				strJ4A += tl::var_to_str(dVal, g_iPrec);
			}
			else
			{
				if(strJ4a != "") strJ4a += "; ";
				strJ4a += tl::var_to_str(dVal, g_iPrec);
			}
		}

		vecFfacts.emplace_back(std::move(ffact));
	}


	// sort elements
	std::stable_sort(vecFfacts.begin(), vecFfacts.end(),
		[](const mag_ffact& ff1, const mag_ffact& ff2) -> bool
		{
			std::string strName1 = tl::remove_chars(ff1.strName, std::string("+-0123456789"));
			std::string strName2 = tl::remove_chars(ff2.strName, std::string("+-0123456789"));

			auto iter1 = g_mapElems.find(strName1);
			auto iter2 = g_mapElems.find(strName2);

			if(iter1 == g_mapElems.end())
				tl::log_err("Element ", strName1, " not in table!");
			if(iter2 == g_mapElems.end())
				tl::log_err("Element ", strName2, " not in table!");

			return iter1->second < iter2->second;
		});


	// write database
	std::size_t iNucl = 0;
	for(const mag_ffact& ffact : vecFfacts)
	{
		std::ostringstream ostr;
		ostr << "atom_" << iNucl;
		std::string strAtom = ostr.str();

		propOut.Add("magffacts.j0." + strAtom + ".name", ffact.strName);
		propOut.Add("magffacts.j2." + strAtom + ".name", ffact.strName);
		propOut.Add("magffacts.j4." + strAtom + ".name", ffact.strName);

		propOut.Add("magffacts.j0." + strAtom + ".A", ffact.strJ0_A);
		propOut.Add("magffacts.j0." + strAtom + ".a", ffact.strJ0_a);

		propOut.Add("magffacts.j2." + strAtom + ".A", ffact.strJ2_A);
		propOut.Add("magffacts.j2." + strAtom + ".a", ffact.strJ2_a);

		propOut.Add("magffacts.j4." + strAtom + ".A", ffact.strJ4_A);
		propOut.Add("magffacts.j4." + strAtom + ".a", ffact.strJ4_a);

		++iNucl;
	}

	propOut.Add("magffacts.num_atoms", tl::var_to_str(vecNuclei.size()));

	propOut.Add("magffacts.source", "Magnetic form factor coefficients extracted from NeutronPy (by D. Fobes).");
	propOut.Add("magffacts.source_url", "https://github.com/neutronpy/neutronpy/blob/master/neutronpy/database/magnetic_form_factors.json");

	if(!propOut.Save("res/data/magffacts.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/magffacts.xml.gz\".");
		return false;
	}

	return true;
}


// ============================================================================


int main()
{
#ifdef NO_TERM_CMDS
	tl::Log::SetUseTermCmds(0);
#endif

	std::cout << "Generating periodic table of elements ... ";
	bool bHasElems = gen_elements();
	if(bHasElems) std::cout << "OK" << std::endl;

#ifndef NO_CLP
	std::cout << "Generating atomic form factor coefficient table ... ";
	if(gen_formfacts_clp()) std::cout << "OK" << std::endl;
#endif

	std::cout << "Generating scattering length table ... ";
	if(gen_scatlens_npy())
	{
		std::cout << "OK" << std::endl;
	}
	else
	{
		std::cout << "FAILED.\nGenerating scattering length table (alternative) ... ";
		if(gen_scatlens()) std::cout << "OK" << std::endl;
	}

	std::cout << "Generating space group table ... ";
	if(gen_spacegroups())
	{
		std::cout << "OK" << std::endl;
	}
	else
	{
		std::cout << "FAILED.\n";
#ifndef NO_CLP
		std::cout << "Generating space group type table (alternative) ... ";
		if(gen_spacegroups_clp()) std::cout << "OK" << std::endl;
#endif
	}

	//std::cout << "Generating magnetic form factor coefficient table ... ";
	//if(gen_magformfacts()) std::cout << "OK" << std::endl;

	if(bHasElems)
	{
		std::cout << "Generating magnetic form factor coefficient table ... ";
		if(gen_magformfacts_npy()) std::cout << "OK" << std::endl;
	}
	else
	{
		tl::log_err("Cannot create magnetic form factor coefficient table, because required periodic table is invalid.");
	}

	return 0;
}
