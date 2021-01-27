/**
 * interface for S(q,w) models
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015, 2016
 * @license GPLv2
 */

#include "sqwbase.h"
#include "tlibs/log/log.h"


/**
 * if the variable "strKey" is known, update it with the value "strNewVal"
 */
bool SqwBase::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	std::vector<t_var> vecVars = GetVars();
	for(const t_var& var : vecVars)
	{
		if(strKey == std::get<0>(var))
		{
			t_var varNew = var;
			std::get<2>(varNew) = strNewVal;
			SetVars(std::vector<t_var>({varNew}));
			return true;
		}
	}

	return false;
}


/**
 * if the variable "strKey" is known, update its error with the value "strNewErr"
 */
bool SqwBase::SetErrIfAvail(const std::string& strKey, const std::string& strNewErr)
{
	for(t_var_fit& var : m_vecFit)
	{
		if(strKey == std::get<0>(var))
		{
			std::get<1>(var) = strNewErr;
			return true;
		}
	}

	return false;
}


/**
 * if the variable "strKey" is known, update its range with the value "strNewRange"
 */
bool SqwBase::SetRangeIfAvail(const std::string& strKey, const std::string& strNewRange)
{
	for(t_var_fit& var : m_vecFit)
	{
		if(strKey == std::get<0>(var))
		{
			std::get<3>(var) = strNewRange;
			return true;
		}
	}

	return false;
}


/**
 * replaces all fit variables
 */
void SqwBase::InitFitVars(const std::vector<t_var_fit>& vecFit)
{
	m_vecFit = vecFit;
}


/**
 * keeps the current fit variable vector and only replaces changed values
 */
void SqwBase::SetFitVars(const std::vector<t_var_fit>& vecFitVars)
{
	for(const t_var_fit& var : vecFitVars)
	{
		const std::string& name = std::get<0>(var);
		auto iter = std::find_if(m_vecFit.begin(), m_vecFit.end(), [&name](const t_var_fit& fitvar)
		{
			return std::get<0>(fitvar) == name;
		});

		if(iter != m_vecFit.end())
		{
			std::get<1>(*iter) = std::get<1>(var);
			std::get<2>(*iter) = std::get<2>(var);
			std::get<3>(*iter) = std::get<3>(var);
		}
		else
		{
			tl::log_err("Tried to update non-existing fit variable \"", name, "\".");
		}
	}
}


const SqwBase& SqwBase::operator=(const SqwBase& sqw)
{
	this->m_bOk = sqw.m_bOk;
	this->m_vecFit = sqw.m_vecFit;

	return *this;
}
