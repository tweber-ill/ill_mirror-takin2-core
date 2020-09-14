/**
 * monte carlo convolution tool -> convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015, 2016
 * @license GPLv2
 */

#include "ConvoDlg.h"

#include <QMutex>
#include <QMessageBox>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <string>
#include <exception>

using t_real = t_real_reso;
using t_real_min = double;


namespace minuit = ROOT::Minuit2;


class StopRequestedEx : public std::runtime_error
{
public:
	StopRequestedEx(const char* msg) : runtime_error{msg} {}
};



/**
  * interface between ConvoDlg and Minuit
  */
class MinuitFunc : public minuit::FCNBase
{
public:
	MinuitFunc(ConvoDlg* convodlg, const ConvoDlg::t_sqwparams* sqwparams)
		: m_convodlg{convodlg}, m_sqwparams{sqwparams}, m_seed{tl::get_rand_seed()}
	{}

	MinuitFunc(const MinuitFunc&) = delete;
	const MinuitFunc& operator=(const MinuitFunc&) = delete;

	virtual ~MinuitFunc() = default;

	
	virtual t_real_min operator()(const std::vector<t_real_min>& params) const override
	{
		std::lock_guard<QMutex> _lock{m_mtxMinuit};

		// set model parameters, [ident, val, err]
		std::vector<std::tuple<std::string, std::string, std::string>> sqwparams;
		for(std::size_t paramidx=0; paramidx<params.size(); ++paramidx)
		{
			const std::string& name = std::get<0>((*m_sqwparams)[paramidx]);
			sqwparams.emplace_back(std::make_tuple(name, tl::var_to_str(params[paramidx]), ""));
		}
		m_convodlg->SetSqwParams(sqwparams);

		// start convolution simulator for new parameters and wait for it to finish
		m_convodlg->StartSim1D(true, m_seed);

		// if a stop is requested, we have no other way of getting out of here than throwing an exception...
		// if not in StartFit(), this exception will be handled at the latest by TakAppl::notify()
		if(m_convodlg->StopRequested())
			throw StopRequestedEx{"Convolution fitter stop requested."};

		return m_convodlg->GetChi2();
	}


	virtual t_real_min Up() const override
	{
		// sigma^2
		return 1.;
	}


private:
	ConvoDlg* m_convodlg = nullptr;
	const ConvoDlg::t_sqwparams* m_sqwparams = nullptr;
	unsigned int m_seed{1234};

	mutable QMutex m_mtxMinuit{};
};



/**
 * start 1d or 2d convolution fits
 */
void ConvoDlg::StartFit()
{
	if(check2dMap->isChecked())
	{
		QMessageBox::critical(this, "Error", "2D fitting is not yet implemented.");
		return;
	}


	// [ ident, type, value, error, fit? ]
	t_sqwparams sqwparams = GetSqwParams(true);
	if(sqwparams.size() == 0)
	{
		QMessageBox::critical(this, "Error", "No fit parameters defined.");
		return;
	}


	// stop any previous fits
	Stop();
	m_atStop.store(false);


	minuit::MnUserParameters params;
	for(const auto& sqwparam : sqwparams)
	{
		t_real val = tl::str_to_var<t_real_min>(std::get<2>(sqwparam));
		t_real err = tl::str_to_var<t_real_min>(std::get<3>(sqwparam));
		params.Add(std::get<0>(sqwparam), val, err);
	}


	// minimise
	bool mini_valid = false;
	std::unique_ptr<minuit::FunctionMinimum> mini;

	try
	{

		MinuitFunc fkt{this, &sqwparams};

		minuit::MnStrategy strat(spinStrategy->value());
		std::unique_ptr<minuit::MnApplication> minimiser;
		if(comboFitter->currentIndex() == 0)
		{
			minimiser.reset(new minuit::MnSimplex(fkt, params, strat));
		}
		else if(comboFitter->currentIndex() == 1)
		{
			minimiser.reset(new minuit::MnMigrad(fkt, params, strat));
		}
		else
		{
			QMessageBox::critical(this, "Error", "Invalid minimiser.");
			return;
		}

		mini.reset(new minuit::FunctionMinimum((*minimiser)(spinMaxCalls->value(), spinTolerance->value())));
		mini_valid = mini->IsValid() && mini->HasValidParameters() && mini->UserState().IsValid();
	}
	catch(const StopRequestedEx& req)
	{
		tl::log_info(req.what());
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	if(!mini_valid || !mini)
		QMessageBox::critical(this, "Error", "Convolution fitter did not converge.");


	// get back minimised parameters, [ident, val, err]
	if(mini)
	{
		tl::log_debug("Final fit results:\n", *mini);

		std::vector<std::tuple<std::string, std::string, std::string>> newsqwparams;
		for(std::size_t paramidx=0; paramidx<sqwparams.size(); ++paramidx)
		{
			const std::string& name = std::get<0>(sqwparams[paramidx]);
			std::string val = tl::var_to_str(mini->UserState().Value(name));
			std::string err = tl::var_to_str(mini->UserState().Error(name));

			newsqwparams.emplace_back(std::make_tuple(name, val, err));
		}
		this->SetSqwParams(newsqwparams);
	}
}
