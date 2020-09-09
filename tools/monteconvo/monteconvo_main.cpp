/**
 * convolution simulation -- CLI program
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2020
 * @license GPLv2
 */

#include <boost/asio/io_service.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/scope_exit.hpp>
#include <boost/program_options.hpp>

#include "sqwfactory.h"

#include "libs/version.h"
#include "libs/globals.h"

#include "tlibs/file/file.h"
#include "tlibs/log/debug.h"
#include "tlibs/time/stopwatch.h"
#include "tlibs/helper/thread.h"


namespace asio = boost::asio;
namespace sys = boost::system;
namespace opts = boost::program_options;

using t_real = t_real_reso;



// ----------------------------------------------------------------------------
// configuration
struct ConvoConfig
{
	t_real h_from{}, k_from{}, l_from{}, E_from{};
	t_real h_to{}, k_to{}, l_to{}, E_to{};
	t_real h_to_2{}, k_to_2{}, l_to_2{}, E_to_2{};
	t_real kfix{};
	t_real tolerance{};
	t_real S_scale{1}, S_slope{0}, S_offs{0};

	int neutron_count{500};
	int sample_step_count{1};
	int step_count{256};

	bool scan_2d{false};
	bool recycle_neutrons{true};
	bool normalise{true};
	bool flip_coords{false};
	bool has_scanfile{false};

	int algo{1};
	int fixedk{1};
	int mono_foc{1}, ana_foc{1};
	int scanaxis{4}, scanaxis2{0};

	std::string crys{}, instr{};
	std::string sqw{}, sqw_conf{};
	std::string scanfile{};
	std::string counter{}, monitor{};
	std::string temp_override{}, field_override{};
	std::string autosave{};
};


/**
 * loads the configuration for the convolution from a job file
 */
static ConvoConfig load_config(tl::Prop<std::string>& xml)
{
	ConvoConfig cfg;
	const std::string strXmlRoot("taz/");

	// real values
	boost::optional<t_real> odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/h_from"); if(odVal) cfg.h_from = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/k_from"); if(odVal) cfg.k_from = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/l_from"); if(odVal) cfg.l_from = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/E_from"); if(odVal) cfg.E_from = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/h_to"); if(odVal) cfg.h_to = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/k_to"); if(odVal) cfg.k_to = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/l_to"); if(odVal) cfg.l_to = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/E_to"); if(odVal) cfg.E_to = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/h_to_2"); if(odVal) cfg.h_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/k_to_2"); if(odVal) cfg.k_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/l_to_2"); if(odVal) cfg.l_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/E_to_2"); if(odVal) cfg.E_to_2 = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/kfix"); if(odVal) cfg.kfix = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"convofit/tolerance"); if(odVal) cfg.tolerance = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/S_scale"); if(odVal) cfg.S_scale = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/S_slope"); if(odVal) cfg.S_slope = *odVal;
	odVal = xml.QueryOpt<t_real>(strXmlRoot+"monteconvo/S_offs"); if(odVal) cfg.S_offs = *odVal;

	// int values
	boost::optional<int> oiVal;
	oiVal = xml.QueryOpt<int>(strXmlRoot+"monteconvo/neutron_count"); if(oiVal) cfg.neutron_count = *oiVal;
	oiVal = xml.QueryOpt<int>(strXmlRoot+"monteconvo/sample_step_count"); if(oiVal) cfg.sample_step_count = *oiVal;
	oiVal = xml.QueryOpt<int>(strXmlRoot+"monteconvo/step_count"); if(oiVal) cfg.step_count = *oiVal;
	//oiVal = xml.QueryOpt<int>(strXmlRoot+"convofit/strategy"); if(oiVal) cfg.strategy = *oiVal;
	//oiVal = xml.QueryOpt<int>(strXmlRoot+"convofit/max_calls"); if(oiVal) cfg.max_calls = *oiVal;

	// bool values
	boost::optional<int> obVal;
	obVal = xml.QueryOpt<int>(strXmlRoot+"monteconvo/scan_2d"); if(obVal) cfg.scan_2d = *obVal != 0;
	obVal = xml.QueryOpt<int>(strXmlRoot+"convofit/recycle_neutrons"); if(obVal) cfg.recycle_neutrons = *obVal != 0;
	obVal = xml.QueryOpt<int>(strXmlRoot+"convofit/normalise"); if(obVal) cfg.normalise = *obVal != 0;
	obVal = xml.QueryOpt<int>(strXmlRoot+"convofit/flip_coords"); if(obVal) cfg.flip_coords = *obVal != 0;
	obVal = xml.QueryOpt<int>(strXmlRoot+"monteconvo/has_scanfile"); if(obVal) cfg.has_scanfile = *obVal != 0;

	// indices for gui comboboxes
	boost::optional<int> oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"monteconvo/algo"); if(oCmb) cfg.algo = *oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"monteconvo/fixedk"); if(oCmb) cfg.fixedk = *oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"monteconvo/mono_foc"); if(oCmb) cfg.mono_foc = *oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"monteconvo/ana_foc"); if(oCmb) cfg.ana_foc = *oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"convofit/scanaxis"); if(oCmb) cfg.scanaxis = *oCmb;
	oCmb = xml.QueryOpt<int>(strXmlRoot+"convofit/scanaxis2"); if(oCmb) cfg.scanaxis2 = *oCmb;
	//oCmb = xml.QueryOpt<int>(strXmlRoot+"convofit/minimiser"); if(oCmb) cfg.minimiser = *oCmb;

	// string values
	boost::optional<std::string> osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/crys"); if(osVal) cfg.crys = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/instr"); if(osVal) cfg.instr = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/sqw"); if(osVal) cfg.sqw = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/sqw_conf"); if(osVal) cfg.sqw_conf = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/scanfile"); if(osVal) cfg.scanfile = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"convofit/counter"); if(osVal) cfg.counter = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"convofit/monitor"); if(osVal) cfg.monitor = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"convofit/temp_override"); if(osVal) cfg.temp_override = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"convofit/field_override"); if(osVal) cfg.field_override = *osVal;
	osVal = xml.QueryOpt<std::string>(strXmlRoot+"monteconvo/autosave"); if(osVal) cfg.autosave = *osVal;
	//osVal = xml.QueryOpt<std::string>(strXmlRoot+"convofit/sqw_params"); if(osVal) cfg.sqw_params = *osVal;

	return cfg;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
static bool start_convo_1d(const ConvoConfig& cfg)
{
	//createSqwModel(qstrSqwConf);
	return true;
}


static bool start_convo_2d(const ConvoConfig& cfg)
{
	//createSqwModel(qstrSqwConf);
	return true;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// main program

int monteconvo_main(int argc, char** argv)
{
	try
	{
#ifdef MONTECONVO_STANDALONE	// only show copyright banner if not already displayed from Takin main program
		tl::log_info("--------------------------------------------------------------------------------");
		tl::log_info("This is the Takin command-line convolution simulator (monteconvo), version " TAKIN_VER ".");
		tl::log_info("Written by Tobias Weber <tweber@ill.fr>, 2014 - 2020.");
		tl::log_info(TAKIN_LICENSE("Takin/Monteconvo"));
		tl::log_debug("Resolution calculation uses ", sizeof(t_real_reso)*8, " bit ", tl::get_typename<t_real_reso>(), "s.");
		tl::log_info("--------------------------------------------------------------------------------");
#endif

		load_sqw_plugins();


		// --------------------------------------------------------------------
		// get job files and program options
		std::vector<std::string> vecJobs;

		// normal args
        opts::options_description args("monteconvo options (overriding config file settings)");
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("job-file",
			opts::value<decltype(vecJobs)>(&vecJobs),
			"convolution config file")));
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("max-threads",
			opts::value<decltype(g_iMaxThreads)>(&g_iMaxThreads),
			"maximum number of threads")));
		// dummy arg if launched from takin executable
#ifndef MONTECONVO_STANDALONE
		bool bStartedFromTakin = 0;
		args.add(boost::shared_ptr<opts::option_description>(
			new opts::option_description("monteconvo-cli",
			opts::bool_switch(&bStartedFromTakin),
			"launch convofit from takin")));
#endif


		// positional args
		opts::positional_options_description args_pos;
		args_pos.add("job-file", -1);

		opts::basic_command_line_parser<char> clparser(argc, argv);
		clparser.options(args);
		clparser.positional(args_pos);
		opts::basic_parsed_options<char> parsedopts = clparser.run();

		opts::variables_map opts_map;
		opts::store(parsedopts, opts_map);
		opts::notify(opts_map);


		int args_to_ignore = 1;	// started with "monteconvo-cli"
		if(bStartedFromTakin)
			++args_to_ignore;	// started with "takin --monteconvo-cli"
		if(argc <= args_to_ignore)
		{
			std::ostringstream ostrHelp;
			ostrHelp << "Usage: ";
			for(int argidx=0; argidx<args_to_ignore; ++argidx)
				ostrHelp << argv[argidx] << " ";
			ostrHelp << "[options] <config file>\n";
			ostrHelp << args;
			tl::log_info(ostrHelp.str());
			return -1;
		}

		if(vecJobs.size() == 0)
		{
			tl::log_err("No config files given.");
			return -1;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// load convolution job file
		const std::string& strJobFile = vecJobs[0];
		if(!tl::file_exists(strJobFile.c_str()))
		{
			tl::log_err("Convolution config file \"", strJobFile, "\" does not exist.");
			return -1;
		}

		// add the location of the convo file as a possible global path
		std::string strGlobDir = tl::get_dir(strJobFile);
		clear_global_paths();
		if(strGlobDir != "")
			add_global_path(strGlobDir);


		tl::Prop<std::string> xml;
		if(!xml.Load(strJobFile, tl::PropType::XML))
		{
			tl::log_err("Convolution config file \"", strJobFile, "\" could not be loaded.");
			return -1;
		}


		ConvoConfig cfg = load_config(xml);
		// --------------------------------------------------------------------



		tl::Stopwatch<t_real> watch;
		watch.start();

		bool ok = 0;
		if(cfg.scan_2d)
		{
			tl::log_info("Performing a 2d convolution simulation.");
			ok = start_convo_2d(cfg);
		}
		else
		{
			tl::log_info("Performing a 1d convolution simulation.");
			ok = start_convo_1d(cfg);
		}

		if(!ok)
			tl::log_err("Simulation failed!");


		watch.stop();
		tl::log_info("================================================================================");
		tl::log_info("Start time:     ", watch.GetStartTimeStr());
		tl::log_info("Stop time:      ", watch.GetStopTimeStr());
		tl::log_info("Execution time: ", tl::get_duration_str_secs<t_real>(watch.GetDur()));
		tl::log_info("================================================================================");

		return ok ? 0 : -1;
    }
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

    return 0;
}
// ----------------------------------------------------------------------------
