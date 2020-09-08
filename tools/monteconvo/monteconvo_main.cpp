/**
 * convolution simulation (CLI)
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

#include "tlibs/log/debug.h"
#include "tlibs/time/stopwatch.h"
#include "tlibs/helper/thread.h"


namespace asio = boost::asio;
namespace sys = boost::system;
namespace opts = boost::program_options;

using t_real = t_real_reso;


// ----------------------------------------------------------------------------
// main program

int monteconvo_main(int argc, char** argv)
{
	try
	{
		tl::log_info("--------------------------------------------------------------------------------");
		tl::log_info("This is the Takin command-line convolution simulator (monteconvo), version " TAKIN_VER ".");
		tl::log_info("Written by Tobias Weber <tweber@ill.fr>, 2014 - 2020.");
		tl::log_info(TAKIN_LICENSE("Takin/Monteconvo"));
		tl::log_debug("Resolution calculation uses ", sizeof(t_real_reso)*8, " bit ", tl::get_typename<t_real_reso>(), "s.");
		tl::log_info("--------------------------------------------------------------------------------");


		load_sqw_plugins();


        // TODO
    }
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

    return 0;
}
// ----------------------------------------------------------------------------
