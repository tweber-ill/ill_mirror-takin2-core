/**
 * crystal libs
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2017
 * @license GPLv2
 */

#include "version.h"
#include <string>

namespace xtl {


const char* get_libcrystal_version()
{
	return LIBCRYSTAL_VERSION;
}

const char* get_libcrystal_infos()
{
	return "This is the crystal library.\n"
		"Written by Tobias Weber <tobias.weber@tum.de>, 2017.\n"
		"License: GPLv2.";
}

/**
 * Check if supplied string matches with compiled one
 * (to check if .so library and header files match)
 */
bool check_libcrystal_version(const char* pcHdrVer)
{
	return std::string(LIBCRYSTAL_VERSION) == std::string(pcHdrVer);
}


}
