/**
 * crystal libs
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2017
 * @license GPLv2
 */

#ifndef __CRYSTALLIBS_VER_H__
#define __CRYSTALLIBS_VER_H__

#define LIBCRYSTAL_VERSION "0.5"

namespace xtl {

extern const char* get_libcrystal_version();
extern const char* get_libcrystal_infos();
extern bool check_libcrystal_version(const char* pcHdrVer);

}
#endif
