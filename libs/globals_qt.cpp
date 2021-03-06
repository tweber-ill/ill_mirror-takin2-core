/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
 * @license GPLv2
 */

#include "globals.h"
#include "globals_qt.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"


QFont g_fontGen("DejaVu Sans",10);
QFont g_fontGfx("DejaVu Sans",10);
//QFont g_fontGL("DejaVu Sans Mono",10);
std::string g_strFontGL("res/fonts/DejaVuSansMono.ttf");
int g_iFontGLSize = 24;


QIcon load_icon(const std::string& strIcon)
{
	std::string strFile = find_resource(strIcon);
	if(strFile != "")
		return QIcon(strFile.c_str());

	return QIcon();
}


QPixmap load_pixmap(const std::string& strIcon)
{
	std::string strFile = find_resource(strIcon);
	if(strFile != "")
		return QPixmap(strFile.c_str());

	return QPixmap();
}
