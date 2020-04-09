#
# Takin
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2013-2020
# @license GPLv2
#

project(takin)
cmake_minimum_required(VERSION 3.0)

# release or debug build?
set(DEBUG FALSE)

# native or mingw build?
set(BUILD_FOR_MINGW FALSE)
#set(BUILD_FOR_MINGW TRUE)

set(USE_PLUGINS TRUE)
set(USE_TR1_FUNCS FALSE)


if(DEBUG)
	message("Doing a DEBUG build.")

	set(CMAKE_VERBOSE_MAKEFILE TRUE)
	set(CMAKE_BUILD_TYPE Debug)

	set(Boost_DEBUG TRUE)
	set(Boost_USE_DEBUG_LIBS TRUE)
	set(Boost_USE_RELEASE_LIBS FALSE)

	add_definitions(-DDEBUG)
	add_definitions(-ggdb)
	add_definitions(-Wall -Wpedantic -Wextra -Wno-unused-local-typedefs -Wno-unused-parameter)
else()
	message("Doing a RELEASE build.")

	set(CMAKE_VERBOSE_MAKEFILE FALSE)
	set(CMAKE_BUILD_TYPE Release)

	set(Boost_DEBUG FALSE)
	set(Boost_USE_DEBUG_LIBS FALSE)
	set(Boost_USE_RELEASE_LIBS TRUE)

	add_definitions(-DNDEBUG -DQT_NO_DEBUG)
	#add_definitions(-march=native -s)
endif()


set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++11 -D_GLIBCXX_USE_CXX11_ABI)
#add_definitions(-fPIC)

# for application plugins
add_definitions("-DBUILD_APPLI")


if(BUILD_FOR_MINGW)
	# stuff needed by mingw
	set(SOCK2 "ws2_32")
	set(USE_GENTAB FALSE)
	set(USE_RtMp FALSE)
	set(SEARCH_PY FALSE)
else()
	set(SOCK2 "")
	set(USE_GENTAB TRUE)
	set(USE_RtMp TRUE)
	set(SEARCH_PY TRUE)
endif()


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "." "${PROJECT_SOURCE_DIR}/tlibs")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/plugins")

set(Boost_FIND_QUIETLY FALSE)
set(BOOST_COMPS program_options iostreams filesystem system)

if(USE_TR1_FUNCS)
	list(APPEND BOOST_COMPS regex)
endif()

if(SEARCH_PY)
	find_package(Python3 COMPONENTS Interpreter Development)
endif()

if(Python3_Interpreter_FOUND AND Python3_Development_FOUND)
	list(APPEND BOOST_COMPS python37)
endif()


find_package(Boost REQUIRED COMPONENTS ${BOOST_COMPS} REQUIRED)
find_package(Threads REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui OpenGL Svg REQUIRED)
find_package(PNG REQUIRED)
find_package(Qwt REQUIRED)
find_package(Julia)
find_package(Minuit2)
find_package(ZLIB)
find_package(BZip2)

if(USE_RtMp)
	find_package(RtMp REQUIRED)
else()
	set(Rt_LIBRARIES "")
	set(Mp_LIBRARIES "")
endif()

if(USE_PLUGINS)
	find_package(Dl REQUIRED)
else()
	set(Dl_LIBRARIES "")
endif()

if(USE_GENTAB)
	find_package(Clipper REQUIRED)
endif()

set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL)
find_package(Freetype)
find_package(X11)


if("${QT_MOC_EXECUTABLE}" STREQUAL "")
	get_target_property(QT_MOC_EXECUTABLE Qt5::moc IMPORTED_LOCATION)
endif()

if("${QT_UIC_EXECUTABLE}" STREQUAL "")
	get_target_property(QT_UIC_EXECUTABLE Qt5::uic IMPORTED_LOCATION)
endif()


message("Project source dir: ${PROJECT_SOURCE_DIR}")
message("CMake source dir: ${CMAKE_SOURCE_DIR}")
message("CMake binary dir: ${CMAKE_BINARY_DIR}")
message("Install path prefix: ${CMAKE_INSTALL_PREFIX}")
message("Using Qt version ${Qt5Core_VERSION}")
if(NOT "${QWT_MAJOR_VERSION}" STREQUAL "")
	message("Using Qwt version ${QWT_MAJOR_VERSION}")
endif()
message("Using Boost version ${Boost_VERSION}")
message("Using Zlib version ${ZLIB_VERSION_STRING}")
message("Using Bzip2 version ${BZIP2_VERSION_STRING}")
if(NOT "${Python3_VERSION}" STREQUAL "")
	message("Using Py version ${Python3_VERSION}")
endif()
message("Qt moc: ${QT_MOC_EXECUTABLE}\nQt uic: ${QT_UIC_EXECUTABLE}")


add_definitions(-DUSE_NET -DUSE_CIF)
add_definitions(-DNO_LAPACK)
add_definitions(-DNO_JPEG -DNO_TIFF)
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
#add_definitions(${Boost_CXX_FLAGS} -DBOOST_FILESYSTEM_NO_DEPRECATED -DBOOST_ERROR_CODE_HEADER_ONLY)
add_definitions(-DTLIBS_INC_HDR_IMPLS)
add_definitions(-DHAS_COMPLEX_ERF)

if(USE_PLUGINS)
	add_definitions(-DUSE_PLUGINS)
endif()

if(USE_TR1_FUNCS)
	add_definitions(-DUSE_BOOST_REX)
endif()

if(NOT BUILD_FOR_MINGW)
	add_definitions(-DUSE_MULTI_TEXTURES)
else()
	add_definitions(-DNO_HELP_ASSISTANT)
	add_definitions(-DNO_TERM_CMDS)
	add_definitions(-DUSE_WIDE_STR)
endif()

# hack to prevent math_fwd.h inclusion
add_definitions(-DBOOST_MATH_FWD_HPP)


#add_definitions("${Qt5Core_DEFINITIONS}")
#add_definitions("${Qt5Widgets_DEFINITIONS}")
#add_definitions("${QtGui_DEFINITIONS}")
#add_definitions("${Qt5OpenGL_DEFINITIONS}")

include_directories("${PROJECT_SOURCE_DIR}" ".")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${PNG_INCLUDE_DIRS}")
#include_directories("${JPEG_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS}")
include_directories("${Qt5Core_INCLUDE_DIRS}" "${Qt5Gui_INCLUDE_DIRS}" "${Qt5Widgets_INCLUDE_DIRS}")
include_directories("${QWT_INCLUDE_DIRS}")
include_directories("${ZLIB_INCLUDE_DIRS}")
include_directories("${BZIP2_INCLUDE_DIR}")



# py scripting support
if(Python3_Interpreter_FOUND AND Python3_Development_FOUND)
	message("Enabling py scripting support.")
	add_definitions(-DUSE_PY)
	include_directories("${Python3_INCLUDE_DIRS}")
else()
	message("Disabling py scripting support.")
	add_definitions(-DNO_PY)
endif()



# jl scripting support
if(Julia_FOUND)
	message("Enabling jl scripting support.")
	add_definitions(-DUSE_JL)
	include_directories("${Julia_INCLUDE_DIRS}")
else()
	message("Disabling jl scripting support.")
	add_definitions(-DNO_JL)
endif()



# fitting support
if(Minuit2_FOUND)
	message("Enabling fitting support.")
	include_directories("${Minuit2_INCLUDE_DIRS}")
	set(LIBS_FIT ${Minuit2_LIBRARIES})
else()
	message("Disabling fitting support.")
	add_definitions(-DNO_FIT)
endif()



if(X11_FOUND)
	include_directories(${X11_X11_INCLUDE_PATH})
else()
	set(X11_X11_LIB "")
endif()



# 3d support
if(OPENGL_FOUND AND FREETYPE_FOUND)
	message("Enabling 3D support. Inc: \"${OPENGL_INCLUDE_DIR}\".")

	include_directories("${FREETYPE_INCLUDE_DIRS}")
	if(NOT "${OPENGL_INCLUDE_DIR}" STREQUAL "")
		include_directories("${OPENGL_INCLUDE_DIR}")
		include_directories("${OPENGL_INCLUDE_DIR}/GL" "${OPENGL_INCLUDE_DIR}/OpenGL")
	else()
		message("Empty OGL dirs, guessing some defaults.")
		include_directories("/usr/include/GL" "/usr/include/OpenGL" "/usr/x86_64-w64-mingw32/sys-root/mingw/include/GL")
	endif()

	add_definitions(-DUSE_3D)

	set(SRCS_3D
		libs/plotgl.cpp libs/plotgl2.cpp tools/taz/recip3d.cpp
		tools/taz/real3d.cpp tools/taz/bz3d.cpp dialogs/EllipseDlg3D.cpp

		tlibs/gfx/gl.cpp
	)
	set(LIBS_3D
		${FREETYPE_LIBRARY}
		Qt5::OpenGL ${OPENGL_LIBRARIES}
		${X11_X11_LIB}
	)
else()
	message("Disabling 3D support: OGL: ${OPENGL_FOUND}, FT: ${FREETYPE_FOUND}")

	add_definitions(-DNO_3D)
endif()



set(CMAKE_AUTOMOC OFF)
set(CMAKE_AUTOUIC OFF)
set(CMAKE_AUTORCC OFF)



# -----------------------------------------------------------------------------
# gentab
# -----------------------------------------------------------------------------
if(USE_GENTAB)
	add_executable(gentab
		tools/gentab/gentab.cpp
		libs/spacegroups/spacegroup_clp.cpp libs/spacegroups/crystalsys.cpp
		tlibs/log/log.cpp
	)

	target_link_libraries(gentab
		${Clipper_LIBRARIES}
		Qt5::Core Qt5::Gui Qt5::Widgets
		Boost::iostreams Boost::system Boost::filesystem Boost::program_options
		${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
	)

	add_custom_command(TARGET gentab POST_BUILD
		COMMAND cd ${CMAKE_BINARY_DIR} && ./bin/gentab
		MAIN_DEPENDENCY gentab)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# takin
# -----------------------------------------------------------------------------
add_executable(takin
	tools/taz/taz_main.cpp
	tools/taz/taz.cpp tools/taz/tas_layout.cpp tools/taz/tof_layout.cpp
	tools/taz/taz_crys.cpp
	tools/taz/taz_file.cpp tools/taz/taz_export.cpp
	tools/taz/tas_layout.cpp tools/taz/scattering_triangle.cpp
	tools/taz/real_lattice.cpp tools/taz/proj_lattice.cpp
	tools/taz/nicos.cpp tools/taz/sics.cpp tools/taz/taz_net.cpp

	tools/res/r0.cpp
	tools/res/cn.cpp tools/res/pop.cpp tools/res/eck.cpp tools/res/viol.cpp tools/res/simple.cpp
	tools/res/ResoDlg.cpp tools/res/ResoDlg_file.cpp

	tools/monteconvo/ConvoDlg.cpp tools/monteconvo/ConvoDlg_file.cpp
	tools/monteconvo/ConvoDlg_sim.cpp tools/monteconvo/ConvoDlg_fit.cpp
	tools/convofit/convofit_import.cpp
	tools/monteconvo/SqwParamDlg.cpp tools/monteconvo/TASReso.cpp
	tools/monteconvo/sqw.cpp tools/monteconvo/sqw_uniform_grid.cpp
	tools/monteconvo/sqwbase.cpp tools/monteconvo/sqwfactory.cpp

	tools/convofit/scan.cpp
	tools/scanviewer/scanviewer.cpp tools/scanviewer/FitParamDlg.cpp
	tools/scanpos/ScanPosDlg.cpp
	tools/powderfit/PowderFitDlg.cpp

	dialogs/SettingsDlg.cpp dialogs/FilePreviewDlg.cpp
	dialogs/GotoDlg.cpp dialogs/DWDlg.cpp dialogs/DynPlaneDlg.cpp
	dialogs/NeutronDlg.cpp dialogs/TOFDlg.cpp
	dialogs/SpurionDlg.cpp dialogs/PowderDlg.cpp
	dialogs/FavDlg.cpp
	dialogs/RecipParamDlg.cpp dialogs/RealParamDlg.cpp
	dialogs/SrvDlg.cpp dialogs/NetCacheDlg.cpp dialogs/ScanMonDlg.cpp
	dialogs/EllipseDlg.cpp tools/ffact/FormfactorDlg.cpp
	dialogs/AtomsDlg.cpp dialogs/AboutDlg.cpp
	dialogs/DeadAnglesDlg.cpp dialogs/LogDlg.cpp

	libs/spacegroups/spacegroup.cpp tools/sglist/SgListDlg.cpp
	libs/globals.cpp libs/globals_qt.cpp libs/spacegroups/crystalsys.cpp
	libs/formfactors/formfact.cpp libs/qt/qthelper.cpp libs/qt/qwthelper.cpp

	${SRCS_3D}

	#tlibs/file/loadinstr.cpp
	libs/qt/recent.cpp  tlibs/file/x3d.cpp
	tlibs/log/log.cpp tlibs/log/debug.cpp
	tlibs/math/rand.cpp
	tlibs/string/spec_char.cpp
	tlibs/helper/flags.cpp  tlibs/version.cpp

	libcrystal/version.cpp

	./3rdparty/Faddeeva.cc
)

target_link_libraries(takin
	${LIBS_FIT}
	${SOCK2}
	Qt5::Core Qt5::Gui Qt5::Svg
	${QWT_LIBRARIES}
	${LIBS_3D}
	${PNG_LIBRARIES}
	${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
	Boost::iostreams Boost::system Boost::filesystem Boost::program_options
	Threads::Threads ${Mp_LIBRARIES} ${Rt_LIBRARIES} ${Dl_LIBRARIES}
)


# pre-build
add_custom_target(prebuild-script
	COMMAND ${PROJECT_SOURCE_DIR}/prebuild.sh ${QT_UIC_EXECUTABLE} ${QT_MOC_EXECUTABLE} ${PROJECT_SOURCE_DIR}
)


# post-build
#add_custom_command(TARGET takin POST_BUILD
#	COMMAND qcollectiongenerator ${PROJECT_SOURCE_DIR}/doc/takin.qhcp -o ${PROJECT_SOURCE_DIR}/doc/takin.qhc
#		&& cp -v ${PROJECT_SOURCE_DIR}/doc/takin.qhc ${PROJECT_SOURCE_DIR}/res/
#		&& cp -v ${PROJECT_SOURCE_DIR}/doc/takin.qch ${PROJECT_SOURCE_DIR}/res/
#	MAIN_DEPENDENCY takin)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET takin POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takin>
		MAIN_DEPENDENCY takin
	)
endif()


add_dependencies(takin prebuild-script)

if(USE_GENTAB)
	add_dependencies(takin gentab)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# build plugins
# -----------------------------------------------------------------------------

if(Python3_Interpreter_FOUND AND Python3_Development_FOUND)
#add_library(takinmod_py SHARED
add_executable(takinmod_py
	tools/monteconvo/sqw_py.cpp tools/monteconvo/sqwbase.cpp
	tlibs/log/log.cpp
	tlibs/math/rand.cpp
)

target_link_libraries(takinmod_py
	Boost::iostreams Boost::system Boost::filesystem
	Boost::python37 ${Python3_LIBRARIES}
	Threads::Threads ${Rt_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET takinmod_py POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takinmod_py>
		MAIN_DEPENDENCY takinmod_py
	)
endif()
endif()



if(Julia_FOUND)
#add_library(takinmod_jl SHARED
add_executable(takinmod_jl
	tools/monteconvo/sqw_jl.cpp tools/monteconvo/sqwbase.cpp
	tlibs/log/log.cpp
	tlibs/math/rand.cpp
)

target_link_libraries(takinmod_jl
	Boost::iostreams Boost::system Boost::filesystem
	${Julia_LIBRARIES}
	Threads::Threads ${Rt_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET takinmod_jl POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takinmod_jl>
		MAIN_DEPENDENCY takinmod_jl
	)
endif()
endif()
# -----------------------------------------------------------------------------



if(Minuit2_FOUND)
# -----------------------------------------------------------------------------
# convofit, convoseries
# -----------------------------------------------------------------------------
add_executable(convofit
	tools/res/r0.cpp
	tools/res/cn.cpp tools/res/pop.cpp tools/res/eck.cpp tools/res/viol.cpp
	#tools/res/simple.cpp

	tools/monteconvo/TASReso.cpp
	tools/monteconvo/sqw.cpp tools/monteconvo/sqw_uniform_grid.cpp
	tools/monteconvo/sqwbase.cpp tools/monteconvo/sqwfactory.cpp

	tools/convofit/convofit.cpp tools/convofit/convofit_import.cpp
	tools/convofit/model.cpp tools/convofit/scan.cpp
	tools/convofit/convofit_main.cpp


	# statically link tlibs externals
	tlibs/log/log.cpp
	tlibs/math/rand.cpp
	tlibs/file/tmp.cpp
	#tlibs/file/loadinstr.cpp
	libs/globals.cpp
)

set_target_properties(convofit PROPERTIES COMPILE_FLAGS "-DNO_QT")

target_link_libraries(convofit
	${LIBS_FIT} ${SOCK2}
	Threads::Threads ${Mp_LIBRARIES} ${Rt_LIBRARIES} ${Dl_LIBRARIES}
	Boost::iostreams Boost::system Boost::filesystem Boost::program_options
	${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET convofit POST_BUILD
		COMMAND strip -v $<TARGET_FILE:convofit>
		MAIN_DEPENDENCY convofit
	)
endif()

# -----------------------------------------------------------------------------

add_executable(convoseries
	tools/convofit/scanseries.cpp

	# statically link tlibs externals
	tlibs/log/log.cpp
)

set_target_properties(convoseries PROPERTIES COMPILE_FLAGS "-DNO_QT")

target_link_libraries(convoseries
	Boost::iostreams Boost::system Boost::filesystem Boost::program_options
	${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_custom_command(TARGET convoseries POST_BUILD
		COMMAND strip -v $<TARGET_FILE:convoseries>
		MAIN_DEPENDENCY convoseries
	)
endif()
# -----------------------------------------------------------------------------

endif()




# -----------------------------------------------------------------------------
# install
# -----------------------------------------------------------------------------
install(TARGETS takin DESTINATION bin COMPONENT apptakin)

if(Minuit2_FOUND)
	install(TARGETS convofit convoseries DESTINATION bin COMPONENT appconvofit)
endif()

if(USE_PLUGINS)
	install(DIRECTORY ${PROJECT_SOURCE_DIR}/plugins DESTINATION share/takin COMPONENT takinplugins)
endif()

if(Python3_Interpreter_FOUND AND Python3_Development_FOUND)
	install(TARGETS takinmod_py DESTINATION bin COMPONENT apptakinmod_py)
endif()

if(Julia_FOUND)
	install(TARGETS takinmod_jl DESTINATION bin COMPONENT apptakinmod_jl)
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples ${PROJECT_SOURCE_DIR}/doc
	DESTINATION share/takin COMPONENT takindata1)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/res/data ${PROJECT_SOURCE_DIR}/res/doc ${PROJECT_SOURCE_DIR}/res/icons
	DESTINATION share/takin/res COMPONENT takindata2)
install(FILES ${PROJECT_SOURCE_DIR}/COPYING ${PROJECT_SOURCE_DIR}/AUTHORS
	${PROJECT_SOURCE_DIR}/LICENSES ${PROJECT_SOURCE_DIR}/LITERATURE
	DESTINATION share/takin COMPONENT takininfo)
install(FILES ${PROJECT_SOURCE_DIR}/takin.sh
	PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE OWNER_READ GROUP_READ WORLD_READ
	DESTINATION bin COMPONENT takinstarter)
install(FILES ${PROJECT_SOURCE_DIR}/takin.desktop DESTINATION share/applications COMPONENT takindesk)
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------
# TODO: fix package
# -----------------------------------------------------------------------------
set(CPACK_PACKAGE_NAME "takin")
set(CPACK_PACKAGE_CONTACT "tweber@ill.fr")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "inelastic neutron scattering package")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_EXECUTABLES "takin;Takin")
set(CPACK_PACKAGE_VERSION "2.0.0")
#set(CPACK_COMPONENTS_ALL apptakin appconvofit takindata1 takindata2 takininfo takinstarter takindesk)
set(CPACK_SET_DESTDIR TRUE)

set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "n/a")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libboost, libqtcore5 (>=5.8.0), libqtgui5 (>=5.8.0), libqtopengl5 (>=5.8.0), libqtsvg5 (>=5.8.0), libqtxml5 (>=5.8.0), libqwt, python, freetype6")

execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_ARCHITECTURE)
message("Deb Arch: ${CPACK_DEBIAN_ARCHITECTURE}")


#include(CPack)
# -----------------------------------------------------------------------------
