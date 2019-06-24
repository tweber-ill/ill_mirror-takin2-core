#
# crystal libs
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2014-2017
# @license GPLv2
#

project(xtl)
cmake_minimum_required(VERSION 3.0)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ".")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

find_package(Boost REQUIRED COMPONENTS iostreams system filesystem REQUIRED)
#find_package(Threads REQUIRED)


#message("Project source dir: ${PROJECT_SOURCE_DIR}")
#message("CMake source dir: ${CMAKE_SOURCE_DIR}")
#message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install path prefix: ${CMAKE_INSTALL_PREFIX}")


add_definitions(-DNDEBUG)
add_definitions(-DNO_QT)
add_definitions(-DUSE_IOSTR -DUSE_BOOST_REX)
add_definitions(-DNO_LAPACK)
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(${Boost_CXX_FLAGS})


include_directories("${PROJECT_SOURCE_DIR}" ".")
include_directories("${Boost_INCLUDE_DIRS}/..")


add_definitions(-std=c++11)
#add_definitions(-march=native -s)



# -----------------------------------------------------------------------------
# build
# -----------------------------------------------------------------------------
add_library(xtl SHARED
	version.cpp libs/libxtl.cpp libs/globals.cpp
	libs/formfactors/formfact.cpp
	libs/spacegroups/crystalsys.cpp libs/spacegroups/spacegroup.cpp
	tlibs/log/log.cpp
)


target_link_libraries(xtl ${Boost_LIBRARIES} #Threads::Threads
)
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------
# install
# -----------------------------------------------------------------------------
install(TARGETS xtl  DESTINATION lib)
# -----------------------------------------------------------------------------
