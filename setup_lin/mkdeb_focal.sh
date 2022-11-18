#!/bin/bash
#
# creates a DEB distro
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#


# installation directory
INSTDIR="$1"

if [ "${INSTDIR}" = "" ]; then
	INSTDIR=~/tmp/takin
fi


# directories
mkdir -p ${INSTDIR}/usr/local/bin
mkdir -p ${INSTDIR}/usr/local/lib
mkdir -p ${INSTDIR}/usr/local/share/takin/res
mkdir -p ${INSTDIR}/usr/local/share/takin/3rdparty_licenses
mkdir -p ${INSTDIR}/usr/share/applications
mkdir -p ${INSTDIR}/DEBIAN


# control file
echo -e "Package: takin\nVersion: 2.3.0" > ${INSTDIR}/DEBIAN/control
echo -e "Description: inelastic neutron scattering software" >> ${INSTDIR}/DEBIAN/control
echo -e "Maintainer: n/a" >> ${INSTDIR}/DEBIAN/control
echo -e "Architecture: $(dpkg --print-architecture)" >> ${INSTDIR}/DEBIAN/control
echo -e "Section: base\nPriority: optional" >> ${INSTDIR}/DEBIAN/control
echo -e "Depends:" \
	"libstdc++6," \
	"libboost-system1.71.0 (>=1.71.0)," \
	"libboost-filesystem1.71.0 (>=1.71.0)," \
	"libboost-iostreams1.71.0 (>=1.71.0)," \
	"libboost-regex1.71.0 (>=1.71.0)," \
	"libboost-program-options1.71.0 (>=1.71.0)," \
	"libboost-python1.71.0 (>=1.71.0)," \
	"libqt5core5a (>=5.9.5)," \
	"libqt5gui5 (>=5.9.5)," \
	"libqt5opengl5 (>=5.9.5)," \
	"libqt5svg5 (>=5.9.5)," \
	"libqt5xml5 (>=5.9.5)," \
	"qt5-assistant," \
	"libqwt-qt5-6 (>=6.1.3)," \
	"libpython3.8 (>=3.8.0)," \
	"python3.8 (>=3.8.0)," \
	"python3-numpy," \
	"python3-scipy," \
	"libfreetype6," \
	"gnuplot," \
	"gnuplot-qt," \
	"libopengl0," \
	"liblapacke," \
	"libqhull-r7," \
	"libhdf5-103" \
	"libhdf5-cpp-103" \
	"libqcustomplot2.0\n" \
		>> ${INSTDIR}/DEBIAN/control

# libqcustomplot2.0 is only needed by the external programs


# copy program files
cp -v bin/takin			${INSTDIR}/usr/local/bin
cp -v bin/convofit		${INSTDIR}/usr/local/bin
cp -v bin/convoseries		${INSTDIR}/usr/local/bin
cp -v bin/sfact			${INSTDIR}/usr/local/bin
cp -v bin/takinmod_py		${INSTDIR}/usr/local/bin
cp -v bin/takinmod_jl		${INSTDIR}/usr/local/bin

cp -rv res/*			${INSTDIR}/usr/local/share/takin/res/
cp -rv doc/* 			${INSTDIR}/usr/local/share/takin/res/doc/
cp -v *.txt			${INSTDIR}/usr/local/share/takin
cp -rv 3rdparty_licenses/*	${INSTDIR}/usr/local/share/takin/3rdparty_licenses/
cp -v setup_lin/takin.desktop	${INSTDIR}/usr/share/applications
cp -v /usr/local/lib/libMinuit2.so ${INSTDIR}/usr/local/lib

# if we have the minuit so file (i.e. if it's not statically linked),
# create some symbolic links
pushd ${INSTDIR}/usr/local/lib
if [ -e libMinuit2.so ]; then
	ln -sf libMinuit2.so libMinuit2.so.0
	ln -sf libMinuit2.so libMinuit2.so.0.0
	ln -sf libMinuit2.so libMinuit2.so.0.0.0
fi
popd


# copy optional external programs
cp -v bin/takin_cif2xml		${INSTDIR}/usr/local/bin
cp -v bin/takin_findsg		${INSTDIR}/usr/local/bin
cp -v bin/takin_pol		${INSTDIR}/usr/local/bin
cp -v bin/takin_bz		${INSTDIR}/usr/local/bin
cp -v bin/takin_structfact	${INSTDIR}/usr/local/bin
cp -v bin/takin_magstructfact	${INSTDIR}/usr/local/bin
cp -v bin/takin_scanbrowser	${INSTDIR}/usr/local/bin
cp -v bin/takin_magsgbrowser	${INSTDIR}/usr/local/bin
cp -v bin/takin_magdyn		${INSTDIR}/usr/local/bin
cp -v bin/takin_moldyn		${INSTDIR}/usr/local/bin


# permissions
chmod a+x ${INSTDIR}/usr/local/bin/*

# stripping
strip -v ${INSTDIR}/usr/local/bin/*
strip -v ${INSTDIR}/usr/local/lib/*


# startup script
cp -v takin.sh			${INSTDIR}/usr/local/bin


# remove unnecessary files
find ${INSTDIR} -type f -name ".dir" -exec rm -fv {} \; -print
rm -v ${INSTDIR}/usr/local/share/takin/CMakeLists.txt


# build package
cd ${INSTDIR}/..
chmod -R 775 ${INSTDIR}
dpkg --build takin
