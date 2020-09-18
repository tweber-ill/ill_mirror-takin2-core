#!/bin/bash
#
# creates a DEB distro
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
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
mkdir -p ${INSTDIR}/usr/share/applications
mkdir -p ${INSTDIR}/DEBIAN


# control file
echo -e "Package: takin\nVersion: 2.0.0" > ${INSTDIR}/DEBIAN/control
echo -e "Description: inelastic neutron scattering software" >> ${INSTDIR}/DEBIAN/control
echo -e "Maintainer: n/a" >> ${INSTDIR}/DEBIAN/control
echo -e "Architecture: $(dpkg --print-architecture)" >> ${INSTDIR}/DEBIAN/control
echo -e "Section: base\nPriority: optional" >> ${INSTDIR}/DEBIAN/control
echo -e "Depends: libstdc++6, libboost-system1.71.0 (>=1.71.0), libboost-filesystem1.71.0 (>=1.71.0), libboost-iostreams1.71.0 (>=1.71.0), libboost-regex1.71.0 (>=1.71.0), libboost-program-options1.71.0 (>=1.71.0), libboost-python1.71.0 (>=1.71.0), libqt5core5a (>=5.9.5), libqt5gui5 (>=5.9.5), libqt5opengl5 (>=5.9.5), libqt5svg5 (>=5.9.5), libqt5xml5 (>=5.9.5), libqwt-qt5-6 (>=6.1.3), libpython3.8 (>=3.8.0), libfreetype6, python3.8 (>=3.8.0), gnuplot, gnuplot-qt, libopengl0\n" >> ${INSTDIR}/DEBIAN/control


# copy program files
cp -v bin/takin			${INSTDIR}/usr/local/bin
cp -v bin/convofit		${INSTDIR}/usr/local/bin
cp -v bin/convoseries		${INSTDIR}/usr/local/bin
cp -v bin/sfact			${INSTDIR}/usr/local/bin
cp -v bin/takinmod_py		${INSTDIR}/usr/local/bin
cp -v bin/takinmod_jl		${INSTDIR}/usr/local/bin

cp -rv res/*			${INSTDIR}/usr/local/share/takin/res/
cp -rv doc/* 			${INSTDIR}/usr/local/share/takin/res/doc/
cp -v COPYING			${INSTDIR}/usr/local/share/takin
cp -v LICENSES			${INSTDIR}/usr/local/share/takin
cp -v LITERATURE		${INSTDIR}/usr/local/share/takin
cp -v AUTHORS			${INSTDIR}/usr/local/share/takin
cp -v setup_lin/takin.desktop		${INSTDIR}/usr/share/applications
cp -v /usr/local/lib/libMinuit2.so ${INSTDIR}/usr/local/lib

pushd ${INSTDIR}/usr/local/lib
ln -sf libMinuit2.so libMinuit2.so.0
ln -sf libMinuit2.so libMinuit2.so.0.0
ln -sf libMinuit2.so libMinuit2.so.0.0.0
popd


# copy optional external programs
cp -v bin/takin_cif2xml		${INSTDIR}/usr/local/bin
cp -v bin/takin_findsg		${INSTDIR}/usr/local/bin
cp -v bin/takin_pol		${INSTDIR}/usr/local/bin


# permissions
chmod a+x ${INSTDIR}/usr/local/bin/*

# stripping
strip -v ${INSTDIR}/usr/local/bin/*
strip -v ${INSTDIR}/usr/local/lib/*


# startup script
cp -v takin.sh			${INSTDIR}/usr/local/bin



# build package
cd ${INSTDIR}/..
chmod -R 775 ${INSTDIR}
dpkg --build takin
