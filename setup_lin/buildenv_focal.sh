#!/bin/bash
#
# install all packages needed for building
# @author Tobias Weber <tobias.weber@tum.de>
# @date 28-jul-20
# @license GPLv2
#

if [[ $(id -u) -gt 0 ]]; then
	echo -e "Please run this script as root."
	exit -1
fi



# -----------------------------------------------------------------------------
# install packages
# -----------------------------------------------------------------------------
if ! apt-get install cmake clang build-essential \
	libboost-all-dev libclipper-dev \
	qt5-default qttools5-dev-tools libqt5svg5-dev qt5-assistant \
	libqwt-qt5-dev libpython3-dev \
	libfreetype6-dev libbz2-dev wget coreutils
then
	echo -e "Error: Could not install packages necessary for building."
	exit -1
fi
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# install minuit
# -----------------------------------------------------------------------------
minuit_web=http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz
minuit_arch=$(basename ${minuit_web})
minuit_archdir="${minuit_arch%.tar.gz}"

minuit_dir=$(mktemp -d)
pushd $minuit_dir

if ! wget $minuit_web
then
	echo -e "Error: Could not download Minuit."
	exit -1
fi

if ! tar xzvf $minuit_arch; then
	echo -e "Error: Could not extract Minuit."
	exit -1
fi

cd  $minuit_archdir

if ! ./configure --disable-openmp; then
	echo -e "Error: Could not configure Minuit."
	exit -1
fi

maxprocs=$(nproc)

if ! make -j $maxprocs; then
	echo -e "Error: Could not build Minuit."
	exit -1
fi

if ! make install; then
	echo -e "Error: Could not install Minuit."
	exit -1
fi

popd
# -----------------------------------------------------------------------------
