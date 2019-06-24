#!/bin/bash
#
# gets external files
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#


GTAR="$(which gtar)"
if [ $? -ne 0 ]; then
	GTAR="$(which gnutar)"

	if [ $? -ne 0 ]; then
		GTAR="$(which tar)"
	fi
fi


FINDQWT=http://cmake.org/Wiki/images/2/27/FindQwt.cmake
FADD=http://ab-initio.mit.edu/Faddeeva
TANGOICONS=http://tango.freedesktop.org/releases/tango-icon-theme-0.8.90.tar.gz

SCATLENS=https://www.ncnr.nist.gov/resources/n-lengths/list.html
MAGFFACT_J0_1=https://www.ill.eu/sites/ccsl/ffacts/ffactnode5.html
MAGFFACT_J0_2=https://www.ill.eu/sites/ccsl/ffacts/ffactnode6.html
MAGFFACT_J0_3=https://www.ill.eu/sites/ccsl/ffacts/ffactnode7.html
MAGFFACT_J0_4=https://www.ill.eu/sites/ccsl/ffacts/ffactnode8.html
MAGFFACT_J2_1=https://www.ill.eu/sites/ccsl/ffacts/ffactnode9.html
MAGFFACT_J2_2=https://www.ill.eu/sites/ccsl/ffacts/ffactnode10.html
MAGFFACT_J2_3=https://www.ill.eu/sites/ccsl/ffacts/ffactnode11.html
MAGFFACT_J2_4=https://www.ill.eu/sites/ccsl/ffacts/ffactnode12.html

SCATLENS2=https://raw.githubusercontent.com/neutronpy/neutronpy/master/neutronpy/database/scattering_lengths.json
MAGFFACT2=https://raw.githubusercontent.com/neutronpy/neutronpy/master/neutronpy/database/magnetic_form_factors.json

SPACEGROUPS=https://raw.githubusercontent.com/egonw/bodr/master/bodr/crystal/space-groups.xml
ELEMENTS=https://raw.githubusercontent.com/egonw/bodr/master/bodr/elements/elements.xml


#
# qwt finder for cmake
#
function dl_findqwt
{
#	rm -f FindQwt.cmake

	if [ ! -f FindQwt.cmake ]; then
		echo -e "Downloading FindQwt...\n"

		if ! wget ${FINDQWT}; then
			echo -e "Error: Cannot download FindQwt.";
			return -1;
		fi
	fi
}



#
# faddeeva module
#
function dl_fadd
{
	if [ ! -f 3rdparty/Faddeeva.hh  ]; then
		echo -e "Downloading Faddeeva library...\n"

		if ! (wget ${FADD}.hh -O 3rdparty/Faddeeva.hh &&
			wget ${FADD}.cc -O 3rdparty/Faddeeva.cc); then
			echo -e "Error: Cannot download Faddeeva library.";
			return -1;
		fi
	fi
}



#
# tango icons
#
function dl_tangoicons
{
#	rm -f tmp/tango-icon-theme.tar.gz

	if [ ! -f tmp/tango-icon-theme.tar.gz  ]; then
		echo -e "Downloading Tango icons...\n"

		if ! wget ${TANGOICONS} -O tmp/tango-icon-theme.tar.gz; then
			echo -e "Error: Cannot download Tango icons.";
			return -1;
		fi
	fi


	echo -e "Extracting Tango icons...\n"
	cd tmp

	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/document-save.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/document-save-as.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/document-open.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/document-new.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/list-add.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/list-remove.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/system-log-out.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/view-refresh.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/media-playback-start.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/actions/media-playback-stop.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/categories/preferences-system.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/devices/media-floppy.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/devices/drive-harddisk.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/devices/network-wireless.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/mimetypes/image-x-generic.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/mimetypes/x-office-spreadsheet-template.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/apps/accessories-calculator.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/status/network-transmit-receive.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/status/network-offline.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/status/dialog-information.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/status/weather-snow.svg --strip-components=3
	${GTAR} --wildcards -xzvf tango-icon-theme.tar.gz */scalable/apps/help-browser.svg --strip-components=3

	cd ..
	mv -v tmp/*.svg res/icons/
}



#
# periodic table of elements
#
function dl_elements
{
	if [ ! -f tmp/${ELEMENTS##*/} ]; then
		echo -e "Obtaining periodic table of elements...\n"

		if ! wget ${ELEMENTS} -O tmp/${ELEMENTS##*/}; then
			echo -e "Error: Cannot download periodic table.";
			return -1;
		fi

		if [ ! -f tmp/${ELEMENTS##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following file manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${ELEMENTS}\n"
			echo -e "Afterwards, please run this script again.\n"
		fi
	fi
}



#
# space groups
#
function dl_spacegroups
{
	if [ ! -f tmp/${SPACEGROUPS##*/} ]; then
		echo -e "Obtaining space group tables...\n"

		if ! wget ${SPACEGROUPS} -O tmp/${SPACEGROUPS##*/}; then
			echo -e "Error: Cannot download space group table.";
			return -1;
		fi

		if [ ! -f tmp/${SPACEGROUPS##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following file manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${SPACEGROUPS}\n"
			echo -e "Afterwards, please run this script again.\n"
		fi
	fi
}



#
# scattering lengths
#
function dl_scatlens
{
	if [ ! -f tmp/scatlens.html ]; then
		echo -e "Obtaining scattering length list...\n"

#		if ! wget ${SCATLENS} -O tmp/${SCATLENS##*/}; then
#			echo -e "Error: Cannot download scattering length list.";
#			return -1;
#		fi

		if [ ! -f tmp/${SCATLENS##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following file manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${SCATLENS}\n"
			echo -e "Afterwards, please run this script again.\n"
		else
			cp -v tmp/${SCATLENS##*/} tmp/scatlens.html
		fi
	fi
}



#
# scattering lengths (alternative table)
#
function dl_scatlens2
{
	if [ ! -f tmp/${SCATLENS2##*/} ]; then
		echo -e "Obtaining scattering length table...\n"

		if ! wget ${SCATLENS2} -O tmp/${SCATLENS2##*/}; then
			echo -e "Error: Cannot download scattering length table.";
			return -1;
		fi

		if [ ! -f tmp/${SCATLENS2##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following file manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${SCATLENS2}\n"
			echo -e "Afterwards, please run this script again.\n"
		fi
	fi
}



#
# magnetic form factors
#
function dl_magffacts
{
	if [ ! -f tmp/j2_4.html ]; then
		echo -e "Obtaining magnetic form factor lists...\n"

#		if ! (wget ${MAGFFACT_J0_1} -O tmp/${MAGFFACT_J0_1##*/} &&
#			wget ${MAGFFACT_J0_2} -O tmp/${MAGFFACT_J0_2##*/} &&
#			wget ${MAGFFACT_J0_3} -O tmp/${MAGFFACT_J0_3##*/} &&
#			wget ${MAGFFACT_J0_4} -O tmp/${MAGFFACT_J0_4##*/} &&
#
#			wget ${MAGFFACT_J2_1} -O tmp/${MAGFFACT_J2_1##*/} &&
#			wget ${MAGFFACT_J2_2} -O tmp/${MAGFFACT_J2_2##*/} &&
#			wget ${MAGFFACT_J2_3} -O tmp/${MAGFFACT_J2_3##*/} &&
#			wget ${MAGFFACT_J2_4} -O tmp/${MAGFFACT_J2_4##*/}); then
#			echo -e "Error: Cannot download magnetic form factor lists.";
#			return -1;
#		fi

		if [ ! -f tmp/${MAGFFACT_J2_4##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following files manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${MAGFFACT_J0_1}"
			echo -e "\t${MAGFFACT_J0_2}"
			echo -e "\t${MAGFFACT_J0_3}"
			echo -e "\t${MAGFFACT_J0_4}"
			echo -e "\t${MAGFFACT_J2_1}"
			echo -e "\t${MAGFFACT_J2_2}"
			echo -e "\t${MAGFFACT_J2_3}"
			echo -e "\t${MAGFFACT_J2_4}\n"
			echo -e "Afterwards, please run this script again.\n"
		else
			cp -v tmp/${MAGFFACT_J0_1##*/} tmp/j0_1.html
			cp -v tmp/${MAGFFACT_J0_2##*/} tmp/j0_2.html
			cp -v tmp/${MAGFFACT_J0_3##*/} tmp/j0_3.html
			cp -v tmp/${MAGFFACT_J0_4##*/} tmp/j0_4.html
			cp -v tmp/${MAGFFACT_J2_1##*/} tmp/j2_1.html
			cp -v tmp/${MAGFFACT_J2_2##*/} tmp/j2_2.html
			cp -v tmp/${MAGFFACT_J2_3##*/} tmp/j2_3.html
			cp -v tmp/${MAGFFACT_J2_4##*/} tmp/j2_4.html
		fi
	fi
}



#
# magnetic form factors (alternative table)
#
function dl_magffacts2
{
	if [ ! -f tmp/${MAGFFACT2##*/} ]; then
		echo -e "Obtaining magnetic form factors table...\n"

		if ! wget ${MAGFFACT2} -O tmp/${MAGFFACT2##*/}; then
			echo -e "Error: Cannot download magnetic form factors table.";
			return -1;
		fi

		if [ ! -f tmp/${MAGFFACT2##*/} ]; then
			echo -e "Files cannot be automatically downloaded. "
			echo -e "Please download the following file manually and place them in the ./tmp subfolder:\n"
			echo -e "\t${MAGFFACT2}\n"
			echo -e "Afterwards, please run this script again.\n"
		fi
	fi
}



mkdir tmp
echo -e "--------------------------------------------------------------------------------"
dl_elements
dl_spacegroups
echo -e "--------------------------------------------------------------------------------"
#dl_scatlens
dl_scatlens2
echo -e "--------------------------------------------------------------------------------"
#dl_magffacts
dl_magffacts2
echo -e "--------------------------------------------------------------------------------"
dl_fadd
echo -e "--------------------------------------------------------------------------------"
dl_tangoicons
echo -e "--------------------------------------------------------------------------------"
dl_findqwt
echo -e "--------------------------------------------------------------------------------"

echo -e "\nAfter successfully obtaining all resource files, the program ./gentab has to be run!\n"
