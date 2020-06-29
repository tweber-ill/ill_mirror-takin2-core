#!/bin/bash
#
# gets dependent repos
# @author Tobias Weber <tweber@ill.fr>
# @date 9-apr-20
# @license GPLv2
#

if [ -L tlibs ] || [ -d tlibs ]; then
        echo -e "A tlibs directory already exists. Skipping.";
else
	echo -e "Cloning tlibs..."
	git clone https://code.ill.fr/scientific-software/takin/tlibs.git
fi


if [ -L data ] || [ -d data ]; then
        echo -e "A data directory already exists. Skipping.";
else
	echo -e "Cloning data repo..."
	git clone https://code.ill.fr/scientific-software/takin/data.git
fi



if [ -L tlibs2 ] || [ -d tlibs2 ]; then
        echo -e "A tlibs2 directory already exists. Skipping.";
else
	echo -e "Cloning tlibs2..."
	git clone https://code.ill.fr/scientific-software/takin/tlibs2.git
fi


if [ -L mag-core ] || [ -d mag-core ]; then
        echo -e "A mag-core directory already exists. Skipping.";
else
	echo -e "Cloning data mag-core tools repo..."
	git clone https://code.ill.fr/scientific-software/takin/mag-core.git
fi
