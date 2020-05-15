#!/bin/bash
#
# cleans temporary files
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#

echo -e "Cleaning stuff made by themakefile..."
make -f themakefile clean


if [ -f Makefile ]
then
	echo -e "\nCleaning stuff made by Makefile..."
	make clean
fi

rm -f CMakeCache.txt
rm -rf CMakeFiles


rm -f doc/takin.qch
rm -f doc/takin.qhc
rm -rf doc/devel/html
rm -f doc/devel/*.tmp


# restore link
#rm takin
#ln -sf bin/takin
