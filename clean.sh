#!/bin/bash
#
# cleans temporary files
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#

find bin -regex 'bin/[_a-zA-Z0-9]*' | xargs rm -f
rm -f bin/*.exe
rm -f plugins/*.so
rm -f plugins/*.dll
rm -f obj/*.o
rm -f ui/*.h
rm -f *.moc
rm -f libs/*.moc
rm -f tools/taz/*.moc
rm -f tools/res/*.moc
rm -f tools/scanviewer/*.moc
rm -f tools/scanpos/*.moc
rm -f tools/powderfit/*.moc
rm -f tools/monteconvo/*.moc
rm -f tools/sglist/*.moc
rm -f dialogs/*.moc

rm -f doc/takin.qch
rm -f doc/takin.qhc
rm -rf doc/devel/html
rm -f doc/devel/*.tmp

if [ -f Makefile ]
then
	echo -e "\nCleaning stuff made by Makefile..."
	make clean
fi

rm -f CMakeCache.txt
rm -rf CMakeFiles
