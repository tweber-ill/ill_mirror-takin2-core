#!/bin/bash
#
# @date 16-dec-2020
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# cleans up the python distribution
#

PRG="takin.app"

# remove cached files
find ${PRG}/Contents/Frameworks/Python.framework -type d -name "__pycache__" -exec rm -rfv {} \;
find ${PRG}/Contents/Frameworks/Python.framework -type f -name "*.pyc" -exec rm -fv {} \;


# remove tcl/tk
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/libtk*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/libtcl*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl8*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tk8*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/itcl4*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/pkgconfig/tcl.pc
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tk*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/pkgconfig/tk.pc
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/python3.9/tkinter
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/python3.9/lib-dynload/_tkinter*.so
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/Tk*


# remove non-needed site packages
pushd ${PRG}/Contents
ln -sf Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages
popd

rm -rfv ${PRG}/Contents/site-packages/setuptools*
rm -rfv ${PRG}/Contents/site-packages/pip*
rm -rfv ${PRG}/Contents/site-packages/easy*
rm -rfv ${PRG}/Contents/site-packages/pkg_resources


# strip libraries
find ${PRG}/Contents/Frameworks/Python.framework -name "*.so" -exec strip -v {} \; -print
find ${PRG}/Contents/Frameworks/Python.framework -name "*.dylib" -exec strip -v {} \; -print


# app
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/Resources/Python.app
