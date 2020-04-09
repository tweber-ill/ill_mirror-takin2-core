find_path(Clipper_INCLUDE_DIRS
	NAMES clipper.h
	PATH_SUFFIXES clipper
	HINTS /usr/local/include/clipper /usr/local/include /usr/include/clipper /usr/include /opt/local/include/clipper /opt/local/include
	DOC "Clipper include directories"
)

find_library(Clipper_LIBRARIES
	NAMES libclipper-core.so
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "Clipper library"
)

if(Clipper_LIBRARIES)
	set(Clipper_FOUND 1)
	message("Clipper include directories: ${Clipper_INCLUDE_DIRS}")
	message("Clipper library: ${Clipper_LIBRARIES}")
else()
	set(Clipper_FOUND 0)
	set(Clipper_INCLUDE_DIRS "")
	set(Clipper_LIBRARIES "")
	message("Error: Clipper could not be found!")
endif()
