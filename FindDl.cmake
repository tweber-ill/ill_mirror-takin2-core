find_library(Dl_LIBRARIES
	NAMES libdl.so
	HINTS /usr/lib64/ /usr/lib/ /usr/local/lib64 /usr/local/lib /opt/local/lib
	DOC "DL library"
)

if(Dl_LIBRARIES)
	set(Dl_FOUND 1)
	message("DL library: ${Dl_LIBRARIES}")
else()
	set(Dl_FOUND 0)
	message("Error: libdl could not be found!")
endif()
