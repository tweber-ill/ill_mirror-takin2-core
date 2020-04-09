find_library(Dl_LIBRARIES
	NAMES libdl.so
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "DL library"
)

if(Dl_LIBRARIES)
	set(Dl_FOUND 1)
	message("DL library: ${Dl_LIBRARIES}")
else()
	set(Dl_FOUND 0)
	message("Error: libdl could not be found!")
endif()
