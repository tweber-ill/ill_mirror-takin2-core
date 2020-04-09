find_library(Rt_LIBRARIES
	NAMES librt.so
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "RT library"
)

if(Rt_LIBRARIES)
	set(Rt_FOUND 1)
	message("RT library: ${Rt_LIBRARIES}")
else()
	set(Rt_FOUND 0)
	message("Error: librt could not be found!")
endif()



find_library(Mp_LIBRARIES
	NAMES libgomp.so libgomp.so.1
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "MP library"
)

if(Mp_LIBRARIES)
	set(Mp_FOUND 1)
	message("MP library: ${Mp_LIBRARIES}")
else()
	set(Mp_FOUND 0)
	message("Error: libgomp could not be found!")
endif()
