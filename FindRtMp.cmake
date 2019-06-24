find_library(Rt_LIBRARIES
	NAMES librt.so
	HINTS /usr/lib64/ /usr/lib/ /usr/local/lib64 /usr/local/lib /opt/local/lib
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
	HINTS /usr/lib64/ /usr/lib/ /usr/local/lib64 /usr/local/lib /opt/local/lib
	DOC "MP library"
)

if(Mp_LIBRARIES)
	set(Mp_FOUND 1)
	message("MP library: ${Mp_LIBRARIES}")
else()
	set(Mp_FOUND 0)
	message("Error: libgomp could not be found!")
endif()
