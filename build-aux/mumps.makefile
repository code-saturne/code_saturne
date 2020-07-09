include ${topdir}/Makefile.inc

getincludedirs:
	@echo ${IORDERINGSC} ${INCS}

getlibdirs:
	@echo ${LORDERINGS} ${LIBS} -Wl,--enable-new-dtags

getlinklibs:
	@echo ${LORDERINGS} ${LIBS} ${LIBBLAS} ${LIBOTHERS}
