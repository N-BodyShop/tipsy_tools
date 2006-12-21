noinst_LIBRARIES = libsiman.a

libsiman_a_SOURCES = basesimsnap.cpp columngrid.cpp endian.cpp \
		filter.cpp gadgetfile.cpp \
		geometry.cpp grid.cpp metric.cpp \
		particle.cpp rational.cpp scripted.cpp \
		simanexception.cpp simanobject.cpp simsnap.cpp \
		sphkernel.cpp subset.cpp \
		tipsyfile.cpp union.cpp units.cpp value.cpp \
		siman.cpp columnlist.cpp

EXTRA_DIST = endian.hpp

libsiman_adir = $(pkgincludedir)/base/
libsiman_ldadd = $(EXTRALIBS)

libsiman_a_HEADERS = siman.hpp basesimsnap.hpp columngrid.hpp \
		filter.hpp gadgetfile.hpp geometry.hpp grid.hpp \
		metric.hpp particle.hpp rational.hpp scripted.hpp \
		simanexception.hpp simanobject.hpp simsnap.hpp \
		sphkernel.hpp subset.hpp\
		tipsyfile.hpp union.hpp units.hpp value.hpp \
		siman.hpp columnlist.hpp coordinate.hpp 
