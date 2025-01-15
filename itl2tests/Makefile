
CXXFLAGS+=-I../itl2 -I../fftw-3.3.7-linux64/include
LDFLAGS+=-L$(BUILD_ROOT)/../itl2
LDLIBS+=-litl2 -lfftw3f -lfftw3f_threads -lstdc++fs -lpng -ltiff -ljpeg $(OPENCL_LIB)

EXTRA_DEPS = ../intermediate/$(CONFIG)/itl2/libitl2.a

include ../easymake.mk

