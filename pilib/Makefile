
CXXFLAGS+=-fPIC -I. -I../itl2 -I../fftw-3.3.7-linux64/include
CFLAGS+=-fPIC -I.
LDFLAGS+=-shared -L$(BUILD_ROOT)/../itl2 -L./../fftw-3.3.7-linux64/lib
LDLIBS+=-litl2 -lfftw3f -lstdc++fs -lpng -ltiff $(OPENCL_LIB)

all: pilib

pilib: $(BUILD_ROOT)/libpi.so

$(BUILD_ROOT)/libpi.so: $(BUILD_ROOT)/../itl2/libitl2.a

include ../easymake.mk
