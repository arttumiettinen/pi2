
CXXFLAGS+=-fPIC -I. -I../itl2
CFLAGS+=-fPIC -I.
LDFLAGS+=-shared -L$(BUILD_ROOT)/../itl2
LDLIBS+=-litl2 -lfftw3f -lstdc++fs -lpng -ltiff -ljpeg $(OPENCL_LIB)

all: pilib

pilib: $(BUILD_ROOT)/libpi.so

$(BUILD_ROOT)/libpi.so: $(BUILD_ROOT)/../itl2/libitl2.a

include ../easymake.mk
