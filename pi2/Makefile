
CXXFLAGS+=-I../pilib -I../itl2
LDFLAGS+=-L$(BUILD_ROOT)/../pilib -Wl,-rpath=\$$$${ORIGIN}
LDLIBS+=-lpi $(OPENCL_LIB)

all: $(BUILD_ROOT)/pi2

include ../easymake.mk
