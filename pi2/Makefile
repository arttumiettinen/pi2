
CXXFLAGS+=-I../pilib -I../itl2
LDFLAGS+=-L$(BUILD_ROOT)/../pilib -Wl,-rpath=\$$$${ORIGIN}
LDLIBS+=-lpi $(OPENCL_LIB)

$(BUILD_ROOT)/pi2: $(BUILD_ROOT)/../itl2/libitl2.a $(BUILD_ROOT)/../pilib/libpi.so

include ../easymake.mk
