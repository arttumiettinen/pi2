
CXXFLAGS+=-I../itl2 -I../fftw-3.3.7-linux64/include
LDFLAGS+=-L$(BUILD_ROOT)/../itl2 -L./../fftw-3.3.7-linux64/lib
LDLIBS+=-litl2 -lfftw3f -lfftw3f_threads -lstdc++fs -lpng -ltiff $(OPENCL_LIB)

include ../easymake.mk
