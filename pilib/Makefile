
CXXFLAGS+=-fPIC -I. -I../itl2 -I../fftw-3.3.7-linux64/include
CFLAGS+=-fPIC -I.
LDFLAGS+=-shared -L./../itl2/bin -L./../fftw-3.3.7-linux64/lib
LDLIBS=-lfftw3f -litl2 -lstdc++fs -lpng -ltiff

all: bin/libpilib.so

bin/libpilib.so: ./../itl2/bin/libitl2.a

include ../easymake.mk
