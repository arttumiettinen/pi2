
# NOTE: Use "make" or "make all" to do a normal build.
#       Use "make NO_OPENCL=1" or "make all NO_OPENCL=1" to do a build without OpenCL support.

CXXFLAGS := -fopenmp -O3 -std=c++17 -fvisibility=hidden -march=native
LDFLAGS := -fopenmp

OPENCL_LIB := -lOpenCL
CONFIG := release

ifdef NO_OPENCL
    CXXFLAGS += -DNO_OPENCL
    CONFIG = release-nocl
    OPENCL_LIB=
endif


TEMP_DIR := $(shell pwd)/intermediate
BUILD_ROOT = $(TEMP_DIR)/$(CONFIG)/$@

export CXXFLAGS
export LDFLAGS
export OPENCL_LIB
export CONFIG
export TEMP_DIR 
export BUILD_ROOT


.PHONY: all clean itl2 pilib pi2 itl2tests

all: itl2tests itl2 pilib pi2
	mkdir -p bin-linux64/$(CONFIG)
	cp ./intermediate/$(CONFIG)/pilib/libpilib.so ./bin-linux64/$(CONFIG)/
	cp ./intermediate/$(CONFIG)/pi2/pi2 ./bin-linux64/$(CONFIG)/
	cp ./python_scripts/*/*.py ./bin-linux64/$(CONFIG)/
	chmod +x ./bin-linux64/$(CONFIG)/*.py
	cp ./example_config/*.txt ./bin-linux64/$(CONFIG)/
	cp ./example_config/*.sh ./bin-linux64/$(CONFIG)/
	chmod +x ./bin-linux64/$(CONFIG)/*.sh
	cp ./example_config/*.cmd ./bin-linux64/$(CONFIG)/
	cp ./LICENSE.txt ./bin-linux64/$(CONFIG)/

clean: itl2tests itl2 pilib pi2

itl2:
	$(MAKE) -C $@ $(MAKECMDGOALS)

pilib: itl2
	./pilib/create_commit_info.sh
	$(MAKE) -C $@ $(MAKECMDGOALS)

pi2: pilib
	$(MAKE) -C $@ $(MAKECMDGOALS)

itl2tests: itl2
	$(MAKE) -C $@ $(MAKECMDGOALS)

