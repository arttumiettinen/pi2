
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


ifeq ($(CONFIG), release)
   	CS_CONFIG = Release
else ifeq ($(CONFIG), release-nocl)
   	CS_CONFIG = Release no OpenCL
endif


.PHONY: all clean itl2 pilib pi2 itl2tests pi2cs pi2csWinFormsTest

all: itl2tests itl2 pilib pi2 pi2cs pi2csWinFormsTest
	
	# Construct full distribution to bin-linux64/$(CONFIG) folder
	mkdir -p bin-linux64/$(CONFIG)
	cp ./intermediate/$(CONFIG)/pilib/libpi.so ./bin-linux64/$(CONFIG)/
	cp ./intermediate/$(CONFIG)/pi2/pi2 ./bin-linux64/$(CONFIG)/
	cp ./x64/$(CONFIG)/*.exe ./bin-linux64/$(CONFIG)/ | true
	cp ./x64/$(CONFIG)/*.dll ./bin-linux64/$(CONFIG)/ | true
	cp ./x64/$(CONFIG)/*.exe.config ./bin-linux64/$(CONFIG)/ | true
	cp ./x64/$(CONFIG)/*.xml ./bin-linux64/$(CONFIG)/ | true
	cp ./python_scripts/*.py ./bin-linux64/$(CONFIG)/
	chmod +x ./bin-linux64/$(CONFIG)/*.py
	cp ./example_config/*.txt ./bin-linux64/$(CONFIG)/
	cp ./example_config/*.sh ./bin-linux64/$(CONFIG)/
	chmod +x ./bin-linux64/$(CONFIG)/*.sh
	cp ./example_config/*.cmd ./bin-linux64/$(CONFIG)/
	cp ./LICENSE.txt ./bin-linux64/$(CONFIG)/

	# For ease of use of .NET builds using pi2, construct full distribution also to
	# "x64/$(CS_CONFIG)" folder. That way the .NET builds can be done with msbuild as usual,
	# without worrying about final platform-specific output folder name.
	mkdir -p "./x64/$(CS_CONFIG)/"
	cp ./bin-linux64/$(CONFIG)/*.so "./x64/$(CS_CONFIG)/"
	cp ./bin-linux64/$(CONFIG)/pi2 "./x64/$(CS_CONFIG)/"
	cp ./example_config/*.txt "./x64/$(CS_CONFIG)/"

clean: itl2tests itl2 pilib pi2 pi2cs pi2csWinFormsTest

itl2:
	$(MAKE) -C $@ $(MAKECMDGOALS)

pilib: itl2
	./pilib/create_commit_info.sh
	$(MAKE) -C $@ $(MAKECMDGOALS)

pi2: pilib
	$(MAKE) -C $@ $(MAKECMDGOALS)

itl2tests: itl2
	$(MAKE) -C $@ $(MAKECMDGOALS)

pi2cs: pilib
	$(MAKE) -C $@ $(MAKECMDGOALS)

pi2csWinFormsTest: pi2cs
	$(MAKE) -C $@ $(MAKECMDGOALS)

