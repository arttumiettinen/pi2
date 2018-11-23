TOPTARGETS := all clean

SUBDIRS := itl2 pilib pi2 itl2tests


CXXFLAGS := -O3 -fopenmp -std=c++1z 
LDFLAGS := -fopenmp

export CXXFLAGS
export LDFLAGS

$(TOPTARGETS): $(SUBDIRS)
# TODO: This copying should be done using suitable targets...
	mkdir -p bin-linux64/release
	cp ./pilib/bin/libpilib.so ./bin-linux64/release/
	cp ./pi2/bin/pi2 ./bin-linux64/release/
	cp ./python_scripts/generic/*.py ./bin-linux64/release/
	chmod +x ./bin-linux64/release/*.py
	cp ./example_config/*.txt ./bin-linux64/release/
	cp ./LICENSE.txt ./bin-linux64/release/

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)
