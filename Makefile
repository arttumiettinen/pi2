TOPTARGETS := all clean

SUBDIRS := itl2 pilib pi2 itl2tests

CXXFLAGS := -fopenmp -O3 -std=c++17 -fvisibility=hidden -march=native
LDFLAGS := -fopenmp

export CXXFLAGS
export LDFLAGS

$(TOPTARGETS): $(SUBDIRS)
# TODO: This copying should be done using suitable targets...
	mkdir -p bin-linux64/release
	cp ./pilib/bin/libpilib.so ./bin-linux64/release/
	cp ./pi2/bin/pi2 ./bin-linux64/release/
	cp ./python_scripts/*/*.py ./bin-linux64/release/
	chmod +x ./bin-linux64/release/*.py
	cp ./example_config/*.txt ./bin-linux64/release/
	cp ./example_config/*.sh ./bin-linux64/release/
	chmod +x ./bin-linux64/release/*.sh
	cp ./example_config/*.cmd ./bin-linux64/release/
	cp ./LICENSE.txt ./bin-linux64/release/

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)
