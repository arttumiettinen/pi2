
CXXFLAGS := -fopenmp -O3 -std=c++17 -fvisibility=hidden -march=native
LDFLAGS := -fopenmp

export CXXFLAGS
export LDFLAGS


.PHONY: all clean itl2 pilib pi2 itl2tests

all: itl2tests itl2 pilib pi2
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

