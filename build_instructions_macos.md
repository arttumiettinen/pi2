Building pi2 in MacOS
---------------------

**Note**
These instructions have been tested on a newly installed MacOS Monterey system running on Apple silicon.
Macs with Intel processor will require at least changes to homebrew-related installation paths, and possibly more.

**Note**
The instructions below guide through installing `gcc` compiler. The default `clang` compiler does not seem to support OpenMP.


In order to build anything, you will need build tools. These are installed by opening a Terminal window and running
```
git
```
If the build tools are missing, the system will guide you through the installation.

Next, let's clone the pi2 repository using Terminal commands
```
mkdir dev
cd dev
mkdir pi2
git clone https://github.com/arttumiettinen/pi2
cd pi2
git checkout experimental
```

Install [homebrew package manager](https://brew.sh) according to its installation instructions (one Terminal command).

Install gcc compiler using Terminal command
```
brew install gcc
```

Now folder `/opt/homebrew/bin` should be in the `PATH` variable. Check that by running
```
echo $PATH
```

Make `gcc` command point to the correct version of `gcc` by running
```
cd /opt/homebrew/bin
ln -s gcc-12 gcc
ln -s g++-12 g++
ln -s c++-12 c++
ln -s cpp-12 cpp
cd ~/dev/pi2
```

**Note**
You might need to change the `12` suffix to some other number depending on the current gcc version.

Install the required libraries using commands
```
brew install fftw
brew install libpng
brew install libtiff
```

Optionally also
```
brew install opencl-headers
brew install opencl-clhpp-headers
```

Finally, compile pi2
```
make NO_OPENCL=1
```

The output will be placed in sub-folder `bin-macos`.

**NOTE**
Compiling with OpenCL does not seem to work at the moment.
