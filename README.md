
# itl2 and pi2

This repository contains
* itl2 project, a C++ Image analysis Template Library (version 2). It contains a small number of parallelized algorithms for processing 3D volume images e.g. from computed tomography scanners.
* pi2 program (Process Image 2) that exposes functionality of itl2 library for end users. It also adds the possibility to distribute processing of large images to a cluster.
* pi2py and pi2cs libraries that are used to access pi2 functionality from Python and C# environments.
* nr_stitcher Python program that can be used for non-rigid stitching or terabyte-scale volume image mosaics.

Pi2 is not intended to be a replacement for huge libraries like ITK and VTK. Instead, it is meant to be simple, performant, and easily extensible program well suitable for original author's own algorithm development and analysis tasks. It has been very useful for the author during the last years but... your mileage may vary.


## Help

The user-oriented help of pi2, including examples, can be read at [ReadTheDocs](https://pi2-docs.readthedocs.io/en/latest/).


## Binaries

Pre-built binaries for Windows and CentOS Linux can be downloaded from the [releases page](https://github.com/arttumiettinen/pi2/releases).


## Design decisions

* Compactness: Compact package such that the version of the program used to process the data can be archived along the data - at least in principle.
* Simplicity: Everything must be as simple and understandable as possible. Therefore, "C with classes" approach is chosen. Every operation is preferably a simple function call. Only a minimal number of stateful objects is used.
* No black/gray boxes: Where possible, implementations of algorithms are not hidden behind multiple levels of derived classes and template specializations. Implementation must be as readable as possible even for inexperienced developers.
* Performance: All algorithms should be as performant as possible and parallelized, but optimization efforts must be focused towards bottlenecks. OpenMP has been chosen as threading library. Code should be written such that automatic SIMD compilation is easy. For the original author this point is easier to ensure if iteration is explicit and not hidden inside iterator objects. Therefore, iteration is made using explicit loops.
* When there are multiple possibilities, the library should automatically choose the most performant algorithm.
* Algorithms are designed for 2D and 3D images unless nD implementation is simpler than implementation bound to specific image dimensions.
* Complex operations are not hidden inside operator overloads. If a calculation might take long time for a big image, it should be initiated using a function call.



## Building

When building from the source, consider using branch 'experimental'. The most recent, but possibly buggy, changes are there.

### Linux

The overall build process is as follows:
* Make sure that gcc 8.4.0 or newer is installed. Often you also need to install build-essential or corresponding package.
* Make sure that FFTW 3 library and its development packages are installed, or place its source to fftw-3.3.7-src folder and run build_again.sh.
* Make sure that libpng and libtiff libraries and their development packages are installed.
* For Python support make sure that Python 3 is installed.
* For OpenCL support make sure that you have suitable OpenCL development files installed.
* Run "make" to generate OpenCL-enabled build or "make NO_OPENCL=1" if no OpenCL should be used. The output is placed in folder bin-linux64.

Typically in an Ubuntu-like system you would run something like this:
```
sudo apt install build-essential libfftw3-dev libpng-dev libtiff-dev libjpeg-dev 
git clone https://github.com/arttumiettinen/pi2.git
cd pi2
make NO_OPENCL=1
```

You can install the executable and libraries to any standard location, but often it is better to just copy the files along with your project. This guarantees that you know which version of the program you used to generate the results.

NOTE: The default makefile compiles the programs for the processor type of the computer where the compilation is done. E.g. in heterogeneous clusters not all nodes might support the same instruction set.
In those cases you will get an 'Illegal instruction'-runtime error. To fix the problem, please determine suitable value for gcc march parameter (https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html#x86-Options) and enter that into the CXXFLAGS in the main makefile.


### Windows

* FFTW, libpng, zlib, libtiff, and libjpeg are required. It's easiest to download all the dependencies from the binary download page, see above. Either use pre-built binaries and place them to folders fftw-3.3.5-dll64, libpng-1.6.34, zlib-1.2.11, tiff-4.0.10, and jpeg-9e, or build them from sources using default Release Library x64 build settings. In particular, libtiff must be built with nmake from x64 Developer Command Prompt using command
```
[path-to-base-folder]\tiff-4.0.5> nmake /f makefile.vc
```
Before building edit nmake.opt file and change OPTFLAGS value /MD and /MDd to /MT and /MTd, for debug and release builds, respectively.
* Build everything in itl2.sln solution file, selecting either Release or Release no OpenCL configuration depending on whether you have OpenCL available. The output is placed to the x64 folder.


## License

This software is licensed under [GNU General Public License v3.0.](LICENSE.txt)

Some external parts of the software (that might be statically or dynamically linked depending on configuration) are licensed under their respective licenses. Please run `pi2 license` for most up-to-date information.

