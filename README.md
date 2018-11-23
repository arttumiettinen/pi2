
# itl2 and pi2

This repository contains
- itl2 project, a C++ Image analysis Template Library (version 2). It contains a small number of parallelized algorithms for processing 3D volume images e.g. from computed tomography scanners.
- pi2 program (Process Image) that exposes functinality of itl2 library for end users. It also adds the possibility to automatically distribute processing of large images to a computing cluster.
- pi2py and pi2cs libraries that are used to access pi2 functionality from Python and C# environments.
- nr_stitcher Python program that can be used for non-rigid stitching or terabyte-scale volume image mosaics.

Pi2 is not intended to be a replacement for huge libraries like ITK and VTK. Instead, it is meant to be simple, performant, and easily extensible program well suitable for original author's own algorithm development and analysis tasks. It has been very useful for the author during the last years but... your mileage may vary.


## Binaries

Pre-built binaries for Windows and CentOS Linux can be downloaded from the [releases page](https://github.com/arttumiettinen/pi2/releases/latest). The linux binaries probably work on other Linux distributions, too, as long as library versions match.


## Examples

Examples including small test datasets can be downloaded from the [releases page](https://github.com/arttumiettinen/pi2/releases/latest).


## Getting help

The program should be relatively easy to use if you are familiar with command line tools or Python. It's best to start from examples (see above), but there is also some help built-in to the pi2 program and in the Python bindings.

For raw pi2: Open command line or terminal and run `pi2` or `pi2 help` to get a list of commands.

From Python:
```
import pi2py
pi2 = pi2py.Pi2()
pi2.help()
```


## Building

### Linux

- Make sure that FFTW 3 library and its development package is installed (or place its source to fftw-3.3.7-src folder, modify paths in build_again.sh to point to correct folders and run it).
- Make sure that libpng library and its development package is installed.
- Make sure that g++ 7.3 is installed.
- For Python support make sure that Python 3 is installed.
- Run make in the project root directory. The output is placed in folder bin-linux64/release.

You can install the executable and libraries to any standard location, but often it is simpler to copy the files along with your image analysis project. This guarantees that you know which version of the program you used to generate your results.


### Windows

- FFTW, libpng and zlib are required. It's easiest to download all these dependencies from the [releases page](https://github.com/arttumiettinen/pi2/releases/latest). Just unpack the .zip file to the project root directory. Another possibility is to build the dependencies from sources using default Release Library x64 build settings.
- Build everything in itl2.sln solution file using Visual Studio. The output is placed to x64/Release folder.



## License

This software is licensed under [GNU General Public License v3.0.](LICENSE.txt)

Some external parts of the software (that might be statically or dynamically linked depending on configuration) are licensed under their respective licenses. At least libpng, zlib, and FFTW 3 are used. Some code from third parties licensed separately is also used. Please run `pi2 license` for most up-to-date information.


## Design decisions & some reasoning why things are done as they are

- Simplicity: Everything must be as simple and understandable as possible. Therefore, "C with classes" approach is mostly used. Every operation is preferably a simple function call. Only a minimal number of stateful objects is used.
- Performance: All algorithms should be as performant as possible and parallelized. OpenMP has been chosen as threading library. Code should be written such that automatic SIMD compilation is easy. For the original author this point is easier to ensure if iteration is explicit and not hidden inside iterator objects. Therefore, iteration is made using explicit loops.
- When there are multiple possibilities, the library should automatically choose the most performant algorithm.
- Implementation must be as readable as possible even for inexperienced developers.
- Algorithms are designed for 2D and 3D images unless nD implementation is simpler than implementation bound to specific image dimensions. Complexity is never increased by adding nD implementation.
- Complex operations are not hidden inside operator overloads. If a calculation might take long time for a big image, it should be initiated using a function call.


## Contact

Problems or questions? Contact arttu.miettinen@psi.ch.

