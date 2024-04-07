Building pi2 in Windows
-----------------------

**Note**
Typically, in Windows you don't need to build from source. Just download the latest version from the [download page](https://github.com/arttumiettinen/pi2/releases).

In Windows, pi2 can be built with Visual Studio (Community edition).

Before building, FFTW, libpng, zlib, libtiff, and libjpeg libraries are required.
It's easiest to download all the dependencies from the [binary download page](https://github.com/arttumiettinen/pi2/releases).
Either
* use pre-built binaries and place them to folders `fftw-3.3.5-dll64`, `libpng-1.6.34`, `zlib-1.2.11`, `tiff-4.0.10`, and `jpeg-9e`, or
* build them from sources using default Release Library x64 build settings.

In particular, libtiff must be built with `nmake` from the x64 Developer Command Prompt using command
```
[path-to-base-folder]\tiff-4.0.10> nmake /f makefile.vc
```
Before building, edit `nmake.opt` file and change `OPTFLAGS` value `/MD` and `/MDd` to `/MT` and `/MTd`, for debug and release builds, respectively.

Finally, build everything in `itl2.sln` solution file in Visual Studio, selecting either Release or Release no OpenCL configuration depending on whether you have OpenCL available.
The output is placed to the `x64` folder.
