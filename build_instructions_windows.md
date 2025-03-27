Building pi2 in Windows
-----------------------

**Note**
Typically, in Windows you don't need to build from source. Just download the latest version from the [download page](https://github.com/arttumiettinen/pi2/releases).

In Windows, pi2 can be built with Visual Studio (Community edition). Before building, FFTW, libpng, libtiff, libjpeg, and libblosc libraries are required.

It's easiest to download all the dependencies from the [binary download page](https://github.com/arttumiettinen/pi2/releases).
Either
* use pre-built binaries and place them to the `deps` subfolder, or
* build them from sources using default build settings (release, static library, 64-bit). Some notes below.

After the dependencies are sorted out, build everything in `itl2.sln` solution file in Visual Studio, selecting either Release or Release no OpenCL configuration depending on whether you have OpenCL available.
The output is placed to the `x64` folder.




libpng
------

Build libpng from the source downloaded from SourceForge, use the `projects/vstudio/vstudio.sln` file and Release Library build settings.
In some cases it appears to be necessary to set the Target Platform to x64 manually using the Configuration Manager.

The output will be in
projects\vstudio\x64\Release Library




libjpeg
-------

The jpeg library is from https://www.ijg.org/.

To build, Open Developer Command Prompt, go to jpeg library directory, and run
nmake nodebug=1 /f makefile.vs setup-v16

Open jpeg.sln
Set platform to x64
Set Runtime library to "Multi Threaded"
Build all

The output will be in
Release\x64\jpeg.lib

Optionally, to build test apps,
Open apps.sln
Set platform to x64
Build all




libtiff
-------

Libtiff is from https://download.osgeo.org/libtiff/.

To build, Open Developer Command Prompt, go to libtiff library directory, and run
cmake -DBUILD_SHARED_LIBS=OFF .
Open tiff.sln.
Go to project tiff.
Select Release build.
Change Runtime Library to Multi-threaded.
Build.

The output will be in
libtiff\Release\tiff.lib


blosc
-----

c-blosc source is from https://github.com/Blosc/c-blosc.

To build, Open Developer Command Prompt, go to c-blosc library directory, and run
cmake .
Open blosc.sln
Go to project blosc_static
Select Release build.
Change Runtime Library to Multi-threaded.
Build.

The output will be in
blosc\Release


fftw
----

FFTW dependency is the default Windows build from https://www.fftw.org/install/windows.html

Don't forget to make import libraries using commands
lib /def:libfftw3-3.def
lib /def:libfftw3f-3.def
lib /def:libfftw3l-3.def



