Building pi2 in Linux
---------------------

The overall build process is as follows:
* Make sure that gcc 8.4.0 or newer is installed. Often you also need to install build-essential or corresponding package.
* Make sure that FFTW 3 library, libpng, libtiff, libjpeg, and their development packages are installed.
* For Python support make sure that Python 3 is installed.
* For OpenCL support make sure that you have suitable OpenCL development files installed.
* Run `make` to generate OpenCL-enabled build or `make NO_OPENCL=1` if no OpenCL is desired.

Typically in an Ubuntu-like system you would run something like this:
```
sudo apt install build-essential libfftw3-dev libpng-dev libtiff-dev libjpeg-dev 
git clone https://github.com/arttumiettinen/pi2.git
cd pi2
make NO_OPENCL=1
```

The output is placed in folder `bin-linux64`. You can install the executable and libraries to any standard location, but often it is better to just copy the files along with your project. This guarantees that you know which version of the program you used to generate the results.

**Note**
The default makefile compiles the programs for the processor architecture of the computer where the compilation is done.
E.g., in heterogeneous clusters not all nodes might support the same instruction set.
In those cases you will get an 'Illegal instruction'-runtime error. To fix the problem, please determine suitable value for gcc march parameter (https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html#x86-Options) and enter that into the `CXXFLAGS` in the main [Makefile](Makefile) or to a `makefile.local` file.
The project root folder contains some default `makefile.local.*` files for different environments. Those can be taken advantage of by selecting the correct one and renaming it to `makefile.local`.

