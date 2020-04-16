#! /bin/sh

make distclean

# This version works on Crick's cluster.
#./configure --prefix=$(pwd)/../fftw-3.3.7-linux64 --enable-shared --enable-float --enable-threads --with-openmp --disable-fortran

# This version works on PSI Ra cluster.
./configure --prefix=$(pwd)/../fftw-3.3.7-linux64 --enable-shared --enable-float --enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-threads --with-openmp --disable-fortran

# These options seem to make the results of FFTs bad... without warning
#--enable-avx-128-fma --enable-generic-simd128 --enable-generic-simd256

make
make install


