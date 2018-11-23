#! /bin/sh

make distclean
./configure --prefix=/afs/psi/project/miettinen_data/dev/itl2/fftw-3.3.7-linux64/ --enable-shared --enable-float --enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-threads --with-openmp --disable-fortran
# These options seem to make the results of FFTs bad...
#--enable-avx-128-fma --enable-generic-simd128 --enable-generic-simd256
make
make install


