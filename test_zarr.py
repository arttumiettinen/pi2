import os
import shutil
import sys
sys.path.append('bin-linux64/release-nocl')
import pi2py2
import numpy as np
import zarrita

pi2 = pi2py2.Pi2()

def output_file(name):
    return "testoutput/" + name

# remove all file in testoutput
shutil.rmtree("testoutput", ignore_errors=True)

arr = np.arange(5*5*5, dtype=np.float32).reshape(5,5,5)

# pi2 zarr
write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, 5, 5, 5)
write_img.set_data(arr)
pi2.writezarr(write_img, output_file("test.zarr"))

# zarrita
store = zarrita.LocalStore(output_file("zarrita.zarr"))
a = zarrita.Array.create(
    store,
    shape=arr.shape,
    dtype='int32',
    chunk_shape=arr.shape,
)
a[:] = arr

#read_img = pi2.read(output_file("nn5test"))
#read_img = pi2.read(output_file("test.zarr"))
#read_arr = read_img.get_data()

#assert np.array_equal(arr, read_arr)


