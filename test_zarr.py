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
#shutil.rmtree("testoutput", ignore_errors=True)

arr = np.arange(5*5*5, dtype=np.float32).reshape(5,5,5)

def pi2_write(provide_chunk_shape=True):
    name = "test.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, 5, 5, 5)
    write_img.set_data(arr)
    args = write_img, output_file(name)
    if provide_chunk_shape:
        args = write_img, output_file(name), [1, 5, 5]
    pi2.writezarr(*args)

def zarrita_write():
    name = "zarrita.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.create(
        store,
        shape=arr.shape,
        dtype='int32',
        chunk_shape=arr.shape,
    )
    a[:] = arr

def pi2_read(name):
    read_img = pi2.read(output_file(name))
    read_arr = read_img.get_data()
    print(read_arr)
    return read_arr

def zarrita_read(name):
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.open(store)
    read_arr = a[:]
    print(read_arr)
    return read_arr


def test_zarrita_to_pi2():
    zarrita_write()
    read_arr = pi2_read("zarrita.zarr")
    assert np.array_equal(arr, read_arr)
    print("passed")

def test_pi2_to_zarrita():
    pi2_write()
    zarrita_write()
    read_arr = zarrita_read("test.zarr")
    assert np.array_equal(arr, read_arr)
    print("passed")

def test_pi2_no_chunk_shape_to_zarrita():
    pi2_write(provide_chunk_shape=False)
    zarrita_write()
    read_arr = zarrita_read("test.zarr")
    assert np.array_equal(arr, read_arr)
    print("passed")

def test_pi2_to_pi2():
    pi2_write()
    read_arr = pi2_read("test.zarr")
    assert np.array_equal(arr, read_arr)
    print("passed")

# passed
def test_zarrita_to_zarrita():
    zarrita_write()
    read_arr = zarrita_read("zarrita.zarr")
    assert np.array_equal(arr, read_arr)
    print("passed")

# passed
def test_pi2_to_zarrita_transpose_corrected():
    pi2_write()
    read_arr = zarrita_read("test.zarr")
    assert np.array_equal(arr, read_arr) or np.array_equal(arr, read_arr.transpose(1, 2, 0))
    print("passed")

