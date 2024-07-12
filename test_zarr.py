import os
import shutil
import sys
print(os.getcwd())
if os.getcwd().endswith("/venv/bin"):
    os.chdir("../..")
    print(os.getcwd())
sys.path.append('bin-linux64/release-nocl')
import pi2py2
import numpy as np
import zarrita

pi2 = pi2py2.Pi2()

def output_file(name):
    return "testoutput/" + name

# remove all file in testoutput
#shutil.rmtree("testoutput", ignore_errors=True)

arr = np.arange(10*10*10, dtype=np.float32).reshape(10,10,10)
w = 2
h = 3
d = 4
arr = arr[:w, :h, :d]
chunk_shape = [w,h,d]
#chunk_shape = [1,1,1]
def pi2_write():
    name = "test.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, [h, w, d])
    write_img.set_data(arr.transpose(1, 0, 2))
    pi2.writezarr(write_img, output_file(name), chunk_shape)

def zarrita_write():
    name = "zarrita.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.create(
        store,
        shape=arr.shape,
        dtype='int32',
        chunk_shape=tuple(chunk_shape),
        fill_value=42,
        codecs=[
            zarrita.codecs.bytes_codec("little"),
        ],
    )
    a[:] = arr

def pi2_read(name):
    read_img = pi2.read(output_file(name))
    read_arr = read_img.get_data()
    return read_arr

def zarrita_read(name):
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.open(store)
    read_arr = a[:]
    return read_arr





#todo test with no chunk_shape

# passed
def test_pi2_to_pi2():
    pi2_write()
    read_arr = pi2_read("test.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
    print("passed")

# passed
def test_zarrita_to_zarrita():
    zarrita_write()
    read_arr = zarrita_read("zarrita.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
    print("passed")

# # passed with chunk_shape = shape
# def test_pi2_to_zarrita_transpose_corrected():
#     pi2_write()
#     read_arr = zarrita_read("test.zarr")
#     assert np.array_equal(arr, read_arr) or np.array_equal(arr, read_arr.transpose(1, 0, 2)), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
#     print("passed")

# works for chunk_shape = [1,1,1]
def test_pi2_to_zarrita():
    pi2_write()
    zarrita_write()
    read_arr = zarrita_read("test.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
    print("passed")

def test_zarrita_to_pi2():
    zarrita_write()
    read_arr = pi2_read("zarrita.zarr")
    # write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, w, h, d)
    # write_img.set_data(read_arr)
    # pi2.writezarr(write_img, output_file("test.zarr"), chunk_shape)
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
    print("passed")

test_pi2_to_zarrita()