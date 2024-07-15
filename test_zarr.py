import os
import shutil
import sys

import pytest

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
# shutil.rmtree("testoutput", ignore_errors=True)

arr = np.arange(10 * 10 * 10, dtype=np.float32).reshape(10, 10, 10)
w = 3
h = 5
d = 7
arr = arr[:w, :h, :d]


def pi2_write(chunk_shape):
    name = "test.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    write_img = pi2.newimage(pi2py2.ImageDataType.UInt32, [h, w, d])
    write_img.set_data(arr.transpose(1, 0, 2))
    pi2.writezarr(write_img, output_file(name), chunk_shape)


def zarrita_write(chunk_shape=None, codecs=None):
    if chunk_shape is None:
        chunk_shape = [w, h, 1]
    if codecs is None:
        codecs = [
            zarrita.codecs.bytes_codec("little"),
        ]
    name = "zarrita.zarr"
    shutil.rmtree(output_file(name), ignore_errors=True)
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.create(
        store,
        shape=arr.shape,
        dtype='int32',
        chunk_shape=tuple(chunk_shape),
        fill_value=42,
        codecs=codecs,
    )
    a[:] = arr


def pi2_read(name):
    read_img = pi2.read(output_file(name))
    read_arr = read_img.get_data()
    return read_arr.transpose(1, 0, 2)


def zarrita_read(name):
    store = zarrita.LocalStore(output_file(name))
    a = zarrita.Array.open(store)
    read_arr = a[:]
    return read_arr


@pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
def test_pi2_to_pi2(chunk_shape):
    pi2_write(chunk_shape)
    read_arr = pi2_read("test.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)


@pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
def test_zarrita_to_zarrita(chunk_shape):
    zarrita_write(chunk_shape)
    read_arr = zarrita_read("zarrita.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)


@pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
def test_pi2_to_zarrita(chunk_shape):
    pi2_write(chunk_shape)
    read_arr = zarrita_read("test.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)


@pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
def test_zarrita_to_pi2(chunk_shape):
    zarrita_write(chunk_shape)
    read_arr = pi2_read("zarrita.zarr")
    # write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, w, h, d)
    # write_img.set_data(read_arr)
    # pi2.writezarr(write_img, output_file("test.zarr"), chunk_shape)
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)


@pytest.mark.parametrize("order", [[0, 1, 2], [1, 0, 2], [2, 1, 0], [2, 0, 1]])
def test_read_transpose(order):
    zarrita_write(codecs=[zarrita.codecs.transpose_codec(order), zarrita.codecs.bytes_codec("little")])
    read_arr = pi2_read("zarrita.zarr")
    assert np.array_equal(arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(arr)
