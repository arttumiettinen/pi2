import sys
sys.path.append('bin-linux64/release-nocl')
import pi2py2
import numpy as np

pi2 = pi2py2.Pi2()

def output_file(name):
    return "testoutput/" + name

write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, 5, 5, 5)
arr = np.arange(5*5*5, dtype=np.float32).reshape(5,5,5)
write_img.set_data(arr)
# pi2.writenn5(write_img, output_file("nn5test"))
pi2.writezarr(write_img, output_file("test.zarr"))


#read_img = pi2.read(output_file("nn5test"))
read_img = pi2.read(output_file("test.zarr"))
read_arr = read_img.get_data()

assert np.array_equal(arr, read_arr)


