import sys

import numpy as np

sys.path.append("x64/Release no OpenCL")

from python_scripts.pi2py2 import *
pi2 = Pi2()

def output_file(name):
    return "testoutput/" + name

write_img = pi2.newimage(ImageDataType.FLOAT32, 5, 5, 5)
arr = np.arange(5*5*5, dtype=np.float32).reshape(5,5,5)
write_img.set_data(arr)
pi2.writenn5(write_img, output_file("nn5test"))


read_img = pi2.read(output_file("nn5test"))
read_arr = read_img.get_data()

assert np.array_equal(arr, read_arr)


