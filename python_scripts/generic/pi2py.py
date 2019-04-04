#!/usr/bin/env python3

from ctypes import *
import os
import atexit
import string
import random
import numpy as np

class Pi2:
    """
    Wraps Pi2 library.
    Use help() and info() methods to get further usage instructions.
    """

    # Dynamic library object
    pilib = 0

    # Handle that is used to access pilib
    piobj = 0


    # Interpolation modes
    NEAREST_INTERPOLATION = "nearest"
    nearest_interpolation = NEAREST_INTERPOLATION
    LINEAR_INTERPOLATION = "linear"
    linear_interpolation = LINEAR_INTERPOLATION
    CUBIC_INTERPOLATION = "cubic"
    cubic_interpolation = CUBIC_INTERPOLATION

    # Boundary conditions
    ZERO_BOUNDARY = "zero"
    zero_boundary = ZERO_BOUNDARY
    NEAREST_BOUNDARY = "nearest"
    nearest_boundary = NEAREST_BOUNDARY

    # Neighbourhood types
    RECTANGULAR = "rectangular"
    rectangular = RECTANGULAR
    rect = RECTANGULAR
    BOX = RECTANGULAR
    box = BOX
    ELLIPSOIDAL = "ellipsoidal"
    ellipsoidal = ELLIPSOIDAL
    SPHERICAL = ELLIPSOIDAL
    spherical = SPHERICAL
    SPHERE = SPHERICAL
    sphere = SPHERICAL

    # Image data types
    UINT8 = "uint8"
    uint8 = UINT8
    UINT16 = "uint16"
    uint16 = UINT16
    UINT32 = "uint32"
    uint32 = UINT32
    UINT64 = "uint64"
    uint64 = UINT64
    FLOAT32 = "float32"
    float32 = FLOAT32



    def getdata(self, img_name):
        """
        Returns given image as a NumPy array.
        The image data is available as long as the image is not cleared from the pi2py system
        using clear command, new image is not created with the same name, and the size of the image
        is not changed.
        If distributed computing mode is enabled, this command causes the image to be
        read to local RAM. Use finishupdate function to flush possible changes to the data to
        files (only relevant in distributed computing mode).
        The image that has been read to local RAM has the same name than the
        distributed image and it can also be discarded using the clear command, but that discards
        both distributed and local versions.
        """

        w = c_longlong(0)
        h = c_longlong(0)
        d = c_longlong(0)
        dt = c_int(0)
        ptr = self.pilib.getImage(self.piobj, f"{img_name}".encode('ASCII'), byref(w), byref(h), byref(d), byref(dt))

        w = w.value
        h = h.value
        d = d.value
        dt = dt.value

        if w == 0 and h == 0 and d == 0:
            err = self.pilib.lastErrorMessage(self.piobj).decode('ASCII')
            raise RuntimeError(err)

        if dt == 0:
            raise RuntimeError("Unable to retrieve image data because pixel data type is not supported.")
        elif dt == 1:
            ptr = cast(ptr, POINTER(c_uint8))
        elif dt == 2:
            ptr = cast(ptr, POINTER(c_uint16))
        elif dt == 3:
            ptr = cast(ptr, POINTER(c_uint32))
        elif dt == 4:
            ptr = cast(ptr, POINTER(c_uint64))
        elif dt == 5:
            ptr = cast(ptr, POINTER(c_float))
        elif dt == 6:
            ptr = cast(ptr, POINTER(c_float))
            w = 2 * w
            print("Warning: complex32 image is output as float32 image with twice the width of the original.")
        else:
            raise RuntimeError("pilib returned unsupported image data type.")

        arr = np.ctypeslib.as_array(ptr, shape=(d, h, w))
        if d > 1:
            return np.moveaxis(arr, 0, 2).squeeze() # moveaxis should always return view of original data
        elif h > 1:
            return np.moveaxis(arr, 0, 1).squeeze() # moveaxis should always return view of original data
        elif w > 1:
            return arr.squeeze();

        return np.reshape(arr, (1))



    def getvalue(self, img_name):
        """
        Returns value of the first pixel in the image whose name is given as argument.
        """

        M = self.getdata(img_name)

        return float(M[0])



    def setdata(self, img_name, numpy_array):
        """
        Sets pixels and size of image in pi2 from a numpy array.
        The image will be in the same format than the numpy array.
        Not all data types supported in numpy work.
        Only 1-, 2-, and 3-dimensional arrays are supported.
        If the target image in pi2 does not exist, it is created.
        If the target image exists, the old image is discarded.
        """

        w = 1
        h = 1
        d = 1

        w = numpy_array.shape[0]
        if len(numpy_array.shape) >= 2:
            h = numpy_array.shape[1]
        if len(numpy_array.shape) >= 3:
            d = numpy_array.shape[2]
        if len(numpy_array.shape) >= 4:
            raise RuntimeError("Maximum 3-dimensional arrays can be transferred to pilib.")


        self.newimage(img_name, numpy_array.dtype, h, w, d)
        target_array = self.getdata(img_name)
        if numpy_array.size > 0:
            target_array.squeeze()[:] = numpy_array.squeeze()[:]
        self.finishupdate(img_name)


    def finishupdate(self, img_name):
        """
        In distributed computing mode, flushes local changes to image data to disk.
        Local changes can be made through objects returned by getdata method.
        """

        if not self.pilib.finishUpdate(self.piobj, f"{img_name}".encode('ASCII')):
            err = self.pilib.lastErrorMessage(self.piobj).decode('ASCII')
            raise RuntimeError(err)


    def id_generator(self, size=10, chars=string.ascii_uppercase):
        """
        Generates random string of given length, containing characters from given list.
        """

        return ''.join(random.choice(chars) for _ in range(size))



    def run_command(self, cmd_name, args):
        """
        Runs pilib command.
        """

        numpy_to_pi = {}
        arg_line = ""
        for arg in args:
            if isinstance(arg, np.ndarray):
                # Argument is numpy array. Copy it to PI system.

                # Create temporary image name
                img_name = self.id_generator()
                numpy_to_pi[img_name] = arg

                # Set data of that image to the numpy array
                # TODO: Make possibility to avoid copying the data and construct temporary image directly from numpy data pointer.
                self.setdata(img_name, arg)

                # Add temporary image name to argument list
                if len(arg_line) > 0:
                    arg_line = arg_line + ', '
                arg_line = arg_line + img_name

            else:
                # Argument is something else than numpy array. Just convert it to string.
                if len(arg_line) > 0:
                    arg_line = arg_line + ', '
                arg_line = arg_line + str(arg)

        # Old version that does not support numpy array input/output
        #arg_line = ', '.join(map(str, args))

        cmd_line = f"{cmd_name}({arg_line})"

        if not self.pilib.run(self.piobj, cmd_line.encode('ASCII')):
            err = self.pilib.lastErrorMessage(self.piobj).decode('ASCII')
            raise RuntimeError(err)

        # Get updated image data back from pi system and delete the temporary images from the pi system.
        for img_name, nparray in numpy_to_pi.items():
            data = self.getdata(img_name).squeeze()

            nparray.resize(data.shape, refcheck=False)
            if len(nparray.shape) > 0:
                nparray[:] = data[:]
            else:
                nparray[()] = data[()]

            #if data.shape == nparray.shape:
            #    nparray[:] = data[:]
            #else:
            #    nparray = data

            self.clear(img_name)



    def add_method(self, cmd_name):
        """
        Adds a pilib method to this class.
        cmd_name - Name of command to add
        """

        doc = self.pilib.help(self.piobj, f"{cmd_name}".encode('ASCII')).decode('ASCII')

        func = lambda *args: self.run_command(cmd_name, args)
        func.__doc__ = doc
        setattr(self, cmd_name, func)


    def __init__(self):

        # (No docstring as it is visible in IPython class documentation.)
        # Finds commands that are available from the pilib and adds them as methods of the class instance.


        # Create pilib instance
        if os.name == 'nt':
            # Windows
            os.environ['PATH'] = os.path.dirname(__file__) + ';' + os.environ['PATH']
            self.pilib = WinDLL("pi.dll")
        else:
            # Linux
            #os.environ['PATH'] = os.path.dirname(__file__) + ':' + os.environ['PATH']
            self.pilib = CDLL(f"{os.path.dirname(__file__)}/libpilib.so")



        self.pilib.createPI.restype = c_void_p

        self.pilib.destroyPI.argtypes = [c_void_p]

        self.pilib.run.restype = c_uint8
        self.pilib.run.argtypes = [c_void_p, c_char_p]

        self.pilib.lastErrorMessage.restype = c_char_p
        self.pilib.lastErrorMessage.argtypes = [c_void_p]

        self.pilib.lastErrorLine.restype = c_int32
        self.pilib.lastErrorLine.argtypes = [c_void_p]

        self.pilib.clearLastError.argtypes = [c_void_p]

        self.pilib.commandList.argtypes = [c_void_p]
        self.pilib.commandList.restype = c_char_p

        self.pilib.help.argtypes = [c_void_p, c_char_p]
        self.pilib.help.restype = c_char_p

        self.pilib.getImage.restype = c_void_p
        self.pilib.getImage.argtypes = [c_void_p, c_char_p, POINTER(c_int64), POINTER(c_int64), POINTER(c_int64), POINTER(c_int32)]
        self.pilib.getImage.paramflags = [1, 1, 2, 2, 2, 2]

        self.pilib.finishUpdate.restype = c_uint8
        self.pilib.finishUpdate.argtypes = [c_void_p, c_char_p]


        def cleanup(ptr):
            ptr.closepi()

        atexit.register(cleanup, self)
        self.piobj = self.pilib.createPI()

        # Get list of all commands
        commands = self.pilib.commandList(self.piobj).decode('ASCII')
        commands = commands.splitlines()

        # For each command, create corresponding method to this class
        for command in commands:
            command = command.replace('(', ',')
            command = command.replace(')', ',')
            parts = command.split(',')
            cmd_name = parts[0]

            self.add_method(cmd_name)


    def closepi(self):
        """
        Closes the PI system if it was initialized.
        """
        if self.piobj != 0:
            self.pilib.destroyPI(self.piobj)
            self.piobj = 0


    def __del__(self):
        self.closepi()
