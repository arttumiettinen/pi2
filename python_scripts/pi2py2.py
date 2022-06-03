#!/usr/bin/env python3

from ctypes import *
import os
import atexit
import string
import random
import numpy as np
from enum import Enum




class Connectivity(Enum):
    """
    Defines connectivity of pixels.
    NEAREST indicates only nearest neighbours (4 in 2D, 6 in 3D).
    ALL indicates all neighbours (8 in 2D, 26 in 3D).
    """

    NEAREST = "nearest"
    Nearest = NEAREST
    ALL = "all"
    All = ALL

    def __str__(self):
        return str(self.value)


class BoundaryCondition(Enum):
    """
    Defines boundary conditions.
    ZERO declares that values outside the image are assumed to be zero.
    NEAREST declares that values outside of the image are assumet to equal to the nearest value inside the image.
    """

    ZERO = "zero"
    Zero = ZERO
    NEAREST = "nearest"
    Nearest = NEAREST

    def __str__(self):
        return str(self.value)


class Interpolation(Enum):
    """
    Defines interpolation mode.
    NEAREST corresponds to nearest neighbour interpolation.
    LINEAR corresponds to linear interpolation.
    CUBIC corresponds to cubic interpolation.
    """

    NEAREST = "nearest"
    Nearest = NEAREST
    LINEAR = "linear"
    Linear = LINEAR
    CUBIC = "cubic"
    Cubic = CUBIC

    def __str__(self):
        return str(self.value)


class NeighbourhoodType(Enum):
    """
    Defines type of neighbourhood.
    RECTANGULAR or BOX defines neighbourhood that is rectangular.
    ELLIPSOIDAL or SPHERE defines neighbourhood that is ellipsoidal or spherical.
    """

    RECTANGULAR = "rectangular"
    Rectangular = RECTANGULAR
    BOX = RECTANGULAR
    Box = BOX
    ELLIPSOIDAL = "ellipsoidal"
    Ellipsoidal = ELLIPSOIDAL
    SPHERE = ELLIPSOIDAL
    Sphere = SPHERE

    def __str__(self):
        return str(self.value)


class ImageDataType(Enum):
    """
    Enumerates supported pixel data types.
    UINT8, UINT16, UINT32 and UINT64 correspond to 1-, 2-, 4-, and 8- byte integer value.
    FLOAT32 corresponds to 4-byte floating point value.
    COMPLEX32 corresponds to complex number consisting of two 4-byte floating point values (i.e. FLOAT32 values).
    """

    UNKNOWN = "unknown"
    Unknown = UNKNOWN
    UINT8 = "uint8"
    UInt8 = UINT8
    UINT16 = "uint16"
    UInt16 = UINT16
    UINT32 = "uint32"
    UInt32 = UINT32
    UINT64 = "uint64"
    UInt64 = UINT64
    INT8 = "int8"
    Int8 = INT8
    INT16 = "int16"
    Int16 = INT16
    INT32 = "int32"
    Int32 = INT32
    INT64 = "int64"
    Int64 = INT64
    FLOAT32 = "float32"
    Float32 = FLOAT32
    COMPLEX32 = "complex32"
    Complex32 = COMPLEX32

    def __str__(self):
        return str(self.value)



class Distributor(Enum):
    """
    Enumerates supported distributed processing modes.
    NONE corresponds to normal non-distributed processing.
    LOCAL mode is used to divide large job into smaller chunks, and to run each of the chunks sequentially on the local computer.
    This mode can be used to process images that do not fit into the RAM of the local computer.
    SLURM and LSF modes are similar to the LOCAL mode but the individual jobs are submitted to a computing cluster using Slurm or IBM LSF job management system.
    These modes are usually used when the code runs on the login node of the cluster.
    """

    NONE = ""
    No = NONE
    OFF = NONE
    Off = NONE

    LOCAL = "local"
    Local = LOCAL
    SLURM = "slurm"
    Slurm = SLURM
    LSF = "lsf"
    Lsf = LSF

    def __str__(self):
        return str(self.value)



class Direction(Enum):
    """
    Enumerates directions up or down.
    """

    UP = "up"
    Up = UP
    DOWN = "down"
    Down = DOWN
    
    def __str__(self):
        return str(self.value)


class ResliceDirection(Enum):
    """
    Enumerates directions used in reslice command.
    """

    TOP = "top"
    Top = TOP
    BOTTOM = "bottom"
    Bottom = BOTTOM
    LEFT = "left"
    Left = LEFT
    RIGHT = "right"
    Right = RIGHT

    def __str__(self):
        return str(self.value)


class EllipsoidType(Enum):
    """
    Enumerates ellipsoid types that are supported by the drawellipsoids function.
    """

    PRINCIPAL = "Principal"
    Principal = PRINCIPAL
    BOUNDING = "BoundingEllipsoid"
    Bounding = BOUNDING
    BOUNDINGSPHERE = "BoundingSphere"
    BoundingSphere = BOUNDINGSPHERE
    VOLUME = "CorrectVolume"
    Volume = VOLUME

    def __str__(self):
        return str(self.value)


class AutoThresholdMethod(Enum):
    """
    Enumerates auto-threshold methods.
    """

    OTSU = "Otsu"
    Otsu = OTSU
    HUANG = "Huang"
    Huang = HUANG
    INTERMODES = "Intermodes"
    Intermodes = INTERMODES
    ISODATA = "IsoData"
    Isodata = ISODATA
    LI = "Li"
    Li = LI
    MAXENTROPY = "MaxEntropy"
    MaxEntropy = MAXENTROPY
    MEAN = "Mean"
    Mean = MEAN
    MINERROR = "MinError"
    MinError = MINERROR
    MINIMUM = "Minimum"
    Minimum = MINIMUM
    MOMENTS = "Moments"
    Moments = MOMENTS
    PERCENTILE = "Percentile"
    Percentile = PERCENTILE
    RENYI = "Renyi"
    Renyi = RENYI
    SHANBHAG = "Shanbhag"
    Shanbhag = SHANBHAG
    TRIANGLE = "Triangle"
    Triangle = TRIANGLE
    YEN = "Yen"
    Yen = YEN
    MEDIAN = "Median"
    Median = MEDIAN
    MIDGREY = "MidGrey"
    MidGrey = MIDGREY
    NIBLACK = "Niblack"
    Niblack = NIBLACK
    PHANSALKAR = "Phansalkar"
    Phansalkar = PHANSALKAR
    SAUVOLA = "Sauvola"
    Sauvola = SAUVOLA
    BERNSEN = "Bernsen"
    Bernsen = BERNSEN


    def __str__(self):
        return str(self.value)




class Pi2Object:
    """
    Base class for objects (values and images) stored in the Pi2 system.
    """

    def __init__(self, pi2, name):
        """
        Do not use the constructor directly, instead use Pi2.newimage etc.
        """
        self.pi2 = pi2
        self.name = name


    def clear(self):
        """
        Clears this image from Pi2 system.
        After a call to clear this image cannot be used in pi2 commands.
        """
        if not self.pi2.piobj is None:
            self.pi2.clear(self.name)


    def __del__(self):
        """
        Destructor.
        Clears this image from Pi2 system.
        """
        self.clear()


    def __str__(self):
        
        return "Pi2 object: " + self.name


class Pi2Value(Pi2Object):
    """
    Represents any non-image value stored in the Pi2 system.
    """

    def as_string(self):
        """
        Retrieves the value of the object as string.
        Raises error if the value is not a string.
        """

        val = self.pi2.pilib.getString(self.pi2.piobj, self.name.encode('UTF-8'))
        if val is None:
            self.pi2.raise_last_error()
        return val.decode('UTF-8')
        





class Pi2Image(Pi2Object):
    """
    Represents image stored in the Pi2 system.
    """


    def get_raw_pointer(self):
        """
        Gets raw pointer to the image data, image dimensions, and data type as integer.
        """

        w = c_longlong(0)
        h = c_longlong(0)
        d = c_longlong(0)
        dt = c_int(0)
        ptr = self.pi2.pilib.getImage(self.pi2.piobj, self.name.encode('UTF-8'), byref(w), byref(h), byref(d), byref(dt))

        w = w.value
        h = h.value
        d = d.value
        dt = dt.value

        if w == 0 and h == 0 and d == 0:
            self.pi2.raise_last_error()

        return ptr, w, h, d, dt

    def get_info(self):
        """
        Gets size and data type of the image without reading data to memory in distributed mode.
        """

        w = c_longlong(0)
        h = c_longlong(0)
        d = c_longlong(0)
        dt = c_int(0)
        self.pi2.pilib.getImageInfo(self.pi2.piobj, self.name.encode('UTF-8'), byref(w), byref(h), byref(d), byref(dt))

        w = w.value
        h = h.value
        d = d.value
        dt = dt.value

        if w == 0 and h == 0 and d == 0:
            self.pi2.raise_last_error()

        return w, h, d, dt

    def info(self):
        """
        Gets size and data type of the image without reading data to memory in distributed mode.
        """

        return self.get_info()

    def get_data_pointer(self):
        """
        Gets pointer to the image data stored in the Pi2 system as NumPy array.
        The data is not copied so it might become unavailable when next Pi2 calls are made.
        """

        ptr, w, h, d, dt = self.get_raw_pointer()

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
        elif dt == 7:
            ptr = cast(ptr, POINTER(c_int8))
        elif dt == 8:
            ptr = cast(ptr, POINTER(c_int16))
        elif dt == 9:
            ptr = cast(ptr, POINTER(c_int32))
        elif dt == 10:
            ptr = cast(ptr, POINTER(c_int64))
        else:
            raise RuntimeError("pilib returned unsupported image data type.")

        arr = np.ctypeslib.as_array(ptr, shape=(d, h, w))
        if d > 1:
            return np.moveaxis(arr, 0, 2).squeeze() # moveaxis should always return view of original data
        elif h > 1:
            return np.moveaxis(arr, 0, 1).squeeze() # moveaxis should always return view of original data
        elif w > 1:
            return arr.squeeze()

        return np.reshape(arr, (1))


    def flush_pointer(self):
        """
        In distributed computing mode, flushes local changes to image data to disk.
        Local changes can be made through objects returned by get_data_pointer method.
        """

        if not self.pi2.pilib.finishUpdate(self.pi2.piobj, self.name.encode('UTF-8')):
            self.pi2.raise_last_error()


    def get_data(self):
        """
        Gets a copy of the pixel data of this image as a NumPy array.
        """

        return np.copy(self.get_data_pointer())



    def get_value(self):
        """
        Returns the value of the first pixel in the image.
        """

        M = self.get_data_pointer()
        return float(M[0])


    def set_data(self, numpy_array):
        """
        Sets pixels, size and pixel data type of this image from a NumPy array.
        The image will be in the same format than the NumPy array.
        Not all data types supported in NumPy work.
        Only 1-, 2-, and 3-dimensional NumPy arrays are supported.
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
            raise RuntimeError("Maximum 3-dimensional arrays can be transferred to pi2.")

        dtype = numpy_array.dtype
        if numpy_array.dtype == np.float64:
            dtype = ImageDataType.FLOAT32

        self.pi2.run_script(f"newimage({self.name}, {dtype}, {h}, {w}, {d})")
        target_array = self.get_data_pointer()
        if numpy_array.size > 0:
            target_array.squeeze()[:] = numpy_array.squeeze()[:]

        self.flush_pointer()


    def get_dimensions(self):
        """
        Gets a NumPy array describing dimensions of this image.
        """

        w, h, d, dt = self.get_info()
        return np.array([w, h, d])

    def get_width(self):
        """
        Gets width of this image.
        """
        return self.get_dimensions()[0]

    def get_height(self):
        """
        Gets height of this image.
        """
        return self.get_dimensions()[1]

    def get_depth(self):
        """
        Gets depth of this image.
        """
        return self.get_dimensions()[2]

    def dimensions(self):
        """
        Gets a NumPy array describing dimensions of this image.
        """

        w, h, d, dt = self.get_info()
        return np.array([w, h, d])

    def width(self):
        """
        Gets width of this image.
        """
        return self.get_dimensions()[0]

    def height(self):
        """
        Gets height of this image.
        """
        return self.get_dimensions()[1]

    def depth(self):
        """
        Gets depth of this image.
        """
        return self.get_dimensions()[2]

    def get_data_type(self):
        """
        Gets data type of pixels as one of values of ImageDataType enumeration.
        """

        w, h, d, dt = self.get_info()

        if dt == 0:
            raise RuntimeError("Unable to retrieve image data because pixel data type is not supported.")
        elif dt == 1:
            return ImageDataType.UINT8
        elif dt == 2:
            return ImageDataType.UINT16
        elif dt == 3:
            return ImageDataType.UINT32
        elif dt == 4:
            return ImageDataType.UINT64
        elif dt == 5:
            return ImageDataType.FLOAT32
        elif dt == 6:
            return ImageDataType.COMPLEX32
        elif dt == 7:
            return ImageDataType.INT8
        elif dt == 8:
            return ImageDataType.INT16
        elif dt == 9:
            return ImageDataType.INT32
        elif dt == 10:
            return ImageDataType.INT64
        else:
            raise RuntimeError("pilib returned unsupported image data type.")

    def data_type(self):
        """
        Gets data type of pixels as one of values of ImageDataType enumeration.
        """
        return self.get_data_type()

    def __str__(self):
        """
        Converts this image to a string.
        """

        return str(self.dimensions()) + ", " + str(self.data_type())





class Pi2:
    """
    Wraps Pi2 library.
    Use help() and info() methods to get further usage instructions.
    """


    def add_method(self, cmd_name):
        """
        Adds a pilib command to this class as a member function.
        cmd_name is the Pi2 name of the command to add
        """

        if hasattr(self, cmd_name):
            # Do not add if the class already has an attribute with the same name.
            # This means that there is a specialized implementation for this functionality.
            return

        doc = self.pilib.help(self.piobj, f"{cmd_name}".encode('UTF-8')).decode('UTF-8')

        func = lambda *args: self.run_command(cmd_name, args)
        func.__doc__ = doc
        setattr(self, cmd_name, func)


    def __init__(self, library_path = ''):
        """
        Wraps Pi2 library.
        Pass path to pi.dll or libpi.so as an argument library_path, if the file is not located in
        the same directory than pi2py2.py.
        """

        # NOTE: Docstring is visible in IPython.
        # This method finds commands that are available from the pilib and adds them as methods of the class instance.

        self.piobj = None

        # Default path to folder containing pi.dll or libpi.so is the folder of this file.
        if library_path == '':
            library_path = os.path.dirname(__file__)

        # Create pilib instance
        if os.name == 'nt':
            # Windows
            os.environ['PATH'] = library_path + ';' + os.environ['PATH']
            self.pilib = WinDLL(os.path.join(library_path, "pi.dll"))
        else:
            # Linux
            #os.environ['PATH'] = library_path + ':' + os.environ['PATH']
            self.pilib = CDLL(os.path.join(library_path, "libpi.so"))

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

        self.pilib.getImageInfo.argtypes = [c_void_p, c_char_p, POINTER(c_int64), POINTER(c_int64), POINTER(c_int64), POINTER(c_int32)]
        self.pilib.getImageInfo.paramflags = [1, 1, 2, 2, 2, 2]

        self.pilib.finishUpdate.restype = c_uint8
        self.pilib.finishUpdate.argtypes = [c_void_p, c_char_p]

        self.pilib.getString.restype = c_char_p
        self.pilib.getString.argtypes = [c_void_p, c_char_p]


        def cleanup(ptr):
            if isinstance(ptr, Pi2):
                ptr.closepi()

        atexit.register(cleanup, self)
        self.piobj = self.pilib.createPI()

        # Get list of all commands
        commands = self.pilib.commandList(self.piobj).decode('UTF-8')
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
        if not self.piobj is None:
            self.pilib.destroyPI(self.piobj)
            self.piobj = None


    def __del__(self):
        self.closepi()


    def raise_last_error(self):
        """
        Raises an RuntimeError with message corresponding to the last error in the pi2 system.
        """

        err = self.pilib.lastErrorMessage(self.piobj).decode('UTF-8')
        raise RuntimeError(err)



    def run_script(self, script):
        """
        Runs any pi2 script given as string.
        """

        if not self.pilib.run(self.piobj, script.encode('UTF-8')):
            self.raise_last_error()


    def make_safe(self, argument):
        """
        Escapes bad characters from a string that is to be an argument for run_script command.
        """

        argument = str(argument)
        argument = argument.replace('\n', '\\n')
        return argument


    def run_command(self, cmd_name, args):
        """
        Runs pilib command, given its name and arguments.
        Arguments can be any objects that can be converted to string.
        NumPy arrays are copied to pilib as images, but modifications are not copied back as
        the shape of the NumPy arrays cannot necessarily be changed accordingly to image shape in Pi2.
        Therefore you should use Pi2Image objects to store output data.
        """

        temp_images = []
        arg_line = ""
        for arg in args:
            arg_as_string = ""

            if isinstance(arg, np.ndarray):
                # Argument is numpy array. Copy it to PI system.

                # Create temporary image name
                temp_image = self.newimage()
                temp_images.append(temp_image)

                # Set data of that image to the numpy array
                temp_image.set_data(arg)
                
                arg_as_string = temp_image.name

            elif isinstance(arg, Pi2Object):
                arg_as_string = arg.name

            else:
                # Argument is something else... Just convert it to string.
                arg_as_string = str(arg)

            # Remove newlines
            arg_as_string = self.make_safe(arg_as_string)
            
            # Add quotes
            arg_as_string = f"'{arg_as_string}'"

            if len(arg_line) > 0:
                arg_line = arg_line + ', '
            arg_line = arg_line + arg_as_string

        cmd_line = f"{cmd_name}({arg_line})"

        self.run_script(cmd_line)

        # Temporary images are automatically cleared when they go out of scope.




    ############ Utility functions ############


    def id_generator(self, prefix = "", size=10, chars=string.ascii_uppercase):
        """
        Generates random string of given length, containing characters from given list.
        """

        return prefix + ''.join(random.choice(chars) for _ in range(size))


    def generate_image_name(self):
        return self.id_generator("image_")

    def generate_value_name(self):
        return self.id_generator("value_")


    ############ Special Python interoperability functions ############


    def get_column_index(self, column_name, analyzers):
        """
        Get the index of a column from a particle analysis result table as a NumPy array.
        :param column_name: Column name without unit. An exception is thrown if the column name is not found.
        :param analyzers: Analyzer list string (previously passed to analyzeparticles function).
        :return: column index as integer, and unit as a string, if available; empty string if no unit has been specified in the header.
        """
        unit = ""

        headers = self.newstring()
        self.headers(analyzers, headers)
        headers = headers.as_string().split(', ')

        index = -1
        for i in range(0, len(headers)):

            hdr = headers[i]
            if hdr.lower().startswith(column_name.lower()):
                remaining = hdr[len(column_name):].strip()
                if not remaining:
                    index = i
                    break
                if remaining.startswith('[') and remaining.endswith(']'):
                    index = i
                    unit = remaining[1:-1]
                    break

        if index < 0:
            raise SystemError(f"Column {column_name} not found.")

        return index, unit


    def get_column(self, column_name, particle_analysis_result, analyzers):
        """
        Get a column from a particle analysis result table as a NumPy array.
        :param column_name: Column name without unit.
        :param particle_analysis_result: Particle analysis result image (previously passed to analyzeparticles function).
        :param analyzers: Analyzer list string (previously passed to analyzeparticles function).
        :return: column data as a Numpy array, and unit as a string, if available; empty string if no unit has been specified in the header.
        """

        index, unit = self.get_column_index(column_name, analyzers)
        data = particle_analysis_result.get_data_pointer()
        column = data[:, index]
        return column, unit



    ############ Implementations of functions that require special arrangements ############


    def newimage(self, data_type = ImageDataType.UINT8, width = 1, height = 1, depth = 1):
        """
        Creates new image to the Pi2 system and returns it as a Pi2Image object.
        """

        image_name = self.generate_image_name()
        if np.isscalar(width):
            self.run_script(f"newimage({image_name}, {data_type}, {width}, {height}, {depth})")
        else:
            self.run_script(f"newimage({image_name}, {data_type}, {width})")

        return Pi2Image(self, image_name)

    def newstring(self, value = ""):
        """
        Creates new string object.
        """

        name = self.generate_value_name()
        value = self.make_safe(value.replace("\"", "\\\""))
        self.run_script(f"newvalue({name}, \"string\", \"{value}\")")
        return Pi2Value(self, name)


    def newlike(self, template_image, data_type = ImageDataType.UNKNOWN, width = 0, height = 0, depth = 0):
        """
        Creates new image to the Pi2 system and returns it as a Pi2Image object.
        Copies image data type and other properties from the template image if they are not provided.
        """

        image_name = self.generate_image_name()
        self.run_script(f"newlike({image_name}, {template_image.name}, {data_type}, {width}, {height}, {depth})")
        
        return Pi2Image(self, image_name)


    def isimagefile(self, filename):
        """
        Checks if the given file is a readable image.
        """
        
        temp_image = self.newimage()
        filename = self.make_safe(filename)
        self.run_script(f"isimagefile({filename}, {temp_image.name})")

        return temp_image.get_value() != 0
         


    def read(self, filename, data_type_override = "", target_image = None):
        """
        Read image file to given image or to a new image.
        Returns the image object.
        Parameter data_type_override can be used to specify pixel data type in the file
        if the default is wrong. This is useful for some .raw files containing e.g. signed integer
        pixels.
        Parameter target_image can be used to load the image to an existing image instead of
        creating a new one.
        """

        if target_image == None:
            image_name = self.generate_image_name()
            target_image = Pi2Image(self, image_name)
        filename = self.make_safe(filename)
        self.run_script(f"read({target_image.name}, \"{filename}\", {data_type_override})")
        return target_image
        


    def mapraw(self, filename, data_type = ImageDataType.UNKNOWN, width = 0, height = 0, depth = 0):
        """
        Create a disk-mapped image from existing or newly created .raw file.
        Returns Pi2Image object.
        """

        image_name = self.generate_image_name()
        filename = self.make_safe(filename)
        self.run_script(f"mapraw({image_name}, {filename}, {data_type}, [{width}, {height}, {depth}])")
        return Pi2Image(self, image_name)



