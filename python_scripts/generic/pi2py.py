#!/usr/bin/env python3

from ctypes import *
import os
import numpy
import atexit

class Pi2:
    """
    Wraps Pi2 library.
    Use help() and info() methods to get further usage instructions.
    """

    # Dynamic library object
    pilib = 0

    # Handle that is used to access pilib
    piobj = 0



    def run_command(self, cmd_name, args):
        """
        Runs pilib command.
        """
        arg_line = ', '.join(map(str, args))
        cmd_line = f"{cmd_name}({arg_line})"
        
        if not self.pilib.run(self.piobj, cmd_line.encode('ASCII')):
            err = self.pilib.lastErrorMessage(self.piobj).decode('ASCII')
            raise RuntimeError(err)


    def getdata(self, img_name):
        """
        Returns given image as a NumPy array.
        The image data is available as long as the image is not cleared from the pi2py system
        using clear command, new image is not created with the same name, and the size of the image
        is not changed.
        The axes in the output array are ordered (x, y, z).
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
        
        arr = numpy.ctypeslib.as_array(ptr, shape = (d, h, w))
        if (w > 1) or (h > 1) or (d > 1):
            #return numpy.transpose(arr, (1, 2, 0)) # This will copy data IF numpy does not return view
            return numpy.moveaxis(arr, 0, 2).squeeze() # moveaxis should always return view of original data
        
        return numpy.reshape(arr, (1))
       


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

        #self.pilib.lastErrorMessage.restype = c_char_p
        #self.pilib.commandList.restype = c_char_p
        #self.pilib.help.restype = c_char_p
        #self.pilib.getImage.restype = c_void_p
        #self.pilib.getImage.argtypes = [c_void_p, c_char_p, POINTER(c_longlong), POINTER(c_longlong), POINTER(c_longlong), POINTER(c_int)]
        #self.pilib.getImage.paramflags = [1, 1, 2, 2, 2, 2]

        

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
        
        
        def cleanup(ptr):
            ptr.closepi();
        
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
        self.closepi();

