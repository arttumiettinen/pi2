using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
// Naming Styles warning is disabled as lower/uppercase is used to separate function names in dynamic library and in wrapper class.
#pragma warning disable IDE1006 

    /// <summary>
    /// Enumerates supported pixel data types.
    /// </summary>
    public enum ImageDataType
    { 
        /**
		Unknown pixel data type.
		*/
		Unknown = 0,
		/**
		Unsigned 8-bit integer pixel data type.
		*/
		UInt8 = 1,
		/**
		Unsigned 16-bit integer pixel data type.
		*/
		UInt16 = 2,
		/**
		Unsigned 32-bit integer pixel data type.
		*/
		UInt32 = 3,
		/**
		Unsigned 64-bit integer pixel data type.
		*/
		UInt64 = 4,
		/**
		32-bit floating point pixel data type.
		*/
		Float32 = 5,
		/**
		Complex value consisting of two 32-bit floating point values.
		*/
		Complex32 = 6,
        /**
        Signed 8-bit integer pixel data type.
        */
        Int8 = 7,
        /**
		Signed 16-bit integer pixel data type.
		*/
        Int16 = 8,
        /**
		Signed 32-bit integer pixel data type.
		*/
        Int32 = 9,
        /**
		Signed 64-bit integer pixel data type.
		*/
        Int64 = 10
    }

    /**
     * Low-level interface to pi.dll.
     */
    public static class PiLib
    {

        /**
	    Create the pi library object.
	    @return Handle that must be passed to all other methods.
	    */
        [DllImport("pi", EntryPoint = "createPI")]
        public static extern IntPtr CreatePI();

        /**
        Close the pi library object.
        */
        [DllImport("pi", EntryPoint = "destroyPI")]
        public static extern void DestroyPI(IntPtr pi);

        
        [DllImport("pi", EntryPoint = "run")]
        private static extern bool run(IntPtr pi, string commands);

        /**
        Run commands.
        @param pi Handle returned by create method.
        @param commands String containig one or more commands that should be executed.
        @return True if command execution was successful; false otherwise.
        */
        public static bool Run(IntPtr pi, string commands)
        {
//#if DEBUG
//            Console.WriteLine("Running command: " + commands);
//#endif
            return run(pi, commands);
        }

        /**
        Run commands, throw ArgumentException if error occurs.
        */
        public static void RunAndCheck(IntPtr pi, string commands)
        {
            if(!Run(pi, commands))
            {
                throw new ArgumentException(LastErrorMessage(pi));
            }
        }
        
        /**
        Get error message identifying the last error that has occured.
        */
        [DllImport("pi")]
        private static extern IntPtr lastErrorMessage(IntPtr pi);

        /**
        Get error message identifying the last error that has occured.
        */
        public static string LastErrorMessage(IntPtr pi)
        {
            IntPtr ptr = lastErrorMessage(pi);
            return Marshal.PtrToStringAnsi(ptr);
        }

        /**
        Get line of script code that caused the last error.
        */
        [DllImport("pi", EntryPoint = "lastErrorLine")]
        public static extern int LastErrorLine(IntPtr pi);

        /**
        Clear last error.
        */
        [DllImport("pi", EntryPoint = "clearLastError")]
        public static extern void ClearLastError(IntPtr pi);

        /**
        Gets list of commands available
        */
        [DllImport("pi")]
        private static extern IntPtr commandList(IntPtr pi);

        /**
        Get error message identifying the last error that has occured.
        */
        public static string CommandList(IntPtr pi)
        {
            IntPtr ptr = commandList(pi);
            return Marshal.PtrToStringAnsi(ptr);
        }

        /**
	    Gets pointer to data storing the given image.
	    Stores the size of the image into the values pointed by the three last arguments.
	    If an error occurs, returns zero and sets width, height and depth to zero, and sets dataType to Unknown (zero).
	    @param dataType The system sets this int to 1 to signify uint8 image, 2 for uint16 image, 3 for float32 image and 4 for complex32 image.
	    */
        [DllImport("pi", EntryPoint = "getImage")]
        public static extern IntPtr GetImage(IntPtr pi, string imgName, out Int64 width, out Int64 height, out Int64 depth, out ImageDataType dataType);

        /**
        Gets a string value from the Pi2 system.
        */
        [DllImport("pi", EntryPoint = "getString")]
        private static extern IntPtr getString(IntPtr pi, string name);

        /**
        Gets a string value from the Pi2 system.
        */
        public static string GetString(IntPtr pi, string name)
        {
            IntPtr ptr = getString(pi, name);
            return Marshal.PtrToStringAnsi(ptr);
        }
    }


#pragma warning restore IDE1006 // Naming Styles
}
