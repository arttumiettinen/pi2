using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
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

        /**
        Run commands.
        @param pi Handle returned by create method.
        @param commands String containig one or more commands that should be executed.
        @return True if command execution was successful; false otherwise.
        */
        [DllImport("pi", EntryPoint = "run")]
        public static extern bool Run(IntPtr pi, string commands);

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
        public static extern IntPtr commandList(IntPtr pi);

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
        public static extern IntPtr GetImage(IntPtr pi, string imgName, out Int64 width, out Int64 height, out Int64 depth, out Int32 dataType);
    }
}
