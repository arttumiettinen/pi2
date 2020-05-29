using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    /// <summary>
    /// Represents image stored in the Pi2 system.
    /// </summary>
    public class Pi2Image : IDisposable
    {
        /// <summary>
        /// Gets the name of this image in the Pi system.
        /// </summary>
        public string ImageName
        {
            get;
            private set;
        }

        /// <summary>
        /// Gets the Pi object that owns this image.
        /// </summary>
        public Pi2 Pi
        {
            get;
            private set;
        }

        /// <summary>
        /// Constructor.
        /// Images should be created with Pi2.newimage function, otherwise the image does not exist in the Pi2 system.
        /// </summary>
        /// <param name="pi"></param>
        /// <param name="imageName"></param>
        public Pi2Image(Pi2 pi, string imageName)
        {
            Pi = pi;
            ImageName = imageName;
        }

        #region Dispose pattern

        private bool disposed = false;

        /// <summary>
        /// Disposes resources.
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Disposes resources.
        /// </summary>
        /// <param name="disposing"></param>
        protected virtual void Dispose(bool disposing)
        {
            if (disposed)
                return;

            if (disposing)
            {
                // Free any other managed objects here.
            }

            if (ImageName != String.Empty)
            {
                // Delete image from Pi ssytem
                PiLib.Run(Pi.Handle, $"clear({ImageName})");
                ImageName = String.Empty;
            }

            disposed = true;
        }

        // This is not needed as Pi system will clear images on exit anyway.
        //~Pi2Image()
        //{
        //    Dispose(false);
        //}

        #endregion


        #region Data getters

        /// <summary>
        /// Get pointer to data of this image.
        /// The pointer is valid until this image object is modified in Pi2.
        /// </summary>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="depth"></param>
        /// <param name="dataType"></param>
        /// <returns></returns>
        public IntPtr GetData(out Int64 width, out Int64 height, out Int64 depth, out ImageDataType dataType)
        {
            IntPtr data = PiLib.GetImage(Pi.Handle, ImageName, out width, out height, out depth, out dataType);
            if (data == IntPtr.Zero || dataType == ImageDataType.Unknown)
                throw new InvalidOperationException("Image " + ImageName + " is inaccessible because it has been deleted from the Pi system.");
            return data;
        }

        /// <summary>
        /// Gets the value of a pixel in the image.
        /// </summary>
        /// <returns></returns>
        public float GetValue(int x = 0, int y = 0, int z = 0)
        {
            Int64 w, h, d;
            ImageDataType dt;
            IntPtr data = GetData(out w, out h, out d, out dt);

            if (x < 0 || x >= w ||
                y < 0 || y >= h ||
                z < 0 || z >= d)
                throw new ArgumentOutOfRangeException("xyz", "Coordinates are outside of image region.");

            unsafe
            {
                switch (dt)
                {
                    case ImageDataType.UInt8: return *((byte*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.UInt16: return *((UInt16*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.UInt32: return *((UInt32*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.UInt64: return *((UInt64*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Float32: return *((float*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Complex32: return *((float*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Int8: return *((sbyte*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Int16: return *((Int16*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Int32: return *((Int32*)data.ToPointer() + z * w * h + y * w + x);
                    case ImageDataType.Int64: return *((Int64*)data.ToPointer() + z * w * h + y * w + x);
                    default:
                        throw new NotImplementedException("Support for image data type " + dt + " is not implemented in Pi2Image.GetValue.");
                }
            }
        }

        #endregion

        #region Dimensions getters

        /// <summary>
        /// Gets dimensions of this image.
        /// </summary>
        /// <returns></returns>
        public Tuple<int, int, int> GetDimensions()
        {
            Int64 w, h, d;
            ImageDataType dt;
            // TODO: Don't call GetData here, call GetInfo as that does not read distributed images to RAM.
            GetData(out w, out h, out d, out dt);
            return Tuple.Create((int)w, (int)h, (int)d);
        }

        /// <summary>
        /// Gets width of the image.
        /// </summary>
        public int Width
        {
            get
            {
                return GetDimensions().Item1;
            }
        }

        /// <summary>
        /// Gets height of the image.
        /// </summary>
        public int Height
        {
            get
            {
                return GetDimensions().Item2;
            }
        }

        /// <summary>
        /// Gets depth of the image.
        /// </summary>
        public int Depth
        {
            get
            {
                return GetDimensions().Item3;
            }
        }

        #endregion

        /// <summary>
        /// Gets pixel data type of this image.
        /// </summary>
        /// <returns></returns>
        public ImageDataType GetDataType()
        {
            Int64 w, h, d;
            ImageDataType dt;
            GetData(out w, out h, out d, out dt);
            return dt;
        }
    }

    /// <summary>
    /// Encapsulates Pi2 image analysis system.
    /// </summary>
    public class Pi2 : IDisposable
    {
        /// <summary>
        /// Gets filter string that can be used in open dialogs to show file types supported by Pi2.
        /// </summary>
        public string SupportedOpenFileFilter
        {
            get
            {
                return "Image files (*.raw, *.tif, *.vol)|*.raw;*.tif;*.tiff;*.vol";
            }
        }

        /// <summary>
        /// PI system handle.
        /// </summary>
        public IntPtr Handle
        {
            get;
            private set;
        }

        /// <summary>
        /// Indicates if this object owns the Handle to the Pi2 system.
        /// </summary>
        public bool OwnsHandle
        {
            get;
            private set;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        public Pi2()
        {
            Handle = PiLib.CreatePI();
            OwnsHandle = true;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="handle">Handle to existing Pi2 object.</param>
        public Pi2(IntPtr handle)
        {
            Handle = handle;
            OwnsHandle = false;
        }

        #region Dispose pattern

        private bool disposed = false;

        /// <summary>
        /// Disposes resources.
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        /// Disposes resources.
        /// </summary>
        /// <param name="disposing"></param>
        protected virtual void Dispose(bool disposing)
        {
            if (disposed)
                return;

            if (disposing)
            {
                // Free any other managed objects here.
            }

            if(OwnsHandle && Handle != IntPtr.Zero)
            {
                PiLib.DestroyPI(Handle);
                Handle = IntPtr.Zero;
                OwnsHandle = false;
            }

            disposed = true;
        }

        /// <summary>
        /// Finalizer.
        /// </summary>
        ~Pi2()
        {
            Dispose(false);
        }

        #endregion

        #region Utilities

        private Random random = new Random();

        /// <summary>
        /// Generates random string.
        /// </summary>
        /// <returns></returns>
        private string RandomString()
        {
            var chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
            var stringChars = new char[8];

            for (int i = 0; i < stringChars.Length; i++)
            {
                stringChars[i] = chars[random.Next(chars.Length)];
            }

            var finalString = new String(stringChars);

            return finalString;
        }

        #endregion

        /// <summary>
        /// Creates new Pi2 image object.
        /// </summary>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="depth"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        public Pi2Image NewImage(ImageDataType dt, int width = 1, int height = 1, int depth = 1)
        {
            string imageName = "image_" + RandomString();

            PiLib.RunAndCheck(Handle, $"newimage({imageName}, {dt}, {width}, {height}, {depth})");
            return new Pi2Image(this, imageName);
        }

        /// <summary>
        /// Makes sure that the size of the image equals the given dimensions.
        /// Reallocates the image in case of size difference.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="depth"></param>
        public void EnsureSize(Pi2Image img, int width = 1, int height = 1, int depth = 1)
        {
            PiLib.RunAndCheck(Handle, $"ensuresize({img.ImageName}, {width}, {height}, {depth})");
        }

        #region I/O commands

        /// <summary>
        /// Tests if the given file name or prefix can be used successfully in Read command.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public bool IsImageFile(string filename)
        {
            using (Pi2Image result = NewImage(ImageDataType.UInt8))
            {
                PiLib.RunAndCheck(Handle, $"isimagefile({filename}, {result.ImageName})");
                return result.GetValue() != 0;
            }
        }

        /// <summary>
        /// Reads dimensions and data type of given image file.
        /// </summary>
        /// <param name="filename">Name of file to examine.</param>
        /// <param name="dimensions">Dimension of the image. Returns zero dimensions if the image file could not be read.</param>
        /// <param name="dataType">Data type of the image. Returns Unknown if the image could not be read.</param>
        public void GetImageInfo(string filename, out Vec3 dimensions, out ImageDataType dataType)
        {
            dimensions = new Vec3();
            dataType = ImageDataType.Unknown;

            using (Pi2Image result = NewImage(ImageDataType.UInt32, 4))
            {
                PiLib.RunAndCheck(Handle, $"fileinfo({filename}, {result.ImageName})");
                dimensions.X = result.GetValue(0);
                dimensions.Y = result.GetValue(1);
                dimensions.Z = result.GetValue(2);
                dataType = (ImageDataType)result.GetValue(3);
            }
        }

        /// <summary>
        /// Read image file to given image.
        /// </summary>
        /// <param name="target"></param>
        /// <param name="filename"></param>
        public void Read(Pi2Image target, string filename)
        {
            PiLib.RunAndCheck(Handle, $"read({target.ImageName}, {filename})");
        }

        /// <summary>
        /// Read image file to a new image.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public Pi2Image Read(string filename)
        {
            string imageName = "image_" + RandomString();
            PiLib.RunAndCheck(Handle, $"read({imageName}, {filename})");
            return new Pi2Image(this, imageName);
        }

        /// <summary>
        /// Read a block of image file to a new image.
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="start">Coordinates of the first pixel to read.</param>
        /// <param name="size">Size of the block to read.</param>
        /// <param name="dataType">Data type of the image. Used only if reading .raw images. Leave empty to guess data type based on file size.</param>
        /// <returns></returns>
        public Pi2Image ReadBlock(string filename, Vec3 start, Vec3 size, ImageDataType dataType = ImageDataType.Unknown)
        {
            string imageName = "image_" + RandomString();
            PiLib.RunAndCheck(Handle, $"readblock({imageName}, {filename}, {(int)start.X}, {(int)start.Y}, {(int)start.Z}, {(int)size.X}, {(int)size.Y}, {(int)size.Z}, {dataType})");
            return new Pi2Image(this, imageName);
        }

        /// <summary>
        /// Create a disk-mapped image from existing .raw file or create a new .raw file if corresponding file does not exist.
        /// If data type and dimensions are not specified, the image file name must contain dimensions of the image, i.e. it must be in format corresponding to image_name_123x456x789.raw.
        /// </summary>
        /// <param name="filename">Name of image file to map.</param>
        /// <param name="dataType">Data type of the image. Specify empty value to infer data type from image dimensions</param>
        /// <param name="width">Width of the image. Omit width, height and depth to infer dimensions from file name.</param>
        /// <param name="height">Height of the image. Omit width, height and depth to infer dimensions from file name.</param>
        /// <param name="depth">Depth of the image. Omit width, height and depth to infer dimensions from file name.</param>
        /// <returns>The mapped image. NOTE: Changes made to the image are IMMEDIATELY reflected on disk.</returns>
        public Pi2Image MapRaw(string filename, ImageDataType dataType = ImageDataType.Unknown, int width = 0, int height = 0, int depth = 0)
        {
            string imageName = "image_" + RandomString();
            Pi2Image img = new Pi2Image(this, imageName);
            MapRaw(img, filename, dataType, width, height, depth);
            return img;
        }

        /// <summary>
        /// Create a disk-mapped image from existing .raw file or create a new .raw file if corresponding file does not exist.
        /// If data type and dimensions are not specified, the image file name must contain dimensions of the image, i.e. it must be in format corresponding to image_name_123x456x789.raw.
        /// </summary>
        /// <param name="image">The mapped image is placed into this image object. NOTE: Changes made to the image are IMMEDIATELY reflected on disk.</param>
        /// <param name="filename">Name of image file to map.</param>
        /// <param name="dataType">Data type of the image. Specify empty value to infer data type from image dimensions</param>
        /// <param name="width">Width of the image. Omit width, height and depth to infer dimensions from file name.</param>
        /// <param name="height">Height of the image. Omit width, height and depth to infer dimensions from file name.</param>
        /// <param name="depth">Depth of the image. Omit width, height and depth to infer dimensions from file name.</param>
        public void MapRaw(Pi2Image image, string filename, ImageDataType dataType = ImageDataType.Unknown, int width = 0, int height = 0, int depth = 0)
        {
            PiLib.RunAndCheck(Handle, $"mapraw({image.ImageName}, {filename}, {dataType}, {width}, {height}, {depth})");
        }

        /// <summary>
        /// Writes the given image to a .raw file.
        /// Appends dimensions and extension to the file name.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="filename"></param>
        public void WriteRaw(Pi2Image img, string filename)
        {
            PiLib.RunAndCheck(Handle, $"writeraw({img.ImageName}, {filename})");
        }

        /// <summary>
        /// Writes the given image to a .tif file.
        /// No extension is appended to the file name.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="filename"></param>
        public void WriteTif(Pi2Image img, string filename)
        {
            PiLib.RunAndCheck(Handle, $"writetif({img.ImageName}, {filename})");
        }

        #endregion

        #region Tomography-related commands

        /// <summary>
        /// Performs FBP preprocessing for given projections.
        /// </summary>
        /// <param name="projections">Flat-field corrected projection data. This image is not modified.</param>
        /// <param name="preprocessed">This image will store the preprocessed projection data.</param>
        /// <param name="settings"></param>
        public void FBPPreprocess(Pi2Image projections, Pi2Image preprocessed, string settings)
        {
            PiLib.RunAndCheck(Handle, $"fbppreprocess({projections.ImageName}, {preprocessed.ImageName}, {settings})");
        }

        /// <summary>
        /// Performs backprojection part of filtered backprojection reconstruction.
        /// </summary>
        /// <param name="preprocessed">Projections processed using FBPPreprocess command.</param>
        /// <param name="output">Reconstruction will be placed in this image.</param>
        /// <param name="settings"></param>
        public void FBP(Pi2Image preprocessed, Pi2Image output, string settings)
        {
            PiLib.RunAndCheck(Handle, $"fbp({preprocessed.ImageName}, {output.ImageName}, {settings})");
        }

        /// <summary>
        /// Create an wx1x1 image that contains values of filter for use in filtered backprojection.
        /// </summary>
        /// <param name="output">Output image.</param>
        /// <param name="size">Desired size or zero to use size of output image.</param>
        /// <param name="filterType">Type of filter. Supported values are Ramp1, Ramp2, Shepp-Logan, Cosine, Hamming, and Hann.</param>
        /// <param name="cutoff">Filter cutoff frequency.</param>
        public void CreateFBPFilter(Pi2Image output, int size = 100, string filterType = "Ramp1", float cutoff = 1.0f)
        {
            PiLib.RunAndCheck(Handle, $"createfbpfilter({output.ImageName}, {size}, {filterType}, {cutoff.ToString(CultureInfo.InvariantCulture)})");
        }

        #endregion

        #region Conversions

        /// <summary>
        /// Converts img to specified data type in-place.
        /// Does not scale pixel values.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="dt"></param>
        public void Convert(Pi2Image img, ImageDataType dt)
        {
            PiLib.RunAndCheck(Handle, $"convert({img.ImageName}, {dt})");
        }

        #endregion

        #region Math

        /// <summary>
        /// Set all pixels to given value.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void Set(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"set({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Make pixels of img equal those of img2.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="img2"></param>
        public void Set(Pi2Image img, Pi2Image img2)
        {
            PiLib.RunAndCheck(Handle, $"set({img.ImageName}, {img2.ImageName})");
        }

        /// <summary>
        /// Add a constant to an image.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void Add(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"add({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Add two images and place the result to the first one.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="img2"></param>
        public void Add(Pi2Image img, Pi2Image img2)
        {
            PiLib.RunAndCheck(Handle, $"add({img.ImageName}, {img2.ImageName})");
        }

        /// <summary>
        /// Subtract a constant from an image.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void Subtract(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"subtract({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Subtract img2 from img and place the result into img.
        /// </summary>
        /// <param name="img">The first image. The result will be placed to this image.</param>
        /// <param name="img2">The seconds image.</param>
        /// <param name="allowBroadcast">Set to true to allow size of parameter image differ from size of input image. If there is a need to access pixel outside of parameter image, the nearest value inside the image is taken instead.If set to false, dimensions of input and parameter images must be equal. If set to true, the parameter image is always loaded in its entirety in distributed processing mode.</param>
        public void Subtract(Pi2Image img, Pi2Image img2, bool allowBroadcast = false)
        {
            PiLib.RunAndCheck(Handle, $"subtract({img.ImageName}, {img2.ImageName}, {allowBroadcast.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Subtract an image from a constant.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void InvSubtract(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"invsubtract({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Subtract img from img2 and place the result into img.
        /// </summary>
        /// <param name="img">The first image. The result will be placed to this image.</param>
        /// <param name="img2">The seconds image.</param>
        /// <param name="allowBroadcast">Set to true to allow size of parameter image differ from size of input image. If there is a need to access pixel outside of parameter image, the nearest value inside the image is taken instead.If set to false, dimensions of input and parameter images must be equal. If set to true, the parameter image is always loaded in its entirety in distributed processing mode.</param>
        public void InvSubtract(Pi2Image img, Pi2Image img2, bool allowBroadcast = false)
        {
            PiLib.RunAndCheck(Handle, $"invsubtract({img.ImageName}, {img2.ImageName}, {allowBroadcast.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Multiply image by a constant.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void Multiply(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"multiply({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Multiply two images and place the result into the first image.
        /// </summary>
        /// <param name="img">The first image. The result will be placed to this image.</param>
        /// <param name="img2">The seconds image.</param>
        /// <param name="allowBroadcast">Set to true to allow size of parameter image differ from size of input image. If there is a need to access pixel outside of parameter image, the nearest value inside the image is taken instead.If set to false, dimensions of input and parameter images must be equal. If set to true, the parameter image is always loaded in its entirety in distributed processing mode.</param>
        public void Multiply(Pi2Image img, Pi2Image img2, bool allowBroadcast = false)
        {
            PiLib.RunAndCheck(Handle, $"multiply({img.ImageName}, {img2.ImageName}, {allowBroadcast.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Divide image by a constant.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="value"></param>
        public void Divide(Pi2Image img, float value)
        {
            PiLib.RunAndCheck(Handle, $"divide({img.ImageName}, {value.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Divide img by img2 and place the result to img.
        /// </summary>
        /// <param name="img">The first image. The result will be placed to this image.</param>
        /// <param name="img2">The seconds image.</param>
        /// <param name="allowBroadcast">Set to true to allow size of parameter image differ from size of input image. If there is a need to access pixel outside of parameter image, the nearest value inside the image is taken instead.If set to false, dimensions of input and parameter images must be equal. If set to true, the parameter image is always loaded in its entirety in distributed processing mode.</param>
        public void Divide(Pi2Image img, Pi2Image img2, bool allowBroadcast = false)
        {
            PiLib.RunAndCheck(Handle, $"divide({img.ImageName}, {img2.ImageName}, {allowBroadcast.ToString(CultureInfo.InvariantCulture)})");
        }

        #endregion

        #region Generation

        /// <summary>
        /// Set value of pixel at specified location.
        /// </summary>
        /// <param name="img"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="value"></param>
        public void Set(Pi2Image img, int x, int y, int z, float value)
        {
            PiLib.RunAndCheck(Handle, String.Format(CultureInfo.InvariantCulture, "set({0}, [{1}, {2}, {3}], {4})", img.ImageName, x, y, z, value));
        }

        /// <summary>
        /// Generate a grayscale ramp.
        /// </summary>
        /// <param name="img">Target image. Existing pixel values will be lost.</param>
        /// <param name="dimension">Dimension where gray values increase.</param>
        public void Ramp(Pi2Image img, int dimension)
        {
            PiLib.RunAndCheck(Handle, $"ramp({img.ImageName}, {dimension})");
        }

        /// <summary>
        /// Generate additive Gaussian noise.
        /// </summary>
        /// <param name="img">Target image where noise will be added.</param>
        /// <param name="mean">Mean of the noise.</param>
        /// <param name="stddev">Standard deviation of the noise.</param>
        public void Noise(Pi2Image img, float mean = 0.0f, float stddev = 50.0f)
        {
            PiLib.RunAndCheck(Handle, $"noise({img.ImageName}, {mean.ToString(CultureInfo.InvariantCulture)}, {stddev.ToString(CultureInfo.InvariantCulture)})");
        }

        #endregion

        #region Transformations

        /// <summary>
        /// Scales input image and stores result to output image.
        /// </summary>
        /// <param name="inImg">Input image.</param>
        /// <param name="outImg">Output image. If scale factor is zero, the input image is scaled to the size of the output image.</param>
        /// <param name="scaleFactor">Scaling factor. If nonzero, the output image size is calculated from scaling factor and input image size. If zero, the input image is scaled to the size of the output image.</param>
        public void Scale(Pi2Image inImg, Pi2Image outImg, float scaleFactor = 0)
        {
            // TODO: Vector scaling factor, interpolation mode, boundary condition.
            PiLib.RunAndCheck(Handle, $"scale({inImg.ImageName}, {outImg.ImageName}, {scaleFactor.ToString(CultureInfo.InvariantCulture)})");
        }

        /// <summary>
        /// Rotates input image 90 degrees clockwise.
        /// </summary>
        /// <param name="inImg"></param>
        /// <param name="outImg"></param>
        public void Rot90CW(Pi2Image inImg, Pi2Image outImg)
        {
            PiLib.RunAndCheck(Handle, $"rot90cw({inImg.ImageName}, {outImg.ImageName})");
        }

        /// <summary>
        /// Copies pixels of source image to target image.
        /// Converts data type if required.
        /// </summary>
        /// <param name="source"></param>
        /// <param name="target"></param>
        public void Copy(Pi2Image source, Pi2Image target)
        {
            PiLib.RunAndCheck(Handle, $"copy({source.ImageName}, {target.ImageName})");
        }

        /// <summary>
        /// Copies pixels of the source image to the target image to the specified location.
        /// Does not change the size of the target image.
        /// </summary>
        /// <param name="target"></param>
        /// <param name="source"></param>
        /// <param name="position"></param>
        public void Copy(Pi2Image source, Pi2Image target, Vec3 position)
        {
            PiLib.RunAndCheck(Handle, $"copy({source.ImageName}, {target.ImageName}, {ToIntVec(position)})");
        }

        /// <summary>
        /// Crops the image into size of output. Set size of output before calling this command or pass size as an argument.
        /// NOTE: If the data type of target image is not correct, this function will change it!
        /// </summary>
        /// <param name="input">Input image.</param>
        /// <param name="output">Output image.</param>
        /// <param name="position">Position in input image where the top-left corner of the cropped image is placed.</param>
        /// <param name="size">Size of output image. Specify zeroes or nothing to crop to current size of output image.</param>
        public void Crop(Pi2Image input, Pi2Image output, Vec3 position, Vec3 size)
        {
            PiLib.RunAndCheck(Handle, $"crop({input.ImageName}, {output.ImageName}, {ToIntVec(position)}, {ToIntVec(size)})");
        }

        #endregion

        #region Projections

        /// <summary>
        /// Mean projection in some coordinate direction.
        /// </summary>
        /// <param name="source"></param>
        /// <param name="target">Output image. Data type must be float32.</param>
        /// <param name="dimension"></param>
        public void MeanProject(Pi2Image source, Pi2Image target, int dimension)
        {
            PiLib.RunAndCheck(Handle, $"meanproject({source.ImageName}, {target.ImageName}, {dimension})");
        }

        #endregion

        /// <summary>
        /// Converts Vec3 to string [x, y, z] where the components are integers.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static string ToIntVec(Vec3 x)
        {
            return "[" + ((int)Math.Round(x.X)).ToString(CultureInfo.InvariantCulture) + ", " + ((int)Math.Round(x.Y)).ToString(CultureInfo.InvariantCulture) + ", " + ((int)Math.Round(x.Z)).ToString(CultureInfo.InvariantCulture) + "]";
        }

        // TODO: Add pilib commands as methods here, translate from Pi2Images to image names etc.
    }
}
