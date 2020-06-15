using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing.Imaging;

namespace pi2cs
{
    /// <summary>
    /// Picture box that retrieves its image from the Pi system.
    /// When the picture should be grabbed from the Pi system, the UpdateImage method must be called.
    /// </summary>
    public partial class Pi2PictureBox : PictureBox
    {
        /// <summary>
        /// Pi image to show.
        /// </summary>
        public Pi2Image PiImage
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets a value indicating if the control should be automatically resized to the size of the image.
        /// </summary>
        public bool AutoResize
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets the zoom level.
        /// </summary>
        public float Zoom
        {
            get;
            set;
        }

        /// <summary>
        /// Sets the position in the Pi2 picture corresponding to the top left of the picture box.
        /// </summary>
        public PointF PicturePosition
        {
            get;
            set;
        }

        /// <summary>
        /// Gray value that is mapped to black.
        /// </summary>
        public float Min
        {
            get;
            set;
        }

        /// <summary>
        /// Gray value that is mapped to white.
        /// </summary>
        public float Max
        {
            get;
            set;
        }

        /// <summary>
        /// Holds version of original without contrast and zoom adjustment.
        /// </summary>
        private List<float> OriginalBitmap
        {
            get;
            set;
        }

        /// <summary>
        /// Width of the original image.
        /// </summary>
        public int OriginalWidth
        {
            get;
            private set;
        }

        /// <summary>
        /// Depth of the original image.
        /// </summary>
        public int OriginalDepth
        {
            get;
            private set;
        }

        /// <summary>
        /// Height of the original image.
        /// </summary>
        public int OriginalHeight
        {
            get;
            private set;
        }

        /// <summary>
        /// Original pixel data type of the image.
        /// </summary>
        public ImageDataType OriginalDataType
        {
            get;
            private set;
        }

        /// <summary>
        /// Minimum pixel value in the original image.
        /// </summary>
        public float GlobalMin
        {
            get;
            private set;
        }

        /// <summary>
        /// Maximum pixel value in the original image.
        /// </summary>
        public float GlobalMax
        {
            get;
            private set;
        }

        /// <summary>
        /// Minimum possible value in the original image.
        /// </summary>
        public float DynamicMin
        {
            get;
            private set;
        }

        /// <summary>
        /// Maximum possible value in the original image.
        /// </summary>
        public float DynamicMax
        {
            get;
            private set;
        }

        /// <summary>
        /// Index of the slice that is currently visible.
        /// After slice is set, one must call UpdateImage to update the view.
        /// </summary>
        public int Slice
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets a list of annotations drawn on top of the image.
        /// </summary>
        public List<Annotation> Annotations
        {
            get;
            set;
        }

        /// <summary>
        /// Gets the histogram of this image.
        /// </summary>
        public IReadOnlyList<PointF> Histogram
        {
            get;
            private set;
        }

        /// <summary>
        /// Backing field for histogram array.
        /// </summary>
        private PointF[] histogram;

        /// <summary>
        /// Constructor
        /// </summary>
        public Pi2PictureBox() : this(null)
        {
            Annotations = new List<Annotation>();
        }

        protected override void OnSizeChanged(EventArgs e)
        {
            base.OnSizeChanged(e);
            UpdateImage();
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="img"></param>
        /// <param name="zoom"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        public Pi2PictureBox(Pi2Image img, float zoom = 1, float min = 0, float max = 255)
        {
            InitializeComponent();
            OriginalBitmap = new List<float>();
            PiImage = img;
            Zoom = zoom;
            Min = min;
            Max = max;
            UpdateHistogram();
        }

        private double ScreenToPictureY(double y)
        {
            return (y - Height / 2) / Zoom + PicturePosition.Y + OriginalHeight / 2;
        }

        private double ScreenToPictureX(double x)
        {
            return (x - Width / 2) / Zoom + PicturePosition.X + OriginalWidth / 2;
        }

        /// <summary>
        /// Converts from screen coordinates to picture coordinates.
        /// </summary>
        /// <param name="s">Point in screen coordinates.</param>
        /// <returns>Point in picture coordinates.</returns>
        public Vec2 ScreenToPicture(Vec2 s)
        {
            return new Vec2(ScreenToPictureX(s.X), ScreenToPictureY(s.Y));
        }

        /// <summary>
        /// Converts from picture coordinates to screen coordinates.
        /// </summary>
        /// <param name="p">Point in picture coordinates</param>
        /// <returns>Point in screen coordinates.</returns>
        public Vec2 PictureToScreen(Vec2 p)
        {
            return new Vec2((p.X - PicturePosition.X - OriginalWidth / 2) * Zoom + Width / 2,
                (p.Y - PicturePosition.Y - OriginalHeight / 2) * Zoom + Height / 2);
        }


        /// <summary>
        /// Gets pixel value at specified location.
        /// Location is given in screen coordinates, and this method accounts for possible non-unity zoom factor.
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public float GetValue(Vec2 p)
        {
            return GetValueImage(ScreenToPicture(p));
        }

        /// <summary>
        /// Gets pixel value at specified location.
        /// The location is given in image coordinates.
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public float GetValueImage(Vec2 p)
        {
            int origX = (int)Math.Round(p.X);
            int origY = (int)Math.Round(p.Y);
            if (origX < 0 || origY < 0 || origX >= OriginalWidth || origY >= OriginalHeight)
                return float.NaN;

            return OriginalBitmap[origX + origY * OriginalWidth];
        }


        /// <summary>
        /// Copies image data from temporary buffer to picture box, obeying zoom and contrast settings.
        /// </summary>
        public void UpdateScreen()
        {
            if (AutoResize)
            {
                Width = (int)Math.Round(OriginalWidth * Zoom);
                Height = (int)Math.Round(OriginalHeight * Zoom);

                if (Width <= 0)
                    Width = 1;
                if (Height <= 0)
                    Height = 1;
            }

            if (Width > 0 && Height > 0)
            {

                if (Image == null || Image.Width != Width || Image.Height != Height)
                {
                    Bitmap b = new Bitmap(Width, Height, PixelFormat.Format32bppArgb);

                    Image = b;
                }

                Bitmap bitmap = (Bitmap)Image;

                //BitmapData bdata = bitmap.LockBits(new Rectangle(0, 0, targetWidth, targetHeight), ImageLockMode.WriteOnly, PixelFormat.Format8bppIndexed);
                BitmapData bdata = bitmap.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);

                try
                {

                    unsafe
                    {
                        byte* po = (byte*)bdata.Scan0.ToPointer();

                        // Assign pixel values
                        for (int y = 0; y < Height; y++)
                        {
                            int srcy = (int)Math.Round(ScreenToPictureY(y));
                            for (int x = 0; x < Width; x++)
                            {
                                int srcx = (int)Math.Round((x - Width / 2) / Zoom + PicturePosition.X + OriginalWidth / 2);

                                if (srcx >= 0 && srcx < OriginalWidth && srcy >= 0 && srcy < OriginalHeight)
                                {
                                    float val = OriginalBitmap[srcx + srcy * OriginalWidth];

                                    // Transform to [0, 1]
                                    val = (val - Min) / (Max - Min);

                                    if (val > 1)
                                        val = 1;
                                    else if (val < 0)
                                        val = 0;

                                    // Output is 32-bit ARGB
                                    *(po + y * bdata.Stride + x * 4) = (byte)(val * 255);
                                    *(po + y * bdata.Stride + x * 4 + 1) = (byte)(val * 255);
                                    *(po + y * bdata.Stride + x * 4 + 2) = (byte)(val * 255);
                                    *(po + y * bdata.Stride + x * 4 + 3) = (byte)255;
                                }
                                else
                                {
                                    const byte DARK_SHADE = 180;
                                    const byte LIGHT_SHADE = 240;

                                    // Find out if we are in even or odd square of size S and
                                    // decide to use dark or light shade of color based on that,
                                    // in order to create a checkerboard pattern background.
                                    const int S = 15;
                                    int xs = (x / S) % 2;
                                    int ys = (y / S) % 2;

                                    byte c = DARK_SHADE;
                                    if ((xs == 0 && ys == 0) ||
                                        (xs == 1 && ys == 1))
                                        c = LIGHT_SHADE;

                                    *(po + y * bdata.Stride + x * 4) = c;
                                    *(po + y * bdata.Stride + x * 4 + 1) = c;
                                    *(po + y * bdata.Stride + x * 4 + 2) = c;
                                    *(po + y * bdata.Stride + x * 4 + 3) = 255;
                                }
                            }
                        }
                    }
                }
                finally
                {
                    bitmap.UnlockBits(bdata);
                }

                Refresh();
            }
        }

        /// <summary>
        /// Re-calculates histogram of the shown image.
        /// </summary>
        private void UpdateHistogram()
        {
            if(histogram == null || histogram.Length != 256)
            {
                histogram = new PointF[256];

                Histogram = histogram;
            }

            if (OriginalBitmap != null && OriginalBitmap.Count > 0)
            {
                // Find global min and max, skip inf and nan pixels.
                GlobalMin = float.PositiveInfinity;
                GlobalMax = float.NegativeInfinity;

                foreach(float pix in OriginalBitmap)
                {
                    if(!float.IsNaN(pix) && !float.IsInfinity(pix))
                    {
                        GlobalMin = Math.Min(GlobalMin, pix);
                        GlobalMax = Math.Max(GlobalMax, pix);
                    }
                }

                if(float.IsInfinity(GlobalMin) && float.IsInfinity(GlobalMax))
                {
                    GlobalMin = 0;
                    GlobalMax = 1;
                }
                else if(float.IsInfinity(GlobalMin))
                {
                    GlobalMin = GlobalMax - 1;
                }
                else if(float.IsInfinity(GlobalMax))
                {
                    GlobalMax = GlobalMin + 1;
                }


                // Init bin centers and set counts to zero.
                for (int n = 0; n < histogram.Length; n++)
                {
                    float binCenter = GlobalMin + (n + 0.5f) * (GlobalMax - GlobalMin) / histogram.Length;
                    histogram[n].X = binCenter;
                    histogram[n].Y = 0;
                }

                // Calculate histogram
                for (int n = 0; n < OriginalBitmap.Count; n++)
                {
                    float pix = OriginalBitmap[n];

                    if (!float.IsNaN(pix) && !float.IsInfinity(pix))
                    {
                        int bin = (int)Math.Round((pix - GlobalMin) / (GlobalMax - GlobalMin) * (histogram.Length - 1));
                        if (bin < 0)
                            bin = 0;
                        else if (bin > histogram.Length - 1)
                            bin = histogram.Length - 1;

                        histogram[bin].Y++;
                    }
                }
            }
        }


        /// <summary>
        /// Copies image data from given pi2 pointer to array in this picture box.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="depth"></param>
        /// <param name="dt"></param>
        private void CopyFromPointerToTemp(IntPtr data, long width, long height, long depth, ImageDataType dt)
        {
            if (width > int.MaxValue)
                throw new InvalidOperationException("Image width is too large to be shown. Maximum supported width is " + int.MaxValue);
            if (height > int.MaxValue)
                throw new InvalidOperationException("Image height is too large to be shown. Maximum supported height is " + int.MaxValue);
            if (depth > int.MaxValue)
                throw new InvalidOperationException("Image depth is too large to be shown. Maximum supported depth is " + int.MaxValue);

            long capacity = width * height;

            if (capacity > int.MaxValue)
                throw new InvalidOperationException("Size of single slice is too large to be shown. Maximum supported count of pixels is " + int.MaxValue);

            OriginalBitmap.Clear();
            OriginalBitmap.Capacity = (int)capacity;
            for (int n = 0; n < width * height; n++)
                OriginalBitmap.Add(0);
            OriginalWidth = (int)width;
            OriginalHeight = (int)height;
            OriginalDepth = (int)depth;

            OriginalDataType = dt;

            if (Slice < 0)
                Slice = 0;
            else if (Slice >= depth)
                Slice = (int)(depth - 1);

            if (data != IntPtr.Zero)
            {
                unsafe
                {
                    // Construct pixel getter function that converts all pixel data types to float.
                    byte* pi8 = (byte*)data.ToPointer();
                    UInt16* pi16 = (UInt16*)pi8;
                    UInt32* pi32 = (UInt32*)pi8;
                    UInt64* pi64 = (UInt64*)pi8;
                    float* pif = (float*)pi8;

                    Func<int, int, float> getPixel;

                    switch (dt)
                    {
                        case ImageDataType.UInt8:
                            getPixel = (int x, int y) => pi8[x + y * width + Slice * width * height];
                            DynamicMin = uint.MinValue;
                            DynamicMax = uint.MaxValue;
                            break;
                        case ImageDataType.UInt16:
                            getPixel = (int x, int y) => pi16[x + y * width + Slice * width * height];
                            DynamicMin = UInt16.MinValue;
                            DynamicMax = UInt16.MaxValue;
                            break;
                        case ImageDataType.UInt32:
                            getPixel = (int x, int y) => pi32[x + y * width + Slice * width * height];
                            DynamicMin = UInt32.MinValue;
                            DynamicMax = UInt32.MaxValue;
                            break;
                        case ImageDataType.UInt64:
                            getPixel = (int x, int y) => pi64[x + y * width + Slice * width * height];
                            DynamicMin = UInt64.MinValue;
                            DynamicMax = UInt64.MaxValue;
                            break;
                        case ImageDataType.Float32:
                            getPixel = (int x, int y) => pif[x + y * width + Slice * width * height];
                            // These are updated later as MinValue and MaxValue are not practically very usable choices.
                            //DynamicMin = float.MinValue;
                            //DynamicMax = float.MaxValue;
                            break;
                        default:
                            getPixel = (int x, int y) => 0;
                            DynamicMin = 0;
                            DynamicMax = 0;
                            break;
                    }

                    // Copy pixel values to the array.
                    int iw = (int)width;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            float val = getPixel(x, y);
                            OriginalBitmap[x + y * iw] = val;
                        }
                    }
                }
            }

            UpdateHistogram();

            if(dt == ImageDataType.Float32)
            {
                float headroom = 0.1f * (GlobalMax - GlobalMin);
                DynamicMin = GlobalMin - headroom;
                DynamicMax = GlobalMax + headroom;
            }
        }

        /// <summary>
        /// Updates the image in the picture box from the Pi system.
        /// Image must be available in the Pi system when this method is executing.
        /// </summary>
        public void UpdateImage()
        {
            if(PiImage != null)
            {
                long w = 0;
                long h = 0;
                long d = 0;
                ImageDataType dt;
                IntPtr data = PiImage.GetData(out w, out h, out d, out dt);

                CopyFromPointerToTemp(data, w, h, d, dt);
            }
            else
            {
                CopyFromPointerToTemp(IntPtr.Zero, 0, 0, 0, ImageDataType.UInt8);
            }
            UpdateScreen();
        }


        #region Picking

        /// <summary>
        /// Buffer where annotation pick information is drawn.
        /// </summary>
        private Bitmap pickBuffer;

        private int maxPickColor;
        private BiDictionary<Tuple<Annotation, int>, Color> pickColors = new BiDictionary<Tuple<Annotation, int>, Color>();


        /// <summary>
        /// Gets color that must be used to map to the pick map.
        /// </summary>
        /// <param name="annotation"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        private Color GetPickColor(Annotation annotation, int point)
        {
            var key = Tuple.Create(annotation, point);
            if (!pickColors.Contains(key))
            {
                int color = maxPickColor + 1;
                maxPickColor++;
                pickColors.Add(key, Color.FromArgb(color));
            }

            return pickColors[key];
        }

        /// <summary>
        /// Gets pen that can be used to draw to pick map.
        /// </summary>
        /// <param name="annotation"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        internal Pen GetPickPen(Annotation annotation, int point = -1)
        {
            Color color = GetPickColor(annotation, point);
            return new Pen(color, 3);
        }

        /// <summary>
        /// Gets brush that can be used to draw to pick map.
        /// </summary>
        /// <param name="annotation"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        internal Brush GetPickBrush(Annotation annotation, int point = -1)
        {
            Color color = GetPickColor(annotation, point);
            return new SolidBrush(color);
        }

        /// <summary>
        /// Gets a tuple identifying annotation and control point index at the specified location.
        /// </summary>
        /// <param name="pos">Position in client coordinates.</param>
        /// <returns></returns>
        public Tuple<Annotation, int> FindAnnotation(Point pos)
        {
            Color c = pickBuffer.GetPixel(pos.X, pos.Y);
            if (pickColors.Contains(c))
                return pickColors[c];

            return null;
        }

        #endregion

        protected override void OnPaint(PaintEventArgs pe)
        {
            base.OnPaint(pe);

            if(pickBuffer != null && (pickBuffer.Width != ClientSize.Width || pickBuffer.Height != ClientSize.Height))
            {
                pickBuffer.Dispose();
                pickBuffer = null;
            }

            if (pickBuffer == null)
                pickBuffer = new Bitmap(ClientSize.Width, ClientSize.Height, PixelFormat.Format32bppRgb);

            using (Graphics pickg = Graphics.FromImage(pickBuffer))
            {
                pickg.Clear(Color.FromArgb(0));
                pickColors.Clear();
                maxPickColor = Color.FromArgb(255, 0, 0).ToArgb();

                // Draw annotations
                foreach (Annotation a in Annotations)
                    a.Draw(this, pe.Graphics, pickg, Slice);

                // For debugging:
                //pe.Graphics.Clear(Color.FromArgb(0));
                //pe.Graphics.DrawImage(pickBuffer, Point.Empty);
            }
        }
    }
}
