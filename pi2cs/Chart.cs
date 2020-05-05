using System;
using System.Collections.Generic;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Windows.Forms;
using System.Drawing.Drawing2D;

namespace pi2cs
{
    /// <summary>
    /// Shows simple chart.
    /// This control is used instead of Windows Forms chart as that one is not supported in Mono.
    /// </summary>
    public partial class Chart : Control
    {
        /// <summary>
        /// List of all data series shown in this control.
        /// </summary>
        public List<DataSeries> DataSeries
        {
            get;
            private set;
        }

        /// <summary>
        /// Properties of the X-axis.
        /// </summary>
        public Axis XAxis
        {
            get;
            private set;
        }

        /// <summary>
        /// Properties of the Y-axis.
        /// </summary>
        public Axis YAxis
        {
            get;
            private set;
        }

        /// <summary>
        /// Gets and sets the location and visibility of the legend box.
        /// </summary>
        public LegendLocation LegendLocation
        {
            get;
            set;
        }

        /// <summary>
        /// Margin between axis labels and other elements
        /// </summary>
        public double AxisMargin
        {
            get;
            set;
        } = 10;

        /// <summary>
        /// Margin between items in the legend box.
        /// </summary>
        public double LegendMargin
        {
            get;
            set;
        } = 5;

        /// <summary>
        /// Font used when drawing axis labels.
        /// The control assumes ownership of the assigned object and will dispose it when the control is disposed.
        /// </summary>
        public Font AxisLabelFont
        {
            get;
            set;
        }

        /// <summary>
        /// Font used to draw axis tick labels.
        /// The control assumes ownership of the assigned object and will dispose it when the control is disposed.
        /// </summary>
        public Font TickLabelFont
        {
            get;
            set;
        }

        /// <summary>
        /// Size of tick marks.
        /// </summary>
        public double TickSize
        {
            get;
            set;
        } = 4;

        /// <summary>
        /// Pen used to draw axis tick marks etc.
        /// </summary>
        public Pen AxisPen
        {
            get;
            set;
        }

        /// <summary>
        /// Formats to be used when drawing labels etc.
        /// </summary>
        private StringFormat AxisLabelFormat;
        private StringFormat YTickFormat;
        private StringFormat LegendFormat;

        /// <summary>
        /// Constructor
        /// </summary>
        public Chart()
        {
            DoubleBuffered = true;

            BackColor = Color.White;
            ColorOrder = new List<Color>{ Color.Black, Color.Red, Color.Blue, Color.Brown, Color.Gray };
            DataSeries = new List<DataSeries>();
            XAxis = new Axis(0, 1, "X [arb.]");
            YAxis = new Axis(0, 1, "Y [arb.]");
            YAxis.Flipped = true;

            
            AxisLabelFont = new Font(SystemFonts.DefaultFont.FontFamily, 1.25f * SystemFonts.DefaultFont.Size, FontStyle.Regular);
            TickLabelFont = new Font(AxisLabelFont.FontFamily, 1.0f * AxisLabelFont.Size, FontStyle.Regular);
            AxisPen = Pens.Black;

            AxisLabelFormat = new StringFormat();
            AxisLabelFormat.Alignment = StringAlignment.Center;
            AxisLabelFormat.LineAlignment = StringAlignment.Near;

            YTickFormat = new StringFormat();
            YTickFormat.Alignment = StringAlignment.Far;
            YTickFormat.LineAlignment = StringAlignment.Center;

            LegendFormat = new StringFormat();
            LegendFormat.Alignment = StringAlignment.Near;
            LegendFormat.LineAlignment = StringAlignment.Center;

            LegendLocation = LegendLocation.TopRight;
        }


        /// <summary>
        /// Order of appearance of colors in auto-created data series.
        /// </summary>
        public List<Color> ColorOrder
        {
            get;
            set;
        }

        private int nextColor = 0;

        /// <summary>
        /// Gets next unused color in the color order.
        /// </summary>
        /// <returns></returns>
        public Color NextColor()
        {
            if (ColorOrder == null || ColorOrder.Count <= 0)
                throw new InvalidOperationException("No colors in ColorOrder.");

            if (nextColor <= ColorOrder.Count - 1)
                return ColorOrder[nextColor++];

            nextColor = 0;
            return ColorOrder.Last();
        }

        /// <summary>
        /// Create a new data series
        /// </summary>
        /// <returns></returns>
        public DataSeries NewSeries(string title = null)
        {
            DataSeries.Add(new DataSeries(NextColor(), title));
            return DataSeries.Last();
        }

        /// <summary>
        /// Draws the data series
        /// </summary>
        private void DrawData(Graphics g, Rectangle dataBounds)
        {
            foreach(DataSeries series in DataSeries)
                series.Draw(g, dataBounds, XAxis, YAxis);
        }

        /// <summary>
        /// Calculates height of tick labels.
        /// </summary>
        /// <param name="g"></param>
        /// <returns></returns>
        private double TickLabelHeight(Graphics g)
        {
            return g.MeasureString("0123456789.", TickLabelFont).Height;
        }
        
        /// <summary>
        /// Measures height of an axis.
        /// ticks list is used only if useTickHeight is false (i.e. measuring width of y-axis)
        /// </summary>
        private double MeasureAxisHeight(Graphics g, Axis axis, bool useTickHeight, List<double> ticks)
        {
            if (!String.IsNullOrEmpty(axis.Label))
            {
                SizeF labelSize = g.MeasureString(axis.Label, AxisLabelFont);

                double tickSize = 0;
                if (useTickHeight)
                {
                    tickSize = TickLabelHeight(g);
                }
                else
                {
                    // Measure each tick label separately
                    foreach(double x in ticks)
                    {
                        tickSize = Math.Max(tickSize, g.MeasureString(GetTickLabel(x), TickLabelFont).Width);
                    }
                }

                return labelSize.Height + 2 * AxisMargin + tickSize + TickSize;
            }
            else
            {
                return AxisMargin;
            }
        }

        /// <summary>
        /// Rounds d to given number of significant digits.
        /// </summary>
        /// <param name="d"></param>
        /// <param name="digits"></param>
        /// <returns></returns>
        private static double RoundToSignificantDigits(double d, int digits)
        {
            if (d == 0)
                return 0;

            double scale = Math.Pow(10, Math.Floor(Math.Log10(Math.Abs(d))) + 1);
            return scale * Math.Round(d / scale, digits);
        }

        /// <summary>
        /// Tests if two number ranges overlap.
        /// Assumes range.X is less than or equal to range.Y for both ranges.
        /// </summary>
        /// <param name="range1"></param>
        /// <param name="range2"></param>
        /// <returns></returns>
        private static bool Overlaps(Vec2 range1, Vec2 range2)
        {
            return range1.X <= range2.Y && range2.X <= range1.Y;
        }

        /// <summary>
        /// Convert tick position to tick label.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private string GetTickLabel(double x)
        {
            return x.ToString();
        }

        /// <summary>
        /// Finds suitable locations for tick marks.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="axis"></param>
        /// <param name="imageBounds"></param>
        /// <param name="horizontal"></param>
        /// <returns></returns>
        private List<double> FindTickLocations(Graphics g, Axis axis, Vec2 imageBounds, bool horizontal)
        {
            List<double> ticks = new List<double>();

            int tickCount = 20;
            while (true)
            {
                // Gap between the ticks
                double delta = RoundToSignificantDigits((axis.Bounds.Y - axis.Bounds.X) / tickCount, 2);
                // Location of the first tick.
                double start;
                
                if(axis.Bounds.X <= 0 && axis.Bounds.Y >= 0)
                {
                    // Ensure that 0 becomes a tick.
                    start = 0;
                    while (start > axis.Bounds.X)
                        start -= delta;

                    start += delta;
                }
                else
                {
                    // Choose whatever tick range
                    start = RoundToSignificantDigits(axis.Bounds.X, 2);

                    // Ensure the first tick is inside the axis range.
                    while (start < axis.Bounds.X)
                        start += delta;
                }

                // Make ticks, round the positions to nice values
                double scale = Math.Max(Math.Abs(axis.Bounds.X), Math.Abs(axis.Bounds.Y)) / 1000.0;
                int digits = -(int)Math.Round(Math.Log10(scale));
                if (digits < 0)
                    digits = 0;

                ticks.Clear();
                for (double x = start; x < axis.Bounds.Y * 1.001; x += delta)
                {
                    ticks.Add(Math.Round(x, digits));
                }

                // Test if the ticks fit
                Vec2 lastRange = new Vec2(Double.NegativeInfinity, Double.NegativeInfinity);
                bool fits = true;
                for(int n = 0; n < ticks.Count; n++)
                {
                    string tickLabel = GetTickLabel(ticks[n]);
                    SizeF tickSize = g.MeasureString(tickLabel, TickLabelFont);
                    double tickX = axis.DataToImage(ticks[n], imageBounds.X, imageBounds.Y);
                    double tickRadius = horizontal ? tickSize.Width / 2.0 : 1.5 * tickSize.Height / 2.0; // NOTE: 1.5 leaves a bit more room for the tick labels
                    double tickLeft = tickX - tickRadius;
                    double tickRight = tickX + tickRadius;

                    Vec2 newRange = new Vec2(Math.Min(tickLeft, tickRight), Math.Max(tickLeft, tickRight));
                    
                    if(Overlaps(lastRange, newRange))
                    {
                        fits = false;
                        break;
                    }
                    
                    lastRange = newRange;
                }

                if (fits)
                    break;

                tickCount--;

                if(tickCount == 0)
                {
                    // No ticks fit ?!
                    ticks.Clear();
                    break;
                }
            }

            return ticks;
        }

        /// <summary>
        /// Draws x- and y-axes.
        /// </summary>
        private void DrawAxes(Graphics g, Rectangle dataBounds, List<double> xTicks, List<double> yTicks)
        {
            double tickLabelHeight = TickLabelHeight(g);

            using (Brush brush = new SolidBrush(AxisPen.Color))
            {
                if (!String.IsNullOrEmpty(XAxis.Label))
                {
                    int cx = dataBounds.Left + dataBounds.Width / 2;
                    int xy = (int)Math.Round(dataBounds.Bottom + TickSize + tickLabelHeight + AxisMargin);
                    g.DrawString(XAxis.Label, AxisLabelFont, brush, cx, xy, AxisLabelFormat);
                }

                int tickY = (int)Math.Round(dataBounds.Bottom + TickSize);
                
                // Draw ticks
                foreach (double x in xTicks)
                {
                    int imgX = (int)XAxis.DataToImage(x, dataBounds.Left, dataBounds.Right);
                    g.DrawLine(AxisPen, imgX, dataBounds.Bottom, imgX, (int)Math.Round(dataBounds.Bottom + TickSize));

                    g.DrawString(GetTickLabel(x), TickLabelFont, brush, imgX, tickY, AxisLabelFormat);
                }

                if (!String.IsNullOrEmpty(YAxis.Label))
                {
                    int cy = dataBounds.Top + dataBounds.Height / 2;
                    int yx = (int)Math.Round(AxisMargin);

                    Matrix f = g.Transform;
                    try
                    {
                        g.TranslateTransform(yx, cy);
                        g.RotateTransform(-90);

                        g.DrawString(YAxis.Label, AxisLabelFont, brush, 0, 0, AxisLabelFormat);
                    }
                    finally
                    {
                        g.Transform = f;
                    }
                }

                int tickX = (int)Math.Round(dataBounds.Left - TickSize);
                foreach (double y in yTicks)
                {
                    int imgY = (int)YAxis.DataToImage(y, dataBounds.Top, dataBounds.Bottom);
                    g.DrawLine(AxisPen, (int)Math.Round(dataBounds.Left - TickSize), imgY, dataBounds.Left, imgY);

                    g.DrawString(GetTickLabel(y), TickLabelFont, brush, tickX, imgY, YTickFormat);
                }
            }
        }

        /// <summary>
        /// Draws legend box.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="dataBounds"></param>
        private void DrawLegend(Graphics g, Rectangle dataBounds)
        {
            if (LegendLocation == LegendLocation.None)
                return;

            // Measure legend size
            double glyphWidth = 25;
            double glyphHeight = 15;
            double w = 0;
            double maxTextHeight = 0;
            double lineCount = 0;
            foreach(DataSeries series in DataSeries)
            {
                if(series.Title != null)
                {
                    SizeF size = g.MeasureString(series.Title, AxisLabelFont);
                    w = Math.Max(w, size.Width + 3 * LegendMargin + glyphWidth);
                    maxTextHeight = Math.Max(maxTextHeight, size.Height);
                    lineCount++;
                }
            }
            double lineHeight = Math.Max(glyphHeight, maxTextHeight);
            double h = LegendMargin + (lineHeight + LegendMargin) * lineCount;

            if (w > 0)
            {
                // Calculate legend box location
                double x = 0, y = 0;
                switch(LegendLocation)
                {
                    case LegendLocation.TopLeft: x = dataBounds.Left + AxisMargin; y = dataBounds.Top + AxisMargin; break;
                    case LegendLocation.TopRight: x = dataBounds.Right - AxisMargin - w; y = dataBounds.Top + AxisMargin; break;
                    case LegendLocation.BottomLeft: x = dataBounds.Left + AxisMargin; y = dataBounds.Bottom - AxisMargin - h; break;
                    case LegendLocation.BottomRight: x = dataBounds.Right - AxisMargin - w; y = dataBounds.Bottom - AxisMargin - h; break;
                }
                Rectangle legendRect = new Rectangle((int)Math.Round(x), (int)Math.Round(y), (int)Math.Round(w), (int)Math.Round(h));

                // Draw box
                using (Brush brush = new SolidBrush(BackColor))
                {
                    g.FillRectangle(brush, legendRect);
                    g.DrawRectangle(AxisPen, legendRect);
                }

                // Draw items
                using (Brush brush = new SolidBrush(AxisPen.Color))
                {
                    int n = 0;
                    foreach (DataSeries series in DataSeries)
                    {
                        if (series.Title != null)
                        {
                            double lineY = legendRect.Top + LegendMargin + n * (lineHeight + LegendMargin) + 0.5 * lineHeight;

                            // Glyph
                            double glyphX = legendRect.Left + LegendMargin;
                            Rectangle glyphRect = new Rectangle((int)Math.Round(glyphX), (int)Math.Round(lineY - 0.5 * glyphHeight), (int)Math.Round(glyphWidth), (int)Math.Round(glyphHeight));
                            //g.DrawRectangle(Pens.Red, glyphRect);
                            g.DrawLine(series.LinePen, (int)Math.Round(glyphX), (int)Math.Round(lineY), (int)Math.Round(glyphX + glyphWidth), (int)Math.Round(lineY));

                            // Text
                            double textX = legendRect.Left + 2 * LegendMargin + glyphWidth;
                            g.DrawString(series.Title, AxisLabelFont, brush, (int)Math.Round(textX), (int)Math.Round(lineY), LegendFormat);
                            n++;
                        }
                    }
                }

                
            }
        }

        /// <summary>
        /// Re-draws the chart after updates to data.
        /// </summary>
        public void Draw()
        {
            if(XAxis.AutoScale)
            {
                double? xmin = DataSeries.Min(series => series.Points.Where(p => !p.IsNan()).Min(v => (double?)v.X));
                double? xmax = DataSeries.Max(series => series.Points.Where(p => !p.IsNan()).Max(v => (double?)v.X));
                if (xmin == null || xmax == null || xmax.Value <= xmin.Value)
                {
                    xmin = 0;
                    xmax = 1;
                }

                XAxis.Bounds = new Vec2(xmin.Value, xmax.Value);
            }

            if(YAxis.AutoScale)
            {
                double? ymin = DataSeries.Min(series => series.Points.Where(p => !p.IsNan()).Where(p => p.X >= XAxis.Bounds.X && p.X <= XAxis.Bounds.Y).Min(v => (double?)v.Y));
                double? ymax = DataSeries.Max(series => series.Points.Where(p => !p.IsNan()).Where(p => p.X >= XAxis.Bounds.X && p.X <= XAxis.Bounds.Y).Max(v => (double?)v.Y));

                if (ymin == null || ymax == null || ymax.Value <= ymin.Value)
                {
                    ymin = 0;
                    ymax = 1;
                }

                YAxis.Bounds = new Vec2(ymin.Value, ymax.Value);
            }

            // Create a back-buffer
            Image img = BackgroundImage;
            if(img == null || img.Size != ClientSize)
            {
                if (img != null)
                    img.Dispose();

                img = new Bitmap(Math.Max(1, ClientSize.Width), Math.Max(1, ClientSize.Height));
                BackgroundImage = img;
            }

            using (Graphics g = Graphics.FromImage(img))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;

                g.Clear(BackColor);

                // Calculate height of X-axis
                double xHeight = MeasureAxisHeight(g, XAxis, true, null);

                int dataTop = (int)Math.Round(AxisMargin);
                int dataBottom = (int)Math.Round(img.Height - xHeight);

                // Find Y ticks that fit into the image
                List<double> yTicks = FindTickLocations(g, YAxis, new Vec2(dataTop, dataBottom), false);

                // Calculate Y-axis width
                double yWidth = MeasureAxisHeight(g, YAxis, false, yTicks);

                // Find X ticks that fit into the image
                int dataLeft = (int)Math.Round(yWidth);
                int dataRight = (int)Math.Round(img.Width - AxisMargin);
                List<double> xTicks = FindTickLocations(g, XAxis, new Vec2(dataLeft, dataRight), true);
                
                Rectangle dataRect = Rectangle.FromLTRB(dataLeft, dataTop, dataRight, dataBottom);
                
                DrawData(g, dataRect);

                DrawAxes(g, dataRect, xTicks, yTicks);

                g.DrawRectangle(AxisPen, dataRect);

                DrawLegend(g, dataRect);
            }

            Invalidate();
        }

        /// <summary>
        /// Updates the control after resizing.
        /// </summary>
        /// <param name="e"></param>
        protected override void OnSizeChanged(EventArgs e)
        {
            base.OnSizeChanged(e);
            Draw();
        }

        /// <summary> 
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                AxisLabelFont.Dispose();
                AxisLabelFont = null;

                TickLabelFont.Dispose();
                TickLabelFont = null;

                AxisLabelFormat.Dispose();
                AxisLabelFormat = null;
                
                YTickFormat.Dispose();
                YTickFormat = null;

                BackgroundImage.Dispose();
                BackgroundImage = null;
            }
            base.Dispose(disposing);
        }
    }

    /// <summary>
    /// Enumerates possible locations for, e.g. legend box.
    /// </summary>
    public enum LegendLocation
    {
        None,
        TopRight,
        BottomRight,
        TopLeft,
        BottomLeft
    }

    /// <summary>
    /// Represents one axis of a chart.
    /// </summary>
    public class Axis
    {
        /// <summary>
        /// Minimum and maximum value shown on the chart.
        /// </summary>
        public Vec2 Bounds
        {
            get;
            set;
        }

        /// <summary>
        /// Axis label.
        /// </summary>
        public string Label
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets a value indicating whether the axis direction should be flipped.
        /// </summary>
        public bool Flipped
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets a value whether Bounds should be set automatically.
        /// </summary>
        public bool AutoScale
        {
            get;
            set;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <param name="label">Axis label.</param>
        public Axis(double min, double max, string label = "")
        {
            Bounds = new Vec2(min, max);
            Label = label;
        }

        /// <summary>
        /// Converts from data coordinates to image coordinates.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="imageLeft"></param>
        /// <param name="imageRight"></param>
        /// <returns></returns>
        internal double DataToImage(double x, double imageLeft, double imageRight)
        {
            if(!Flipped)
                return imageLeft + (x - Bounds.X) / (Bounds.Y - Bounds.X) * (imageRight - imageLeft);

            return imageRight - (x - Bounds.X) / (Bounds.Y - Bounds.X) * (imageRight - imageLeft);
        }
    }

    /// <summary>
    /// Single data series to be shown in the chart control.
    /// </summary>
    public class DataSeries
    {
        private List<Vec2> points;

        /// <summary>
        /// Gets and sets a list of points in this data series.
        /// </summary>
        public List<Vec2> Points
        {
            get
            {
                return points;
            }
            set
            {
                points = value;
                if (points == null)
                    points = new List<Vec2>();
            }
        }

        /// <summary>
        /// Pen that is used to draw the line corresponding to this series.
        /// </summary>
        public Pen LinePen
        {
            get;
            set;
        }

        /// <summary>
        /// Title of the series shown in legend.
        /// Set to null to not show this plot in the legend.
        /// </summary>
        public string Title
        {
            get;
            set;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        public DataSeries(Color color, string title = null)
        {
            LinePen = new Pen(color);
            Points = null;
            Title = title;
        }

        /// <summary>
        /// Draws this data series to given graphics.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="imageBounds">Rectangle where the data is to be drawn, in image coordinates.</param>
        /// <param name="xAxis">X-axis object.</param>
        /// <param name="yAxis">Y-axis object.</param>
        internal void Draw(Graphics g, Rectangle imageBounds, Axis xAxis, Axis yAxis)
        {
            g.SetClip(imageBounds);
            for(int n = 0; n < Points.Count - 1; n++)
            {
                Vec2 p0 = DataToImage(Points[n], imageBounds, xAxis, yAxis);
                Vec2 p1 = DataToImage(Points[n + 1], imageBounds, xAxis, yAxis);
                if(!p0.IsNan() && !p1.IsNan())
                    g.DrawLine(LinePen, p0.ToPointF(), p1.ToPointF());
            }
            g.ResetClip();
        }

        /// <summary>
        /// Converts a point from data coordinates to image coordinates.
        /// </summary>
        /// <param name="data">Point in data coordinates.</param>
        /// <param name="imageBounds">Bounds of data area in image coordinates.</param>
        /// <param name="xAxis">X-axis object.</param>
        /// <param name="yAxis">Y-axis object.</param>
        /// <returns>Point in image coordinates.</returns>
        private Vec2 DataToImage(Vec2 data, Rectangle imageBounds, Axis xAxis, Axis yAxis)
        {
            return new Vec2(xAxis.DataToImage(data.X, imageBounds.Left, imageBounds.Right),
                            yAxis.DataToImage(data.Y, imageBounds.Top, imageBounds.Bottom));
        }
    }
}
