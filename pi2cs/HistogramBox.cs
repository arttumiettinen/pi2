﻿using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using System.Globalization;

namespace pi2cs
{
    /// <summary>
    /// Control that shows histogram of Pi2 image and scroll bars to adjust gray scale.
    /// </summary>
    public partial class HistogramBox : UserControl
    {
        private Pi2PictureBox pb;

        private Chart chartHist;

        /// <summary>
        /// Gets and sets the box whose histogram is shown.
        /// </summary>
        public Pi2PictureBox PictureBox
        {
            get
            {
                return pb;
            }
            set
            {
                pb = value;
                UpdateHistogram();
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        public HistogramBox() : this(null)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="pictureBox"></param>
        public HistogramBox(Pi2PictureBox pictureBox)
        {
            InitializeComponent();
            chartHist = new Chart();
            chartHist.NewSeries();
            chartHist.NewSeries();
            chartHist.NewSeries();
            chartHist.XAxis.Label = "Gray value";
            chartHist.YAxis.Label = "Count";
            chartHist.Dock = DockStyle.Fill;
            Controls.Add(chartHist);
            chartHist.BringToFront();
            
            PictureBox = pictureBox;
        }

        /// <summary>
        /// This is from StackOverflow, see
        /// https://stackoverflow.com/questions/8137391/percentile-calculation
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="excelPercentile"></param>
        /// <returns></returns>
        private static double Percentile(double[] sequence, double excelPercentile)
        {
            Array.Sort(sequence);
            int N = sequence.Length;
            double n = (N - 1) * excelPercentile + 1;
            // Another method: double n = (N + 1) * excelPercentile;
            if (n == 1d)
            {
                return sequence[0];
            }
            else if (n == N)
            {
                return sequence[N - 1];
            }
            else
            {
                int k = (int)n;
                double d = n - k;
                return sequence[k - 1] + d * (sequence[k] - sequence[k - 1]);
            }
        }

        /// <summary>
        /// Range mapped to min and max of scroll bars.
        /// </summary>
        float rangeMin;
        float rangeMax;

        private bool updating = false;

        /// <summary>
        /// Updates histogram window.
        /// </summary>
        public void UpdateHistogram()
        {
            if (PictureBox != null)
            {
                updating = true;
                try
                {
                    IReadOnlyList<PointF> histogram = PictureBox.Histogram;

                    while (chartHist.DataSeries[1].Points.Count < 2)
                        chartHist.DataSeries[1].Points.Add(new Vec2());
                    while (chartHist.DataSeries[2].Points.Count < 2)
                        chartHist.DataSeries[2].Points.Add(new Vec2());

                    var points = chartHist.DataSeries[0].Points;
                    while (points.Count < histogram.Count)
                        points.Add(new Vec2(0, 0));

                    List<double> counts = new List<double>();

                    double max = 0;
                    for (int n = 0; n < histogram.Count; n++)
                    {
                        points[n] = new Vec2(histogram[n]);
                        counts.Add(histogram[n].Y);
                        if (histogram[n].Y > max)
                            max = histogram[n].Y;
                    }
                    
                    chartHist.YAxis.Bounds = new Vec2(0, Percentile(counts.ToArray(), 0.995));

                    if (checkFullRange.Checked)
                    {
                        // Range from DynamicMin to DynamicMax
                        rangeMin = PictureBox.DynamicMin;
                        rangeMax = PictureBox.DynamicMax;
                    }
                    else
                    {
                        // Range from GlobalMin to GlobalMax
                        rangeMin = PictureBox.GlobalMin;
                        rangeMax = PictureBox.GlobalMax;
                    }

                    //float minVal = PictureBox.Min;
                    //float maxVal = PictureBox.Max;

                    //rangeMin = Math.Min(rangeMin, minVal);
                    //rangeMax = Math.Max(rangeMax, maxVal);

                    if (rangeMin >= rangeMax)
                        rangeMax = rangeMin + 1;
                    
                    chartHist.XAxis.Bounds = new Vec2(rangeMin, rangeMax);

                    chartHist.DataSeries[1].Points[0] = new Vec2(PictureBox.Min, chartHist.YAxis.Bounds.X);
                    chartHist.DataSeries[1].Points[1] = new Vec2(PictureBox.Min, chartHist.YAxis.Bounds.Y);

                    chartHist.DataSeries[2].Points[0] = new Vec2(PictureBox.Max, chartHist.YAxis.Bounds.X);
                    chartHist.DataSeries[2].Points[1] = new Vec2(PictureBox.Max, chartHist.YAxis.Bounds.Y);

                    scrollMinimum.Maximum = Math.Max(256, (int)(rangeMax - rangeMin)) + scrollMinimum.LargeChange;
                    scrollMaximum.Maximum = scrollMinimum.Maximum;
                    
                    SetScrollBarValue(scrollMinimum, PictureBox.Min);
                    SetScrollBarValue(scrollMaximum, PictureBox.Max);
                    UpdateLabels();
                    chartHist.Draw();
                    chartHist.Invalidate();
                }
                finally
                {
                    updating = false;
                }
            }
        }

        private void checkFullRange_CheckedChanged(object sender, EventArgs e)
        {
            if (!updating)
                UpdateHistogram();
        }

        /// <summary>
        /// Gets value of given scroll bar as gray value.
        /// </summary>
        /// <param name="bar"></param>
        /// <returns></returns>
        private float GetScrollBarValue(ScrollBar bar)
        {
            return rangeMin + (float)bar.Value / (bar.Maximum - bar.LargeChange) * (rangeMax - rangeMin);
        }

        /// <summary>
        /// Sets value of scroll bar based on gray value.
        /// </summary>
        /// <param name="bar"></param>
        /// <param name="value"></param>
        private void SetScrollBarValue(ScrollBar bar, float value)
        {
            int sval = (int)Math.Round((value - rangeMin) / (rangeMax - rangeMin) * (bar.Maximum - bar.LargeChange));
            if (sval < bar.Minimum)
                sval = bar.Minimum;
            if (sval > bar.Maximum - bar.LargeChange)
                sval = bar.Maximum - bar.LargeChange;
            int oldVal = bar.Value;
            bar.Value = sval;
        }

        /// <summary>
        /// Updates min/max label text.
        /// </summary>
        private void UpdateLabels()
        {
            if (PictureBox.OriginalDataType == ImageDataType.Float32)
                labelMinimum.Text = PictureBox.Min.ToString("0.000", CultureInfo.InvariantCulture);
            else
                labelMinimum.Text = PictureBox.Min.ToString("0", CultureInfo.InvariantCulture);

            if (PictureBox.OriginalDataType == ImageDataType.Float32)
                labelMaximum.Text = PictureBox.Max.ToString("0.000", CultureInfo.InvariantCulture);
            else
                labelMaximum.Text = PictureBox.Max.ToString("0", CultureInfo.InvariantCulture);
        }

        private void scrollMinimum_ValueChanged(object sender, EventArgs e)
        {
            if (PictureBox != null && !updating)
            {
                //PictureBox.Min = scrollMinimum.Value;
                PictureBox.Min = GetScrollBarValue(scrollMinimum);

                //chartHist.ChartAreas[0].AxisX.Minimum = PictureBox.Min;
                PictureBox.UpdateScreen();

                UpdateLabels();

                chartHist.DataSeries[1].Points[0] = new Vec2(PictureBox.Min, chartHist.YAxis.Bounds.X);
                chartHist.DataSeries[1].Points[1] = new Vec2(PictureBox.Min, chartHist.YAxis.Bounds.Y);
                chartHist.Draw();
                chartHist.Refresh();
            }
        }

        private void scrollMaximum_ValueChanged(object sender, EventArgs e)
        {
            if (PictureBox != null && !updating)
            {
                //PictureBox.Max = scrollMaximum.Value;
                PictureBox.Max = GetScrollBarValue(scrollMaximum);

                //chartHist.ChartAreas[0].AxisX.Maximum = PictureBox.Max;
                PictureBox.UpdateScreen();

                UpdateLabels();

                chartHist.DataSeries[2].Points[0] = new Vec2(PictureBox.Max, chartHist.YAxis.Bounds.X);
                chartHist.DataSeries[2].Points[1] = new Vec2(PictureBox.Max, chartHist.YAxis.Bounds.Y);
                chartHist.Draw();
                chartHist.Refresh();
            }
        }
    }
}
