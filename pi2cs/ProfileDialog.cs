using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace pi2cs
{
    /// <summary>
    /// Form that shows profile plot.
    /// </summary>
    public partial class ProfileDialog : Form
    {
        private Chart chart;

        /// <summary>
        /// Gets radius of the profile.
        /// </summary>
        public float ProfileWidth
        {
            get
            {
                return (float)numericUpDownProfileWidth.Value;
            }
        }

        /// <summary>
        /// Gets viewer whose profile is being plotted.
        /// </summary>
        public Pi2PictureViewer Viewer { get; }

        /// <summary>
        /// Clears the plot.
        /// </summary>
        public void Clear()
        {
            chart.DataSeries[0].Points.Clear();
        }

        /// <summary>
        /// Adds a point to the plot.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public void AddPoint(double x, double y)
        {
            chart.DataSeries[0].Points.Add(new Vec2(x, y));
        }

        /// <summary>
        /// Finishes update of the chart and draws it.
        /// </summary>
        public void Draw()
        {
            chart.Draw();
            chart.Refresh();
        }

        /// <summary>
        /// Cconstructor.
        /// </summary>
        /// <param name="viewer"></param>
        public ProfileDialog(Pi2PictureViewer viewer)
        {
            chart = new Chart();
            chart.Dock = DockStyle.Fill;
            chart.XAxis.Label = "Position [pix]";
            chart.XAxis.AutoScale = true;
            chart.YAxis.Label = "Gray value";
            chart.YAxis.AutoScale = true;
            chart.NewSeries();
            chart.LegendLocation = LegendLocation.None;
            Controls.Add(chart);

            Viewer = viewer;
            InitializeComponent();
        }

        private void ProfileDialog_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (e.CloseReason == CloseReason.UserClosing)
            {
                Viewer.EraseProfileAnnotation();
                e.Cancel = true;
                Hide();
            }
        }

        private void numericUpDownProfileWidth_ValueChanged(object sender, EventArgs e)
        {
            if (Viewer.ProfileAnnotation != null)
            {
                Viewer.UpdateProfile();
                Viewer.Refresh();
            }
        }
    }
}
