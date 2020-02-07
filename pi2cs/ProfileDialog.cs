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
            chart.Series[0].Points.Clear();
            chart.ChartAreas[0].AxisX.Minimum = 0;
            chart.ChartAreas[0].AxisX.Maximum = 1;
        }

        /// <summary>
        /// Adds a point to the plot.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public void AddPoint(double x, double y)
        {
            chart.Series[0].Points.AddXY(x, y);
            chart.ChartAreas[0].AxisX.Maximum = Math.Max(chart.ChartAreas[0].AxisX.Maximum, x);
        }

        /// <summary>
        /// Cconstructor.
        /// </summary>
        /// <param name="viewer"></param>
        public ProfileDialog(Pi2PictureViewer viewer)
        {
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
