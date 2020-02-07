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
    /// Form that shows histogram chart.
    /// </summary>
    public partial class HistogramDialog : Form
    {
        private HistogramBox Histogram;

        /// <summary>
        /// Gets and sets the box whose histogram is shown.
        /// </summary>
        public Pi2PictureBox PictureBox
        {
            get
            {
                return Histogram.PictureBox;
            }
            set
            {
                Histogram.PictureBox = value;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="picture">Picture box whose histogram is to be shown.</param>
        public HistogramDialog(Pi2PictureBox picture)
        {
            InitializeComponent();

            Histogram = new HistogramBox(picture);
            Histogram.Dock = DockStyle.Fill;
            Controls.Add(Histogram);
        }

        private void HistogramDialog_Deactivate(object sender, EventArgs e)
        {
            Hide();
        }

        private void HistogramDialog_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (e.CloseReason == CloseReason.UserClosing)
            {
                Hide();
                e.Cancel = true;
            }
        }

        private void HistogramDialog_VisibleChanged(object sender, EventArgs e)
        {
            if(Visible)
                Histogram.UpdateHistogram();
        }
    }
}
