using pi2cs;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace pi2csWinFormsTest
{
    public partial class PictureBoxTestForm : Form
    {
        private Pi2 pi;

        private Pi2PictureBox box;
        private HistogramBox hist;

        public PictureBoxTestForm()
        {
            InitializeComponent();
            pi = new Pi2();

            Pi2Image img = pi.NewImage(ImageDataType.UInt16);
            pi.Read(img, "../../testing/uint16.png");

            box = new Pi2PictureBox(img, 2, 0, 65535/2);
            box.AutoResize = true;
            box.Location = new Point(50, 200);
            Controls.Add(box);

            hist = new HistogramBox();
            hist.Location = new Point(300, 30);
            hist.PictureBox = box;
            Controls.Add(hist);

            box.UpdateImage();
            hist.UpdateHistogram();
        }

        private void buttonZoomIn_Click(object sender, EventArgs e)
        {
            box.Zoom *= 1.25f;
            box.UpdateScreen();
        }

        private void buttonZoomOut_Click(object sender, EventArgs e)
        {
            box.Zoom /= 1.25f;
            if (box.Zoom <= 0.01f)
                box.Zoom = 0.01f;
            box.UpdateScreen();
        }

        private void buttonTranslate_Click(object sender, EventArgs e)
        {
            box.PicturePosition = new PointF(box.PicturePosition.X + 1, box.PicturePosition.Y);
            box.UpdateScreen();
        }


        //private void button1_Click(object sender, EventArgs e)
        //{

        //    IntPtr pi = PiLib.CreatePI();

        //    bool retval = PiLib.Run(pi, "distribute(false)");
        //    labelRetVal.Text = "Return value: " + retval;

        //    int line = PiLib.LastErrorLine(pi);
        //    labelErrorMessage.Text = PiLib.LastErrorMessage(pi);

        //    labelCommandList.Text = PiLib.CommandList(pi);

        //    PiLib.DestroyPI(pi);

        //}
    }
}
