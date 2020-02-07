using pi2cs;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace pi2csWinFormsTest
{
    public partial class PictureViewerTestForm : Form
    {
        private Pi2 pi;

        public PictureViewerTestForm()
        {
            InitializeComponent();

            pi = new Pi2();
            //Pi2Image img = pi.Read("../../testing/uint8.png");
            //Pi2Image img = pi.Read("../../testing/t1-head_slice64.tif");
            //Pi2Image img = pi.Read("../../testing/t1-head");
            //pi.Convert(img, ImageDataType.Float32);
            //pi.Divide(img, 500);
            //pi.Set(img, 10, 10, 0, float.PositiveInfinity);

            Pi2Image img = pi.MapRaw("../../testing/input_data/t1-head_256x256x129.raw");

            Pi2PictureViewer viewer = new Pi2PictureViewer();
            viewer.PiImage = img;
            viewer.Dock = DockStyle.Fill;
            Controls.Add(viewer);

            Pi2PictureToolStrip tools = new Pi2PictureToolStrip();
            tools.Dock = DockStyle.Top;
            Controls.Add(tools);
            tools.PictureViewer = viewer;

            viewer.UpdateUI(true);
        }
    }
}
