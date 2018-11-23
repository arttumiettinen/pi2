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
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {

            IntPtr pi = PiLib.CreatePI();

            bool retval = PiLib.Run(pi, "distribute(false)");
            labelRetVal.Text = "Return value: " + retval;

            int line = PiLib.LastErrorLine(pi);
            labelErrorMessage.Text = PiLib.LastErrorMessage(pi);

            labelCommandList.Text = PiLib.CommandList(pi);

            PiLib.DestroyPI(pi);

        }
    }
}
