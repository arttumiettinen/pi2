using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using pi2cs;

namespace pi2csWinFormsTest
{
    static class Program
    {
        private static void MetadataTest()
        {
            Pi2 pi = new Pi2();
            Pi2Image img = pi.NewImage(ImageDataType.UInt16);

            pi.SetMeta(img, "Key1", "String value 1");
            pi.SetMeta(img, "Key2", "String value 2");
            string[] keys = pi.ListMeta(img);

            foreach (string s in keys)
                Console.WriteLine(s);
        }

        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);

            MetadataTest();
            //Application.Run(new PictureBoxTestForm());
            //Application.Run(new PictureViewerTestForm());
            //Application.Run(new ChartTestForm());
        }
    }
}

