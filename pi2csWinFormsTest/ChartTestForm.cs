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
    public partial class ChartTestForm : Form
    {
        private Chart chart;

        public ChartTestForm()
        {
            InitializeComponent();

            chart = new Chart();
            chart.Dock = DockStyle.Fill;
            Controls.Add(chart);

            DataSeries series = chart.NewSeries("Parabola");
            for (double x = -10; x <= 10; x += 0.5)
                series.Points.Add(new Vec2(x, x * x));

            series = chart.NewSeries("Hyperbola");
            for (double x = -10; x <= 10; x += 0.5)
                series.Points.Add(new Vec2(x, x * x * x));

            chart.XAxis.Bounds = new Vec2(-11.3423, 13.212);
            //chart.XAxis.Bounds = new Vec2(-0.0001, -0.00005);
            //chart.XAxis.Bounds = new Vec2(-10000, 5);
            chart.YAxis.Bounds = new Vec2(-130, 130);
            //chart.YAxis.Bounds = new Vec2(-1300, -1000);

            chart.LegendLocation = LegendLocation.BottomLeft;

            chart.Draw();
        }
    }
}
