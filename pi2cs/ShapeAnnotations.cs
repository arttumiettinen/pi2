using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    public sealed class RectangleAnnotation : Annotation2D
    {
        public override int RequiredPointCount
        {
            get
            {
                return 2;
            }
        }

        [DisplayName("Draw center cross")]
        public bool DrawCenterCross
        {
            get;
            set;
        } = true;

        [DisplayName("Dash style for center cross")]
        public DashStyle CenterCrossStyle
        {
            get;
            set;
        } = DashStyle.Dash;

        protected override void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics)
        {
            if (ControlPoints.Count >= 2)
            {
                Vec2 p1 = box.PictureToScreen(ControlPoints[0].GetXY());
                Vec2 p2 = box.PictureToScreen(ControlPoints[1].GetXY());

                float x = (float)Math.Min(p1.X, p2.X);
                float y = (float)Math.Min(p1.Y, p2.Y);
                float w = (float)Math.Abs(p2.X - p1.X);
                float h = (float)Math.Abs(p2.Y - p1.Y);

                using (Pen pen = new Pen(Color, (float)LineWidth))
                {
                    g.DrawRectangle(pen, x, y, w, h);

                    if (DrawCenterCross)
                    {
                        float s = 20;
                        float cx = x + w / 2.0f;
                        float cy = y + h / 2.0f;
                        pen.DashStyle = CenterCrossStyle;
                        g.DrawLine(pen, cx, y - s, cx, y + h + s);
                        g.DrawLine(pen, x - s, cy, x + w + s, cy);
                    }
                }
                 
                

                using (Pen pen = box.GetPickPen(this))
                    pickGraphics.DrawRectangle(pen, x, y, w, h);
            }
        }
    }
}
