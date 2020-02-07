using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;

namespace pi2cs
{
    /// <summary>
    /// Horizontal line
    /// </summary>
    public sealed class HLineAnnotation : Annotation2D
    {
        public override int RequiredPointCount
        {
            get
            {
                return 1;
            }
        }

        protected override void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics)
        {
            Vec2 p1 = box.PictureToScreen(ControlPoints[0].GetXY());
            p1.X = 0;
            Vec2 p2 = box.PictureToScreen(ControlPoints[0].GetXY());
            p2.X = box.Width;

            using (Pen pen = new Pen(Color, (float)LineWidth))
                g.DrawLine(pen, p1.ToPointF(), p2.ToPointF());

            using (Pen pen = box.GetPickPen(this))
                pickGraphics.DrawLine(pen, p1.ToPointF(), p2.ToPointF());
        }
    }

    /// <summary>
    /// Vertical line
    /// </summary>
    public sealed class VLineAnnotation : Annotation2D
    {
        public override int RequiredPointCount
        {
            get
            {
                return 1;
            }
        }

        protected override void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics)
        {
            Vec2 p1 = box.PictureToScreen(ControlPoints[0].GetXY());
            p1.Y = 0;
            Vec2 p2 = box.PictureToScreen(ControlPoints[0].GetXY());
            p2.Y = box.Height;

            using (Pen pen = new Pen(Color, (float)LineWidth))
                g.DrawLine(pen, p1.ToPointF(), p2.ToPointF());

            using (Pen pen = box.GetPickPen(this))
                pickGraphics.DrawLine(pen, p1.ToPointF(), p2.ToPointF());
        }
    }

    /// <summary>
    /// Line segment
    /// </summary>
    public sealed class LineSegmentAnnotation : Annotation2D
    {
        public override int RequiredPointCount
        {
            get
            {
                return 2;
            }
        }

        protected override void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics)
        {
            if (ControlPoints.Count >= 2)
            {
                Vec2 p1 = box.PictureToScreen(ControlPoints[0].GetXY());
                Vec2 p2 = box.PictureToScreen(ControlPoints[1].GetXY());

                using (Pen pen = new Pen(Color, (float)LineWidth))
                    g.DrawLine(pen, p1.ToPointF(), p2.ToPointF());

                using (Pen pen = box.GetPickPen(this))
                    pickGraphics.DrawLine(pen, p1.ToPointF(), p2.ToPointF());
            }
        }
    }

    /// <summary>
    /// Center cross annotation.
    /// </summary>
    public sealed class CenterCrossAnnotation : Annotation2D
    {
        public override int RequiredPointCount
        {
            get
            {
                return 0;
            }
        }

        protected override void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics)
        {
            Vec2 p1 = box.PictureToScreen(new Vec2(box.OriginalWidth / 2.0, 0));
            Vec2 p2 = box.PictureToScreen(new Vec2(box.OriginalWidth / 2.0, box.OriginalHeight));

            Vec2 p3 = box.PictureToScreen(new Vec2(0, box.OriginalHeight / 2.0));
            Vec2 p4 = box.PictureToScreen(new Vec2(box.OriginalWidth, box.OriginalHeight / 2.0));

            using (Pen pen = new Pen(Color, (float)LineWidth))
            {
                g.DrawLine(pen, p1.ToPointF(), p2.ToPointF());
                g.DrawLine(pen, p3.ToPointF(), p4.ToPointF());
            }

            using (Pen pen = box.GetPickPen(this))
            {
                pickGraphics.DrawLine(pen, p1.ToPointF(), p2.ToPointF());
                pickGraphics.DrawLine(pen, p3.ToPointF(), p4.ToPointF());
            }
        }
    }
}