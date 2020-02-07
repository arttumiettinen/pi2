using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;

namespace pi2cs
{
    /// <summary>
    /// Base class for annotations
    /// </summary>
    public abstract class Annotation
    {
        /// <summary>
        /// Gets and sets color used for drawing.
        /// </summary>
        [DisplayName("Color")]
        public Color Color
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets line width.
        /// </summary>
        [DisplayName("Line width")]
        public double LineWidth
        {
            get;
            set;
        }

        // <summary>
        // Gets and sets a value indicating whether the object is higlighted or not.
        // </summary>
        //[Browsable(false)]
        //public bool IsHighlighted
        //{
        //    get;
        //    set;
        //}
        
        /// <summary>
        /// Gets and sets a modifiable list of control points of this drawing object.
        /// </summary>
        [Browsable(false)]
        public List<Vec3> ControlPoints
        {
            get;
            set;
        }

        /// <summary>
        /// Gets count of control points required for this object to be ready.
        /// </summary>
        [Browsable(false)]
        public abstract int RequiredPointCount
        {
            get;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        protected Annotation()
        {
            ControlPoints = new List<Vec3>();
            Color = Color.Yellow;
            LineWidth = 1;
        }

        /// <summary>
        /// Gets line width adjusted for highlighting.
        /// </summary>
        /// <returns></returns>
        //protected double GetLineWidthAdjustedForHighlighting()
        //{
        //    if (!IsHighlighted)
        //        return LineWidth;
        //    else
        //        return 2 * LineWidth;
        //}

        
        protected void DrawControlPoint(Pi2PictureBox box, Graphics g, Graphics pickGraphics, int pi, int currZ)
        {
            Vec3 p = ControlPoints[pi];
            if (Math.Abs(p.Z - currZ) < 0.5f)
            {
                const float R = 3;
                Vec2 pos = box.PictureToScreen(p.GetXY());
                RectangleF rect = new RectangleF((float)(pos.X - R), (float)(pos.Y - R), 2*R, 2*R);
                g.FillRectangle(Brushes.WhiteSmoke, rect);
                g.DrawRectangle(Pens.Black, rect.X, rect.Y, rect.Width, rect.Height);

                using (Brush b = box.GetPickBrush(this, pi))
                    pickGraphics.FillRectangle(b, rect.X, rect.Y, rect.Width, rect.Height);
            }
        }

        /// <summary>
        /// Draws control points if DrawControlPoints is set.
        /// </summary>
        protected void DoDrawControlPoints(Pi2PictureBox box, Graphics g, Graphics pickGraphics, int currZ)
        {
            for(int n = 0; n < ControlPoints.Count; n++)
            {
                DrawControlPoint(box, g, pickGraphics, n, currZ);
            }
        }

        /// <summary>
        /// Draw this object to the given renderer.
        /// This version of Draw method is called if no data is required by the drawing object.
        /// The default implementation draws control points if appropriate.
        /// </summary>
        public virtual void Draw(Pi2PictureBox box, Graphics g, Graphics pickGraphics, int currZ)
        {
            DoDrawControlPoints(box, g, pickGraphics, currZ);
        }


        ///// <summary>
        ///// Get end point for line that should look like infinite.
        ///// </summary>
        ///// <param name="start"></param>
        ///// <param name="dir"></param>
        ///// <param name="screen"></param>
        ///// <returns></returns>
        //protected static Vec2 GetInfiniteLineEndPoint(Vec2 start, Vec2 dir, RectangleD screen)
        //{
        //    const double EPS = 0.001;
        //    if (dir.X > EPS && start.X < screen.Right)
        //    {
        //        double c = (screen.Right - start.X) / dir.X;
        //        return start + c * dir;
        //    }
        //    else if (dir.X < -EPS && start.X >= screen.Left)
        //    {
        //        double c = (screen.Left - start.X) / dir.X;
        //        return start + c * dir;
        //    }
        //    else if (dir.Y > EPS && start.Y < screen.Bottom)
        //    {
        //        double c = (screen.Bottom - start.Y) / dir.Y;
        //        return start + c * dir;
        //    }
        //    else if (dir.Y < -EPS && start.Y >= screen.Top)
        //    {
        //        double c = (screen.Top - start.Y) / dir.Y;
        //        return start + c * dir;
        //    }

        //    // dir is zero vector or the line is not visible. Return start point.
        //    return start;
        //}

    }

    /// <summary>
    /// 2D annotation.
    /// Provides possibility to draw to designated slice or to all slices.
    /// </summary>
    public abstract class Annotation2D : Annotation
    {
        [DisplayName("Draw to all slices")]
        public bool DrawToAllSlices
        {
            get;
            set;
        }

        public Annotation2D()
        {
            DrawToAllSlices = true;
        }

        /// <summary>
        /// Method that draws the 2D annotation, discarding all Z coordinates of control points.
        /// </summary>
        /// <param name="box">Picture box that owns the annotation.</param>
        /// <param name="g">Graphics where to draw to.</param>
        /// <param name="pickGraphics">Graphics where pick info is going to be drawn to.</param>
        protected abstract void Draw2D(Pi2PictureBox box, Graphics g, Graphics pickGraphics);

        public override void Draw(Pi2PictureBox box, Graphics g, Graphics pickGraphics, int currZ)
        {
            if ((ControlPoints.Count > 0 && (Math.Abs(ControlPoints[0].Z - currZ) < 0.5 || DrawToAllSlices)) ||
                RequiredPointCount <= 0)
            {
                Draw2D(box, g, pickGraphics);
                for (int n = 0; n < ControlPoints.Count; n++)
                {
                    DrawControlPoint(box, g, pickGraphics, n, (int)Math.Round(ControlPoints[n].Z));
                }
            }
        }
    }
}
