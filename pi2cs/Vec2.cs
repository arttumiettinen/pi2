using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    public class Vec2
    {
        public double X;
        public double Y;

        public Vec2()
        {
            X = 0;
            Y = 0;
        }

        public Vec2(double x, double y)
        {
            X = x;
            Y = y;
        }

        public Vec2(Point p)
        {
            X = p.X;
            Y = p.Y;
        }

        public Vec2(PointF p)
        {
            X = p.X;
            Y = p.Y;
        }

        public Vec2(Size p)
        {
            X = p.Width;
            Y = p.Height;
        }

        public Vec2(SizeF p)
        {
            X = p.Width;
            Y = p.Height;
        }

        public PointF ToPointF()
        {
            return new PointF((float)X, (float)Y);
        }

        public static Vec2 operator+(Vec2 l, Vec2 r)
        {
            return new Vec2(l.X + r.X, l.Y + r.Y);
        }

        public static Vec2 operator -(Vec2 l, Vec2 r)
        {
            return new Vec2(l.X - r.X, l.Y - r.Y);
        }

        public static Vec2 operator *(Vec2 l, double r)
        {
            return new Vec2(l.X * r, l.Y * r);
        }

        public static Vec2 operator *(double r, Vec2 l)
        {
            return new Vec2(l.X * r, l.Y * r);
        }

        public static Vec2 operator /(Vec2 l, double r)
        {
            return new Vec2(l.X / r, l.Y / r);
        }

        public double Dot(Vec2 r)
        {
            return X * r.X + Y * r.Y;
        }

        public double NormSquared()
        {
            return this.Dot(this);
        }

        public double Norm()
        {
            return Math.Sqrt(NormSquared());
        }

        public override string ToString()
        {
            return "[" + X.ToString(CultureInfo.InvariantCulture) + ", " + Y.ToString(CultureInfo.InvariantCulture) + "]";
        }
    }
}
