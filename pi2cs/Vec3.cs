using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    /// <summary>
    /// 3-component vector
    /// </summary>
    public class Vec3
    {
        public double X;
        public double Y;
        public double Z;

        public Vec3()
        {
            X = 0;
            Y = 0;
            Z = 0;
        }

        public Vec3(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public Vec3(Vec2 xy, double z)
        {
            X = xy.X;
            Y = xy.Y;
            Z = z;
        }

        public Vec2 GetXY()
        {
            return new Vec2(X, Y);
        }

        public static Vec3 operator +(Vec3 l, Vec3 r)
        {
            return new Vec3(l.X + r.X, l.Y + r.Y, l.Z + r.Z);
        }

        public static Vec3 operator -(Vec3 l, Vec3 r)
        {
            return new Vec3(l.X - r.X, l.Y - r.Y, l.Z - r.Z);
        }

        public override string ToString()
        {
            return "[" + X.ToString(CultureInfo.InvariantCulture) + ", " + Y.ToString(CultureInfo.InvariantCulture) + ", " + Z.ToString(CultureInfo.InvariantCulture) + "]";
        }
    }
}
