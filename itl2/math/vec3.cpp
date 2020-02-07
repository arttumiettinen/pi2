
#include "math/vec3.h"
#include "utilities.h"
#include "test.h"
#include <iostream>

using namespace std;

namespace itl2
{
	namespace tests
	{
		void vectorAngles()
		{
			Vec3d a(1, 0, 0);

			itl2::testAssert(NumberUtils<double>::equals(a.angleTo(Vec3d(0.5, 0.5, 0)), PI / 4, 1e-10), "45 deg angle");
			itl2::testAssert(NumberUtils<double>::equals(a.sharpAngleTo(Vec3d(0.5, 0.5, 0)), PI / 4, 1e-10), "45 deg angle sharp");
			itl2::testAssert(NumberUtils<double>::equals(a.angleTo(Vec3d(-0.5, 0.5, 0)), 3 * PI / 4, 1e-10), "135/45 deg angle");
			itl2::testAssert(NumberUtils<double>::equals(a.sharpAngleTo(Vec3d(-0.5, 0.5, 0)), PI / 4, 1e-10), "135/45 deg angle sharp");

			for (double gtAngle = 0.0; gtAngle < PI; gtAngle += 0.1)
			{
				Vec3d b = a.rotate(Vec3d(0, 0, 1), gtAngle);

				double angle2;
				if (a.equals(b))
					angle2 = 0;
				else
				{
					Vec3d V = a.cross(b);
					double sine = a.cross(b).dot(V) / (a.norm() * b.norm() * V.norm());
					if (sine >= 1)
						angle2 = PI / 2;
					else if (sine <= -1)
						angle2 = -PI / 2;
					else
						angle2 = asin(sine);
				}

				cout << gtAngle << ", " << a.angleTo(b) << ", " << a.sharpAngleTo(b) << ", " << angle2 << endl;

				itl2::testAssert(NumberUtils<double>::equals(a.angleTo(b), gtAngle, 1e-10), string("angle between a and b at theta = ") + itl2::toString(gtAngle));
				itl2::testAssert(NumberUtils<double>::equals(a.sharpAngleTo(b), angle2, 1e-10), string("sharp angle between a and b at theta = ") + itl2::toString(gtAngle));
			}

		}
	}
}