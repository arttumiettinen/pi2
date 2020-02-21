
#include "transform.h"
#include "io/raw.h"
#include "io/itltiff.h"
#include "math/vec3.h"
#include "pointprocess.h"
#include "testutils.h"
#include "generation.h"


using namespace std;

namespace itl2
{
	namespace tests
	{
		void scaleHelper(const Image<uint16_t>& head16, size_t size)
		{
			Image<uint16_t> scaled1(size, size, size/2);
			scale(head16, scaled1, false);
			tiff::writed(scaled1, "./transform/scale_no_avg_" + toString(size));

			Image<uint16_t> scaled2(size, size, size/2);
			scale(head16, scaled2, true);
			tiff::writed(scaled2, "./transform/scale_with_avg_" + toString(size));
		}

		void scale()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");

			scaleHelper(head16, 255);

			for (coord_t size = 50; size < 250; size += 50)
			{
				scaleHelper(head16, size);
			}


			scaleHelper(head16, 127);
			scaleHelper(head16, 129);
		}


		void translate()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");

			add(head16, 1);

			Image<uint16_t> shifted(head16.dimensions());

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), NearestNeighbourInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_nearest");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), LinearInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_linear");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), LinearInvalidValueInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero, 0, 0));
			raw::writed(shifted, "./transform/shift_linear_invalidvalue");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), CubicInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_cubic");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), CubicInvalidValueInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero, 0, 0));
			raw::writed(shifted, "./transform/shift_cubic_invalidvalue");
		}

		void binning()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");

			Image<uint16_t> headb;
			itl2::binning(head16, headb, 2);
			raw::writed(headb, "./transform/binning_2");

			itl2::binning<uint16_t, uint16_t, binningop::max<uint16_t, uint16_t> >(head16, headb, 2);
			raw::writed(headb, "./transform/binning_2_max");

			itl2::binning<uint16_t, uint16_t, binningop::min<uint16_t, uint16_t> >(head16, headb, 2);
			raw::writed(headb, "./transform/binning_2_min");

			itl2::binning(head16, headb, 3);
			raw::writed(headb, "./transform/binning_3");

			itl2::binning(head16, headb, 4);
			raw::writed(headb, "./transform/binning_4");
		}

		void genericTransform()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16, headb;
			raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");

			vector<Vec3f> refPoints, defPoints;
			defPoints.push_back(Vec3f(0, 0, 64));
			refPoints.push_back(Vec3f(0, 50, 64));

			
			defPoints.push_back(Vec3f(256, 256, 64));
			refPoints.push_back(Vec3f(256, 256-50, 64));
			
			headb.ensureSize(head16.dimensions());
			itl2::genericTransform(head16, headb, Vec3c(0, 0, 0), refPoints, defPoints);

			raw::writed(headb, "./transform/generic_transform");
		}

		void scaleLabels()
		{
			// NOTE: No asserts

			// Small test
			{
				Image<uint8_t> img(3, 3);
				img(0, 0) = 128;
				img(2, 2) = 255;

				Image<uint8_t> out;
				itl2::scaleLabels(img, out, true, Vec3c(10, 10, 10));
				raw::writed(out, "./scale_labels/small");
			}

			// Larger test
			{
				Image<uint8_t> img;
				raw::read(img, "./input_data/t1-head_bin");

				Image<uint8_t> out;
				itl2::scaleLabels(img, out, true, Vec3c(5, 5, 5));
				raw::writed(out, "./scale_labels/head");
			}
		}

		void rot90()
		{
			Image<uint8_t> orig(100, 200, 30);
			ramp(orig, 0);
			raw::writed(orig, "./rotations/orig");

			Image<uint8_t> rot;
			Image<uint8_t> gt(200, 100, 30);

			rot90cw(orig, rot);
			ramp(gt, 1);

			raw::writed(rot, "./rotations/rot90cw");
			raw::writed(gt, "./rotations/rot90cw_gt");

			checkDifference(rot, gt, "rotate 90 deg cw");

			rot90ccw(orig, rot);
			ramp(gt, 1);
			flip(gt, 1);

			raw::writed(rot, "./rotations/rot90ccw");
			raw::writed(gt, "./rotations/rot90ccw_gt");

			checkDifference(rot, gt, "rotate 90 deg ccw");
		}

		void rotate()
		{
			Image<uint16_t> orig;
			raw::read(orig, "./input_data/t1-head");

			Image<uint16_t> out(orig.dimensions());
			rotate(orig, out, degToRad(45));
			tiff::writed(out, "./rotations/head_z_+45");

			rotate(orig, out, degToRad(-90));
			tiff::writed(out, "./rotations/head_z_-90");

			rotate(orig, out, degToRad(45), Vec3d(1, 0, 0));
			tiff::writed(out, "./rotations/head_x_+45");

			rotate(orig, out, degToRad(45), Vec3d(0, 1, 0));
			tiff::writed(out, "./rotations/head_y_+45");

			rotate(orig, out, degToRad(45), Vec3d(1, 1, 0));
			tiff::writed(out, "./rotations/head_xy_+45");
		}

		void reslice()
		{
			Image<uint16_t> orig;

			// Create "coordinate system" image
			orig.ensureSize(100, 200, 300);
			draw(orig, AABox(Vec3c(0, 0, 0), Vec3c(orig.width(), 10, 10)), (uint16_t)100);
			draw(orig, AABox(Vec3c(0, 0, 0), Vec3c(10, orig.height(), 10)), (uint16_t)150);
			draw(orig, AABox(Vec3c(0, 0, 0), Vec3c(10, 10, orig.depth())), (uint16_t)200);

			tiff::writed(orig, "./reslice/orig");

			Image<uint16_t> out;
			reslice(orig, out, ResliceDirection::Top);
			tiff::writed(out, "./reslice/top");

			reslice(orig, out, ResliceDirection::Bottom);
			tiff::writed(out, "./reslice/bottom");

			reslice(orig, out, ResliceDirection::Left);
			tiff::writed(out, "./reslice/left");

			reslice(orig, out, ResliceDirection::Right);
			tiff::writed(out, "./reslice/right");

			Image<uint16_t> test1, test2;
			reslice(orig, test1, ResliceDirection::Top);
			reslice(test1, test2, ResliceDirection::Bottom);
			checkDifference(orig, test2, "2x resliced is not original (top-bottom).");

			reslice(orig, test1, ResliceDirection::Left);
			reslice(test1, test2, ResliceDirection::Right);
			checkDifference(orig, test2, "2x resliced is not original (left-right).");
		}

		void singleCropTest(const Vec3c& pos)
		{
			Image<uint8_t> orig(100, 100);
			ramp(orig, 0);
			tiff::writed(orig, "./crop/orig");

			Image<uint8_t> img;
			setValue(img, orig);

			Image<uint8_t> part(20, 30);
			itl2::crop(img, part, pos);
			tiff::writed(part, "./crop/part");

			draw(img, AABox(pos, pos + part.dimensions()), (uint8_t)0);
			tiff::writed(img, "./crop/orig_part_removed");

			itl2::copyValues(img, part, pos);
			tiff::writed(img, "./crop/orig_part_is_back");

			checkDifference(orig, img, string("cropped and back-copied are different, pos = ") + toString(pos));
		}

		void crop()
		{
			singleCropTest(Vec3c(50, 40, 0));
			singleCropTest(Vec3c(90, 75, 0));
			singleCropTest(Vec3c(-10, -5, 0));
			singleCropTest(Vec3c(110, 90, 0));
			singleCropTest(Vec3c(-50, 90, 0));
		}
	}
}
