
#include "stitching.h"
#include <cmath>

namespace itl2
{
	namespace tests
	{
		//void stitchSimple()
		//{
		//	Image<uint16_t> head16;
		//	raw::read(head16, "./t1-head_bright_256x256x129.raw");

		//	Vec3d c(0, 0, 0);
		//	double t = 40.0 / 180.0 * 3.14;
		//	Matrix3x3d R(std::cos(t), -std::sin(t), 0,
		//				std::sin(t), std::cos(t), 0,
		//				0, 0, 1);
		//	double a = 1.0;
		//	Transformation<double> transformation(c, a, R);

		//	Image<uint16_t> output(head16.dimensions());
		//	Vec3c outPos(-128, 0, 0);

		//	vector<Image<uint16_t>* > images;
		//	images.push_back(&head16);

		//	vector<Transformation<double> > transformations;
		//	transformations.push_back(transformation);

		//	stitch(images, transformations, outPos, output);

		//	raw::writed(output, "./stitching/simple_out");
		//}

		//void stitchSimple2()
		//{
		//	Image<uint16_t> head16;
		//	raw::read(head16, "./t1-head_bright_256x256x129.raw");

		//	Vec3d c(0, 0, 0);
		//	double t = 40.0 / 180.0 * 3.14;
		//	Matrix3x3d R(std::cos(t), -std::sin(t), 0,
		//		std::sin(t), std::cos(t), 0,
		//		0, 0, 1);
		//	double a = 1.0;
		//	Transformation<double> transformation1(c, a, R);

		//	Transformation<double> transformation2(c, a, Matrix3x3d());

		//	Image<uint16_t> output(head16.dimensions());
		//	Vec3c outPos(-120, 100, 0);

		//	vector<Image<uint16_t>* > images;
		//	images.push_back(&head16);
		//	images.push_back(&head16);

		//	vector<Transformation<double> > transformations;
		//	transformations.push_back(transformation1);
		//	transformations.push_back(transformation2);

		//	stitch(images, transformations, outPos, output);

		//	raw::writed(output, "./stitching/simple2_out");
		//}

		void stitchFiles()
		{

			//vector<string> images;
			//images.push_back("./elastic_stitching/N0098_test1_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test2_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test3_510x510x512.raw");
			////images.push_back("./elastic_stitching/N0098_test4_510x510x512.raw");
			////images.push_back("./elastic_stitching/N0098_test5_510x510x512.raw");
			////images.push_back("./elastic_stitching/N0098_test6_510x510x512.raw");
			////images.push_back("./elastic_stitching/N0098_test7_510x510x512.raw");
			////images.push_back("./elastic_stitching/N0098_test8_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test9_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test10_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test11_510x510x512.raw");
			///*images.push_back("./elastic_stitching/N0098_test12_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test13_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test14_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test15_510x510x512.raw");
			//images.push_back("./elastic_stitching/N0098_test16_510x510x512.raw");
			//*/

			//vector<string> transf;
			//transf.push_back("./elastic_stitching/1_transformation.txt");
			//transf.push_back("./elastic_stitching/2_transformation.txt");
			//transf.push_back("./elastic_stitching/3_transformation.txt");
			////transf.push_back("./elastic_stitching/4_transformation.txt");
			////transf.push_back("./elastic_stitching/5_transformation.txt");
			////transf.push_back("./elastic_stitching/6_transformation.txt");
			////transf.push_back("./elastic_stitching/7_transformation.txt");
			////transf.push_back("./elastic_stitching/8_transformation.txt");
			//transf.push_back("./elastic_stitching/9_transformation.txt");
			//transf.push_back("./elastic_stitching/10_transformation.txt");
			//transf.push_back("./elastic_stitching/11_transformation.txt");
			///*transf.push_back("./elastic_stitching/12_transformation.txt");
			//transf.push_back("./elastic_stitching/13_transformation.txt");
			//transf.push_back("./elastic_stitching/14_transformation.txt");
			//transf.push_back("./elastic_stitching/15_transformation.txt");
			//transf.push_back("./elastic_stitching/16_transformation.txt");
			//*/


			//Image<uint16_t> output;

			//stitch<uint16_t>(images, transf, Vec3c(0, 0, 0), Vec3c(1600 / 4 * 3, 900 / 2, 1000), output, false);

			//raw::writed(output, "./elastic_stitching/test_output");

		}
	}
}
