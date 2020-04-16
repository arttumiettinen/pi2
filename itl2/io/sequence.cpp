
#include "io/sequence.h"
#include "io/alphanum.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "testutils.h"
#include <iostream>
#include <algorithm>

using namespace std;


namespace itl2
{
	namespace sequence
	{
		namespace internals
		{
			/**
			Tests if given str matches template.
			In the template
			* matches to sequence of any characters
			? matches to any single character
			@ matches to one or more numerical digits
			*/
			bool matches(const string& str, const string& templ)
			{
				size_t strPos = 0;
				size_t templPos = 0;
				while (true)
				{
					if (strPos >= str.length())
						return templPos == templ.length();

					if (templPos >= templ.length())
						return false;

					if (templ[templPos] == '?')
					{
						// Any character is ok, no need to test.
						templPos++;
						strPos++;
					}
					else if (templ[templPos] == '*')
					{
						if (templPos >= templ.length() - 1)
							return true;

						for (size_t n = strPos; n < str.length(); n++)
							if (matches(str.substr(n), templ.substr(templPos + 1)))
								return true;

						return false;
					}
					else if (templ[templPos] == '@')
					{
						if (!isdigit(str[strPos]))
							return false;

						while(strPos < str.length() && isdigit(str[strPos]))
							strPos++;

						templPos++;
					}
					else
					{
						// Characters must be equal.
						if (str[strPos] != templ[templPos])
							return false;
						strPos++;
						templPos++;
					}
				}
			}

			/**
			Separates directory and filename parts of a sequence template.
			*/
			void separatePathAndFileTemplate(const string& templ, fs::path& dir, string& fileTemplate)
			{
				fs::path p(templ);

				if (fs::is_directory(p))
				{
					dir = p;
					fileTemplate = "";
				}
				else
				{
					dir = p.parent_path();
					fileTemplate = p.filename().string();
					if (fileTemplate == ".")
						fileTemplate = "";
				}
			}

			///*
			//Separates directory and filename parts of a sequence template.
			//*/
			//void separatePathAndFileTemplate(const string& templ, fs::path& dir, string& fileTemplate)
			//{
			//	fs::path p;
			//	separatePathAndFileTemplate(templ, p, fileTemplate);
			//	dir = p.string();
			//}

			
			vector<string> buildFileList(const string& templ)
			{
				// Separate directory and file name template
				string fileTemplate;
				fs::path dir;
				separatePathAndFileTemplate(templ, dir, fileTemplate);

				if (dir == "")
					dir = ".";

				if(fileTemplate == "")
					fileTemplate = "*";
				
				//cout << "Directory: " << dir << endl;
				//cout << "Template: " << fileTemplate << endl;

				// Get those files in directory that match the template
				vector<string> filenames;

				if (fs::is_directory(dir)) // Note: This is required in Linux, or otherwise we get an exception for non-existing directories.
				{
					for (auto & p : fs::directory_iterator(dir))
					{
						string filename = p.path().filename().string();
						if (matches(filename, fileTemplate))
							filenames.push_back(p.path().string());
					}

					// Sort to natural order
					sort(filenames.begin(), filenames.end(), doj::alphanum_less<std::string>());
				}

				//for (size_t n = 0; n < filenames.size(); n++)
				//	cout << filenames[n] << endl;
				return filenames;
			}

			vector<string> buildFilteredFileList(const string& templ)
			{
				vector<string> results = buildFileList(templ);

				// Remove non-images
				//coord_t w, h;
				//ImageDataType dt;
				for (size_t n = 0; n < results.size(); n++)
				{
					// This is too slow!
					//if (!getInfo2D(result[n], w, h, dt))
					fs::path p(results[n]);
					fs::path ext = p.extension();
					string exts = ext.string();
					toLower(exts);
					if(ext != ".tif" && ext != ".tiff" && ext != ".png")
					{
						results.erase(results.begin() + n);
						n--;
					}
				}

				return results;
			}
			
			
			
			bool getInfo2D(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType, string& reason)
			{
				// TODO: Add other formats here.

				// All .png files that can be read are supported.
				if (png::getInfo(filename, width, height, dataType, reason))
					return true;

				// .tif files with single slice are supported.
				Vec3c dimensions;
				if (tiff::getInfo(filename, dimensions, dataType, reason))
				{
					width = dimensions.x;
					height = dimensions.y;
					if (dimensions.z == 1)
						return true;

					reason = "The sequence contains 3D tiff files.";
					return false;
				}

				reason = "Sequence slices could not be determined to be in any supported file format.";
				return false;
			}
		}
		
		
		bool getInfo(const string& filename, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dataType, string& reason)
	    {
		    width = 0;
		    height = 0;
		    depth = 0;
		    dataType = ImageDataType::Unknown;
		    
		    vector<string> files = internals::buildFilteredFileList(filename);

			if (files.size() <= 0)
			{
				reason = "The sequence does not match to any files.";
				return false;
			}

		    depth = files.size();
		    
		    return internals::getInfo2D(files[0], width, height, dataType, reason);
	    }
	
		

		namespace tests
		{
			void matchTest(const string& str, const string& templ, bool expected)
			{
				bool result = internals::matches(str, templ);
				cout << str << " / " << templ << " = " << result;
				if (result == expected)
					cout << endl;
				else
					cout << " FAIL" << endl;
				testAssert(result == expected, "matching");
			}

			void match()
			{
				matchTest("abc1", "abc@", true);
				matchTest("abc10", "abc@", true);
				matchTest("abc1x", "abc@", false);
				matchTest("abc1x", "abc@x", true);
				matchTest("abc11x", "abc@x", true);

				matchTest("test_sequence0000.png", "test_sequence@.png", true);
				matchTest("test_sequence_2_0014.png", "test_sequence@.png", false);
				matchTest("test_sequence2_0014.png", "test_sequence@.png", false);

				matchTest("abc", "*", true);
				matchTest("abc", "a*", true);
				matchTest("abc", "b*", false);
				matchTest("abc", "*b", false);
				matchTest("abc", "*bc", true);
				matchTest("abc", "*b*", true);
				matchTest("abc", "***", true);
				matchTest("abc", "*ab", false);

				matchTest("abc", "", false);
				matchTest("abc", "a", false);
				matchTest("abc", "ab", false);
				matchTest("abc", "abc", true);
				matchTest("abc", "a?b", false);
				matchTest("abc", "a?c", true);
				matchTest("abc", "ab?", true);
				matchTest("abc", "?bc", true);
				matchTest("abc", "abcd", false);

				matchTest("filename_001.png", "filename_*.png", true);
				matchTest("filename_002.png", "filename_*.png", true);
				matchTest("filename2_002.png", "filename_*.png", false);
				matchTest("filename_002.tif", "filename_*.png", false);
				matchTest("filename_001.png", "*.png", true);
			}

			void buildFileList()
			{
				internals::buildFileList("C:\\mytemp\\dev\\itl2\\testing\\sequence\\test_seq\\test_sequence@.png");
				internals::buildFileList("C:\\mytemp\\dev\\itl2\\testing\\sequence\\test_seq\\");
				internals::buildFileList("C:\\mytemp\\dev\\itl2\\testing\\sequence\\test_seq");
			}

			void sequence()
			{
				// NOTE: No enough asserts!

				//Image<uint8_t> img;
				//sequence::read(img, "./sequence/test_seq/test_sequence@.png");
				//raw::writed(img, "./sequence/test_sequence");

				Image<uint8_t> img2;
				sequence::read(img2, "./input_data/test_seq/t1-head_bin_@.png");
				raw::writed(img2, "./sequence/test_sequence");
				sequence::write(img2, "./sequence/test_sequence");

				// Read folder containing bad file
				writeText("./sequence/test_sequence/BAD_FILE", "");
				Image<uint8_t> withbad;
				sequence::read(withbad, "./sequence/test_sequence");
				testAssert(equals(withbad, img2), "sequence with non-image file");
				deleteFile("./sequence/test_sequence/BAD_FILE");

				sequence::write(img2, "./sequence/save_test/auto_@(-).png");
				sequence::write(img2, "./sequence/save_test/zero_@.png");
				sequence::write(img2, "./sequence/save_test/ten_@(10).png");
				sequence::write(img2, "./sequence/save_malformed/two_@(2");
				sequence::write(img2, "./sequence/save_to_folder_test/");
				sequence::write(img2, "./sequence/save_to_folder_test2");

				sequence::write(img2, "sequence_save_test");
				sequence::read(img2, "sequence_save_test");
				sequence::read(img2, "./sequence/save_to_folder_test2");

				// Partial read and write
				sequence::read(img2, "./input_data/test_seq/t1-head_bin_@.png", 100, 110);
				sequence::write(img2, "./sequence/partial_100-110/", 100);

			}

			void fileFormats()
			{
				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./input_data/t1-head_256x256x129.raw");

				sequence::write(head, "./sequence/formats/@.png");
				sequence::write(head, "./sequence/formats/@.tif");

				Image<uint16_t> out1;
				sequence::read(out1, "./sequence/formats/@.png");
				testAssert(equals(head, out1), "png sequence");

				Image<uint16_t> out2;
				sequence::read(out2, "./sequence/formats/@.tif");
				testAssert(equals(head, out2), "tif sequence");
			}

			void readWriteBlock()
			{
				// NOTE: No asserts!

				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./input_data/t1-head_256x256x129.raw");
				sequence::write(head, "./sequence/readwriteblock/write_normal/");

				Image<uint16_t> block(100, 100, 100);


				Vec3c outputDimensions = round(1.5 * Vec3d(head.dimensions()));
				string outFile = "./sequence/readwriteblock/head_3D_montage";

				Vec3c blockStart(50, 50, 0);
				Vec3c blockSize = head.dimensions() - blockStart - Vec3c(0, 10, 0);

				for (coord_t z = 0; z < 2; z++)
				{
					for (coord_t y = 0; y < 2; y++)
					{
						for (coord_t x = 0; x < 2; x++)
						{
							cout << x << "/2, " << y << "/2, " << z << "/2" << endl;
							Vec3c pos(x * blockSize.x, y * blockSize.y, z * blockSize.z);
							sequence::writeBlock(head, outFile, pos, outputDimensions, blockStart, blockSize, true);
						}
					}
				}

			}

			void readWriteBlockOptimization()
			{
				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./input_data/t1-head_256x256x129.raw");

				sequence::write(head, "./sequence/readwriteblock/write_normal/");
				sequence::writeBlock(head, "./sequence/readwriteblock/write_block/", Vec3c(0, 0, 0), head.dimensions());

				Image<uint16_t> comp1, comp2;
				sequence::read(comp1, "./sequence/readwriteblock/write_normal/");
				sequence::read(comp2, "./sequence/readwriteblock/write_block/");

				checkDifference(comp1, head, "normal writing and original");
				checkDifference(comp2, head, "optimized block writing and original");
			}
		}
	}
}
