
#include "itldicom.h"
using namespace std;

namespace itl2
{
	namespace dicom
	{
		namespace internals
		{
			const int PIXEL_REPRESENTATION = 0x00280103;
			const int TRANSFER_SYNTAX_UID = 0x00020010;
			const int MODALITY = 0x00080060;
			const int SLICE_THICKNESS = 0x00180050;
			const int SLICE_SPACING = 0x00180088;
			const int IMAGER_PIXEL_SPACING = 0x00181164;
			const int SAMPLES_PER_PIXEL = 0x00280002;
			const int PHOTOMETRIC_INTERPRETATION = 0x00280004;
			const int PLANAR_CONFIGURATION = 0x00280006;
			const int NUMBER_OF_FRAMES = 0x00280008;
			const int ROWS = 0x00280010;
			const int COLUMNS = 0x00280011;
			const int PIXEL_SPACING = 0x00280030;
			const int BITS_ALLOCATED = 0x00280100;
			const int WINDOW_CENTER = 0x00281050;
			const int WINDOW_WIDTH = 0x00281051;
			const int RESCALE_INTERCEPT = 0x00281052;
			const int RESCALE_SLOPE = 0x00281053;
			const int RED_PALETTE = 0x00281201;
			const int GREEN_PALETTE = 0x00281202;
			const int BLUE_PALETTE = 0x00281203;
			const int ACQUISITION_CONTEXT_SEQUENCE = 0x00400555;
			const int VIEW_CODE_SEQUENCE = 0x00540220;
			const int ICON_IMAGE_SEQUENCE = 0x00880200;
			const int ITEM = 0xFFFEE000;
			const int ITEM_DELIMINATION = 0xFFFEE00D;
			const int SEQUENCE_DELIMINATION = 0xFFFEE0DD;
			const int FLOAT_PIXEL_DATA = 0x7FE00008;
			const int PIXEL_DATA = 0x7FE00010;

			const int AE = 0x4145, AS = 0x4153, AT = 0x4154, CS = 0x4353, DA = 0x4441, DS = 0x4453, DT = 0x4454,
				FD = 0x4644, FL = 0x464C, IS = 0x4953, LO = 0x4C4F, LT = 0x4C54, PN = 0x504E, SH = 0x5348, SL = 0x534C,
				SS = 0x5353, ST = 0x5354, TM = 0x544D, UI = 0x5549, UL = 0x554C, US = 0x5553, UT = 0x5554,
				OB = 0x4F42, OW = 0x4F57, SQ = 0x5351, UN = 0x554E, QQ = 0x3F3F,
				OF = 0x4F46, OL = 0x4F4C, OD = 0x4F44, UC = 0x5543, UR = 0x5552, OV = 0x4F56, SV = 0x5356, UV = 0x5556;

			int getByte(ifstream& in)
			{
				int b = 0;
				in.read((char*)&b, 1);
				return b;
			}

			int getShort(ifstream& in, bool littleEndian)
			{
				int b0 = getByte(in);
				int b1 = getByte(in);
				if (littleEndian)
					return ((b1 << 8) + b0);
				else
					return ((b0 << 8) + b1);
			}

			int getInt(ifstream& in, bool littleEndian)
			{
				int b0 = getByte(in);
				int b1 = getByte(in);
				int b2 = getByte(in);
				int b3 = getByte(in);
				if (littleEndian)
					return ((b3 << 24) + (b2 << 16) + (b1 << 8) + b0);
				else
					return ((b0 << 24) + (b1 << 16) + (b2 << 8) + b3);
			}

			string getString(ifstream& in, int length)
			{
				vector<char> buf(length, 0);
				for (size_t n = 0; n < length; n++)
					buf[n] = getByte(in);
				return string(buf.data(), length);
			}

			vector<uint8_t> getLut(ifstream& in, int length, bool littleEndian)
			{
				if ((length & 1) != 0) { // odd
					getString(in, length);
					return vector<uint8_t>();
				}
				length /= 2;
				vector<uint8_t> lut(length);
				for (size_t i = 0; i < length; i++)
					lut[i] = (uint8_t)((unsigned int)getShort(in, littleEndian) >> 8);
				return lut;
			}

			int getLength(ifstream& in, bool littleEndian)
			{
				int b0 = getByte(in);
				int b1 = getByte(in);
				int b2 = getByte(in);
				int b3 = getByte(in);

				// We cannot know whether the VR is implicit or explicit
				// without the full DICOM Data Dictionary for public and
				// private groups.

				// We will assume the VR is explicit if the two bytes
				// match the known codes. It is possible that these two
				// bytes are part of a 32-bit length for an implicit VR.

				int vr = (b0 << 8) + b1;

				switch (vr) {
				case OB: case OW: case SQ: case UN: case UT:
				case OF: case OL: case OD: case UC: case UR:
				case OV: case SV: case UV:
					// Explicit VR with 32-bit length if other two bytes are zero
					if ((b2 == 0) || (b3 == 0))
						return getInt(in, littleEndian);
					// Implicit VR with 32-bit length
					//vr = IMPLICIT_VR;
					if (littleEndian)
						return ((b3 << 24) + (b2 << 16) + (b1 << 8) + b0);
					else
						return ((b0 << 24) + (b1 << 16) + (b2 << 8) + b3);
				case AE: case AS: case AT: case CS: case DA: case DS: case DT:  case FD:
				case FL: case IS: case LO: case LT: case PN: case SH: case SL: case SS:
				case ST: case TM: case UI: case UL: case US: case QQ:
					// Explicit vr with 16-bit length
					if (littleEndian)
						return ((b3 << 8) + b2);
					else
						return ((b2 << 8) + b3);
				default:
					// Implicit VR with 32-bit length...
					//vr = IMPLICIT_VR;
					if (littleEndian)
						return ((b3 << 24) + (b2 << 16) + (b1 << 8) + b0);
					else
						return ((b0 << 24) + (b1 << 16) + (b2 << 8) + b3);
				}
			}

			int getNextTag(ifstream& in, bool bigEndianTransferSyntax, bool& littleEndian, bool oddLocations, int& elementLength)
			{
				int groupWord = getShort(in, littleEndian);
				if (groupWord == 0x0800 && bigEndianTransferSyntax)
				{
					littleEndian = false;
					groupWord = 0x0008;
				}
				int elementWord = getShort(in, littleEndian);
				int tag = groupWord << 16 | elementWord;
				elementLength = getLength(in, littleEndian);

				// hack needed to read some GE files
				// The element length must be even!
				if (elementLength == 13 && !oddLocations)
					elementLength = 10;

				// "Undefined" element length.
				// This is a sort of bracket that encloses a sequence of elements.
				if (elementLength == -1)
				{
					elementLength = 0;
					//inSequence = true;
				}

				return tag;
			}

			double s2d(string s)
			{
				if (startsWith(s, "\\"))
					s = s.substr(1);

				return fromString<double>(s);
			}

			bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, size_t& pixelSizeBytes, uint64_t& rawDataOffset, bool& littleEndian, std::string& reason)
			{
				// NOTE: This code is ported from ImageJ (public domain)

				const int ID_OFFSET = 128;  //location of "DICM"
				const string DICM = "DICM";

				dimensions = Vec3c(0, 0, 0);
				dataType = ImageDataType::Unknown;
				pixelSizeBytes = 0;
				rawDataOffset = 0;
				reason = "";

				ifstream in(filename, ios_base::in | ios_base::binary);

				if (!in)
				{
					reason = string("Unable to open file.");
					return false;
				}

				uint8_t bytes[ID_OFFSET];
				for (int i = 0; i < ID_OFFSET; i++)
					bytes[i] = getByte(in);

				if (!in)
				{
					reason = string("No header preamble.");
					return false;
				}

				string header = getString(in, 4);

				if (!in)
				{
					reason = string("No header.");
					return false;
				}

				if(header != DICM)
				{
					if (!((bytes[0] == 8 || bytes[0] == 2) && bytes[1] == 0 && bytes[3] == 0))
					{
						reason = string("Not a DICOM or ACR/NEMA file.");
						return false;
					}

					in.seekg(0, ios_base::beg);
				}

				bool decodingTags = true;
				littleEndian = true;
				dataType = ImageDataType::UInt16;
				bool bigEndianTransferSyntax = false;
				bool oddLocations = false;
				bool isSigned = false;
				int samplesPerPixel = 0;
				while (decodingTags)
				{
					int elementLength;
					int tag = getNextTag(in, bigEndianTransferSyntax, littleEndian, oddLocations, elementLength);
					uint64_t location = in.tellg();
					if ((location & 1) != 0) // DICOM tags must be at even locations
						oddLocations = true;

					//if (inSequence && !acquisitionSequence) {
					//	addInfo(tag, null);
					//	continue;
					//}

					string s, modality, photoInterpretation, planarConfiguration;
					string spacing, center, width, intercept, slop;
					double frames = 0;
					int bitsAllocated = 16, pixelRepresentation = 0;
					switch (tag)
					{
						case TRANSFER_SYNTAX_UID:
							s = getString(in, elementLength);
							//addInfo(tag, s);
							if (contains(s, "1.2.4") || contains(s, "1.2.5"))
							{
								reason = "Cannot open compressed DICOM images. Transfer Syntax UID = " + s;
								return false;
							}
							if (contains(s, "1.2.840.10008.1.2.2"))
								bigEndianTransferSyntax = true;
							break;
						case MODALITY:
							modality = getString(in, elementLength);
							//addInfo(tag, modality);
							break;
						case NUMBER_OF_FRAMES:
							s = getString(in, elementLength);
							//addInfo(tag, s);
							frames = s2d(s);
							if (frames > 1.0)
								dimensions.z = (coord_t)frames;
							break;
						case SAMPLES_PER_PIXEL:
							samplesPerPixel = getShort(in, littleEndian);
							//addInfo(tag, samplesPerPixel);
							break;
						case PHOTOMETRIC_INTERPRETATION:
							photoInterpretation = getString(in, elementLength);
							//addInfo(tag, photoInterpretation);
							break;
						case PLANAR_CONFIGURATION:
							planarConfiguration = getShort(in, littleEndian);
							//addInfo(tag, planarConfiguration);
							break;
						case ROWS:
							dimensions.y = getShort(in, littleEndian);
							//addInfo(tag, fi.height);
							break;
						case COLUMNS:
							dimensions.x = getShort(in, littleEndian);
							//addInfo(tag, fi.width);
							break;
						case IMAGER_PIXEL_SPACING: case PIXEL_SPACING:
							s = getString(in, elementLength);
							//getSpatialScale(fi, scale);
							//addInfo(tag, scale);
							break;
						case SLICE_THICKNESS: case SLICE_SPACING:
							spacing = getString(in, elementLength);
							//fi.pixelDepth = s2d(spacing);
							//addInfo(tag, spacing);
							break;
						case BITS_ALLOCATED:
							bitsAllocated = getShort(in, littleEndian);
							if (bitsAllocated == 8)
								dataType = ImageDataType::UInt8;
							else if (bitsAllocated == 32)
								dataType = ImageDataType::UInt32;
							pixelSizeBytes = bitsAllocated / 8;
							//addInfo(tag, bitsAllocated);
							break;
						case PIXEL_REPRESENTATION:
							pixelRepresentation = getShort(in, littleEndian);
							if (pixelRepresentation == 1) {
								dataType = ImageDataType::Int16;
								isSigned = true;
							}
							//addInfo(tag, pixelRepresentation);
							break;
						case WINDOW_CENTER:
							center = getString(in, elementLength);
							//int index = center.indexOf('\\');
							//if (index != -1) center = center.substring(index + 1);
							//windowCenter = s2d(center);
							//addInfo(tag, center);
							break;
						case WINDOW_WIDTH:
							width = getString(in, elementLength);
							//index = width.indexOf('\\');
							//if (index != -1) width = width.substring(index + 1);
							//windowWidth = s2d(width);
							//addInfo(tag, width);
							break;
						case RESCALE_INTERCEPT:
							intercept = getString(in, elementLength);
							//rescaleIntercept = s2d(intercept);
							//addInfo(tag, intercept);
							break;
						case RESCALE_SLOPE:
							slop = getString(in, elementLength);
							//rescaleSlope = s2d(slop);
							//addInfo(tag, slop);
							break;
						case RED_PALETTE:
							//fi.reds =
							getLut(in, elementLength, littleEndian);
							//addInfo(tag, elementLength / 2);
							break;
						case GREEN_PALETTE:
							//fi.greens =
							getLut(in, elementLength, littleEndian);
							//addInfo(tag, elementLength / 2);
							break;
						case BLUE_PALETTE:
							//fi.blues =
							getLut(in, elementLength, littleEndian);
							//addInfo(tag, elementLength / 2);
							break;
						case FLOAT_PIXEL_DATA:
							dataType = ImageDataType::Float32;
							// continue without break
							[[fallthrough]];
						case PIXEL_DATA:
							// Start of image data...
							if (elementLength != 0)
							{
								rawDataOffset = location;
								//addInfo(tag, location);
								decodingTags = false;
							}
							//else
							//	addInfo(tag, null);
							break;
						case 0x7F880010:
							// What is this? - RAK
							if (elementLength != 0)
							{
								rawDataOffset = location + 4;
								decodingTags = false;
							}
							break;
						default:
							// Not used, skip over it...
							//addInfo(tag, null);
							break;
					}
				}

				if (dimensions.z <= 0)
					dimensions.z = 1;

				if (dataType == ImageDataType::UInt32 && isSigned)
					dataType = ImageDataType::Int32;

				if (samplesPerPixel > 1)
				{
					reason = "Multi-sample or RGB images are not supported.";
					return false;
				}

				return true;
			}
		}

		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			if (!fs::exists(filename))
				reason = "File not found.";

			// This is required in some Linux systems to differentiate files from directories.
			if (!fs::is_regular_file(filename))
			{
				reason = "Not a file.";
				return false;
			}
			
			size_t bytesPerPixel, offset;
			bool littleEndian;
			return internals::getInfo(filename, dimensions, dataType, bytesPerPixel, offset, littleEndian, reason);
		}

		namespace tests
		{
			void read()
			{
				Image<uint16_t> img;
				dicom::read(img, "../test_input_data/test.dcm");

				raw::writed(img, "../test_output_data/dicom/test_as_raw");
			}
		}
	}
}