#pragma once

#include "commandsbase.h"
#include "generation.h"
#include "distributable.h"
#include "math/vectoroperations.h"
#include "commandlist.h"
#include <numeric>

namespace pilib
{
	inline std::string genSeeAlso()
	{
		return "line, capsule, sphere, ellipsoid, set, get, ramp";
	}

	/**
	Base class for simple generation commands.
	*/
	template<typename pixel_t> class GenerationDistributable : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		GenerationDistributable(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}) :
			OneImageInPlaceCommand<pixel_t>(name, help,
				concat(extraArgs,
					{
					CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current calculation block in coordinates of the full image. This argument is used internally in distributed processing. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0)),
					}),
					genSeeAlso()
			)
		{
		}

	public:

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const
		{
			return 1;
		}
	};

	template<typename pixel_t> class SetPixelCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		SetPixelCommand() : GenerationDistributable<pixel_t>("set", "Sets a pixel in the image.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position to set."),
				CommandArgument<double>(ParameterDirection::In, "value", "Value that the pixel at given position should be set to.")
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3c pos = pop<Vec3c>(args);
			double value = pop<double>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			pos -= origin;

			if (in.isInImage(pos))
				in(pos) = pixelRound<pixel_t>(value);
		}
	};


	template<typename pixel_t> class GetPixelsCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		GetPixelsCommand() : Command("get", "Reads multiple pixels from an input image. In distributed operating mode, the output and positions images must fit into the RAM.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image where the pixels are read from."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "output", "Values of the pixels read from the input image are stored in this image. The size of the image will be set to the number of pixels read."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "positions", "Positions of pixels to read. Each row of this image contains (x, y, z) coordinates of a pixel to read from the input image. The size of the image must be 3xN where N is the count of pixels to read."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current calculation block in coordinates of the full image. This argument is used internally in distributed processing. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0))
			},
			"set, getpointsandlines")
		{
		}

		/**
		Undoes what order function did to values array.
		Orders elements such that
		unordered[index[i]] = values[i]
		Replaces values array by unordered array.
		The index array will become sorted in ascending order.
		NOTE: Works for Nx1x1 image only!
		*/
		static void unorder(Image<pixel_t>& values, vector<coord_t>& index)
		{
			for (size_t i = 0; i < index.size(); i++)
			{
				while (index[i] != i)
				{
					coord_t index_i = index[i];
					std::swap(values(i), values(index_i));
					std::swap(index[i], index[index_i]);
				}
			}
		}

		/**
		Order elements of values according to index such that
		ordered[i] = values[index[i]]
		Replaces values array by ordered array.
		The index array will become sorted in ascending order.
		This is from StackOverflow: https://stackoverflow.com/questions/46775994/reorder-array-according-to-given-index
		NOTE: Works for 3xNx1 image only!
		*/
		static void order(Image<float32_t>& values, vector<coord_t>& index)
		{
			for (size_t i = 0; i < index.size(); i++)
			{
				if (i != index[i])
				{
					size_t j = i;
					size_t k;
					while (i != (k = index[j]))
					{
						std::swap(values(0, j), values(0, k));
						std::swap(values(1, j), values(1, k));
						std::swap(values(2, j), values(2, k));
						index[j] = j;
						j = k;
					}
					index[j] = j;
				}
			}
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& out = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& positions = *pop<Image<float32_t>* >(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			
			out.ensureSize(positions.height());

			for (coord_t n = 0; n < positions.height(); n++)
			{
				Vec3d dpos(positions(0, n), positions(1, n), positions(2, n));
				Vec3c pos = round(dpos);
				pos -= origin;
				if (in.isInImage(pos))
					out(n) = in(pos);
			}
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			DistributedImage<float32_t>& positions = *std::get<DistributedImage<float32_t>* >(args[2]);

			// Simple algorithm:
			// Divide in to blocks
			// Read whole positions
			// Fill only those positions that could be read from in
			// -> out will be sparse, how to write that to disk transparently?
			// -> Does not work
			// --> hence we use the following process:
			// Sort positions according to z (, y, x), keep track of original indices
			// Divide image in to blocks
			// Corresponding out block and positions block = determined by bounds of in block
			//		The positions blocks are always continuous
			// Get pixels in each block and save as usually
			// Restore order of out such that it becomes the same than order of positions before sorting.

			// Create an array of indices
			vector<coord_t> idx(positions.height());
			std::iota(idx.begin(), idx.end(), 0);

			Image<float32_t> pos;
			positions.readTo(pos);

			// Sort the indices
			std::sort(idx.begin(), idx.end(),
				[&pos](size_t i1, size_t i2)
				{
					return itl2::vecComparer(Vec3f(pos(0, i1), pos(1, i1), pos(2, i1)), Vec3f(pos(0, i2), pos(1, i2), pos(2, i2)));
				});


			// Sort pos according to indices
			// TODO: This permutation operation and extra copy of idx array could be avoided by using a custom
			// sort or something like this: https://web.archive.org/web/20120422174751/http://www.stanford.edu/~dgleich/notebook/2006/03/sorting_two_arrays_simultaneou.html
			vector<coord_t> idx2(idx);
			order(pos, idx);
			idx.clear();
			idx.shrink_to_fit();
			positions.setData(pos);

			// Distribute. Flush before distributing so that readTo commands are ok in getCorrespondingBlock.
			distributor.flush();
			distributor.distribute(this, args);

			// Undo sorting so that the out array is in the correct order
			Image<pixel_t> outLocal;
			out.readTo(outLocal);
			unorder(outLocal, idx2);
			out.setData(outLocal);

			return vector<string>();
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const
		{
			// We have to use the input image as reference.
			return 0;
		}

		/**
		This function is given coordinates of a block in reference image (first output image in argument list or first input if there are no outputs)
		and it determines the corresponding block in another argument image.
		If this method does nothing, it is assumed that the argument image can be divided similarly than the reference image.
		@param argIndex Index of argument image.
		@param readStart, readSize File position and size of data that is loaded from disk for the reference output. Relevant only for Input and InOut images.
		@param writeFilePos, writeImPos, writeSize File position, image position and size of valid data generated by the command for the given block. Relevant only for Output and InOut images. Set writeSize to all zeroes to disable writing of output file.
		*/
		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 1 || argIndex == 2)
			{
				// Relevant block of out & positions depends on the z-range of the reference/input block.

				DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
				DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
				DistributedImage<float32_t>& positions = *std::get<DistributedImage<float32_t>* >(args[2]);

				coord_t zStart = readStart.z;
				coord_t zEnd = zStart + readSize.z;

				// TODO: This reading will be slow as this method will be called multiple times.
				Image<float32_t> pos;
				positions.readToNoFlush(pos);

				coord_t firstRow;
				for (firstRow = 0; firstRow < pos.height(); firstRow++)
				{
					coord_t z = itl2::round(pos(2, firstRow));
					if (z >= zStart)
						break;
				}

				coord_t lastRow;
				for (lastRow = firstRow + 1; lastRow < pos.height(); lastRow++)
				{
					coord_t z = itl2::round(pos(2, lastRow));
					if (z >= zEnd)
						break;
				}

				if (argIndex == 1)
				{
					// Output image dimensions are Nx1x1
					readStart = Vec3c(firstRow, 0, 0);
					readSize = Vec3c(lastRow - firstRow, 1, 1);
				}
				else
				{
					// Positions image dimensions are 3xNx1
					readStart = Vec3c(0, firstRow, 0);
					readSize = Vec3c(3, lastRow - firstRow, 1);
				}
				
				writeFilePos = readStart;
				writeSize = readSize;
				writeImPos = Vec3c(0, 0, 0);
			}
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			DistributedImage<float32_t>& positions = *std::get<DistributedImage<float32_t>* >(args[2]);

			out.ensureSize(positions.height());

			return Vec3c(0, 0, 0);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const
		{
			return JobType::Fast;
		}
	};


	template<typename pixel_t> class RampCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		RampCommand() : GenerationDistributable<pixel_t>("ramp", "Fills image with ramp in given dimension, i.e. performs img[r] = r[dimension].",
			{
				CommandArgument<coord_t>(ParameterDirection::In, "dimension", "Dimension of the ramp.", 0)
			})
		{
		}

	public:		
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			coord_t dim = pop<coord_t>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			if (dim < 0 || dim >= (coord_t)origin.size())
				throw ITLException(string("Invalid dimension: ") + itl2::toString(dim));

			ramp(in, dim, origin[dim]);
		}
	};


	template<typename pixel_t> class BoxCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		BoxCommand() : GenerationDistributable<pixel_t>("box",
			"Draws a filled axis-aligned box into the image. "
			"The filling is performed such that the left-, top- and front-faces of the box are filled but the right-, back- and bottom-faces are not. "
			"Notice that this convention is different from what is used to draw a non-axis-aligned box with the version of this command that takes box orientation vectors as arguments.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the left-top corner of the box.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "size", "Size of the box.", Vec3c(10, 10, 10)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the box.", 1)
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3c pos = pop<Vec3c>(args);
			Vec3c size = pop<Vec3c>(args);
			double value = pop<double>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			pos -= origin;

			draw(in, AABox<coord_t>(pos, pos + size), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class GenericBoxCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		GenericBoxCommand() : GenerationDistributable<pixel_t>("box",
			"Draws a filled generic non-axis-aligned box into the image. "
			"The filling is performed such that pixels on the surface of the box are not filled. "
			"Notice that this convention is different than what is used in the version of the command that does not take box orientation vectors as arguments.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "position", "Position of the left-top corner of the box.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "size", "Size of the box.", Vec3d(10, 10, 10)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the box.", 1),
				CommandArgument<Vec3d>(ParameterDirection::In, "direction 1", "Direction vector for the first axis of the box."),
				CommandArgument<Vec3d>(ParameterDirection::In, "direction 2", "Direction vector for the second axis of the box."),
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3d pos = pop<Vec3d>(args);
			Vec3d size = pop<Vec3d>(args);
			double value = pop<double>(args);
			Vec3d dir1 = pop<Vec3d>(args);
			Vec3d dir2 = pop<Vec3d>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			pos -= Vec3d(origin);

			draw(in, Box(pos, size / 2.0, dir1, dir2), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class SphereCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		SphereCommand() : GenerationDistributable<pixel_t>("sphere", "Draws a filled sphere into the image.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "position", "Position of the center point of the sphere.", Vec3d(0, 0, 0)),
				CommandArgument<double>(ParameterDirection::In, "radius", "Radius of the sphere.", 10.0),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the sphere.", 1)
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3d pos = pop<Vec3d>(args);
			double radius = pop<double>(args);
			double value = pop<double>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			pos -= Vec3d(origin);

			draw(in, Sphere(pos, radius), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class EllipsoidCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		EllipsoidCommand() : GenerationDistributable<pixel_t>("ellipsoid",
			"Draws a filled ellipsoid into the image, given position, semi-axis lengths and the directions of the first two semi-axes. "
			"The direction of the third semi-axis is the cross-product of the first two semi-axis direction vectors. "
			"The command ensures that the direction of the second semi-axis is perpendicular to the first and to the third one."
			"The lengths of the direction vectors can be anything. "
			"The filling is performed such that pixels on the surface of the ellipsoid are not filled.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "position", "Position of the center point of the ellipsoid.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "semi-axis lengths", "Lengths of the semi-axes of the ellipsoid.", Vec3d(10, 20, 30)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the ellipsoid.", 1),
				CommandArgument<Vec3d>(ParameterDirection::In, "direction 1", "Direction vector for the first semi-axis.", Vec3d(1, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "direction 2", "Direction vector for the second semi-axis.", Vec3d(0, 1, 0)),
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3d pos = pop<Vec3d>(args);
			Vec3d radius = pop<Vec3d>(args);
			double value = pop<double>(args);
			Vec3d dir1 = pop<Vec3d>(args);
			Vec3d dir2 = pop<Vec3d>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			pos -= Vec3d(origin);

			draw(in, Ellipsoid(pos, radius, dir1, dir2), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class LineCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		LineCommand() : GenerationDistributable<pixel_t>("line", "Draws a single-pixel wide line into the image.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "start", "Start position of the line.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "end", "End position of the line.", Vec3d(1, 1, 1)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels of the line.", 1)
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3d start = pop<Vec3d>(args);
			Vec3d end = pop<Vec3d>(args);
			double value = pop<double>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			start -= Vec3d(origin);
			end -= Vec3d(origin);

			draw(in, Line(start, end), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class CapsuleCommand : public GenerationDistributable<pixel_t>
	{
	protected:
		friend class CommandList;

		CapsuleCommand() : GenerationDistributable<pixel_t>("capsule", "Draws a capsule (with rounded ends) into the image.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "start", "Start position of the line.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "end", "End position of the line.", Vec3d(1, 1, 1)),
				CommandArgument<double>(ParameterDirection::In, "radius", "Radius of the capsule", 5.0),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the capsule.", 1)
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3d start = pop<Vec3d>(args);
			Vec3d end = pop<Vec3d>(args);
			double r = pop<double>(args);
			double value = pop<double>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			start -= Vec3d(origin);
			end -= Vec3d(origin);

			draw(in, Capsule(start, end, r), pixelRound<pixel_t>(value));
		}
	};






	template<typename pixel_t> class DrawGraphCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DrawGraphCommand() : OneImageInPlaceCommand<pixel_t>("drawgraph", "Draws a graph into the image. Vertices are drawn as spheres and edges as lines.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "vertices", "Image where vertex coordinates are stored. The size of the image must be 3xN, where N is the number of vertices in the graph."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "edges", "Image where vertex indices corresponding to each edge will be set. The size of the image must be 2xM where M is the number of edges. Each row of the image consists of a pair of indices to the vertex array."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "measurements", "Image that stores properties of each edge. See output from traceskeleton command."),
				CommandArgument<Image<int32_t> >(ParameterDirection::In, "edge points", "Image that stores some points on each edge. See output from traceskeleton command."),
				CommandArgument<double>(ParameterDirection::In, "vertex radius", "Radius of the spheres corresponding to the vertices.", 2),
				CommandArgument<double>(ParameterDirection::In, "vertex color", "Value that is used to fill spheres corresponding to vertices.", (double)std::numeric_limits<pixel_t>::max()),
				CommandArgument<double>(ParameterDirection::In, "edge color", "Value that is used to draw lines corresponding to edges.", (double)std::numeric_limits<pixel_t>::max()),
				CommandArgument<bool>(ParameterDirection::In, "use measured edge area", "Set to true to draw edges as capsules with cross-sectional area read from edge properties.", true),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image. This argument is used internally in distributed processing to shift the vertices to correct locations when only a part of the image is processed. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0)),
			})
		{
		}

	public:

		void run(Image<pixel_t>& in, const Image<float32_t>& vertices, const Image<uint64_t>& edges, const Image<float32_t>* meas, const Image<int32_t>* points,
			double vertexRadius, double vertexColor, double edgeColor, bool useMeasArea, Distributor::BLOCK_ORIGIN_ARG_TYPE origin) const
		{
			Network net;
			net.fromImage(vertices, edges, meas, points);

			net.vertices -= Vec3f(origin);
			draw(in, net, true, (float32_t)vertexRadius, pixelRound<pixel_t>(vertexColor), true, useMeasArea, pixelRound<pixel_t>(edgeColor));
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>*>(args);
			Image<float32_t>* meas = pop<Image<float32_t>*>(args);
			Image<int32_t>* points = pop<Image<int32_t>*>(args);
			double vertexRadius = pop<double>(args);
			double vertexColor = pop<double>(args);
			double edgeColor = pop<double>(args);
			bool useMeasArea = pop<bool>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			run(in, vertices, edges, meas, points, vertexRadius, vertexColor, edgeColor, useMeasArea, origin);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			// The input/output image is always the reference
			return 0;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			// Load images 1-4 (network definition) always
			if (argIndex == 1)
			{
				DistributedImage<float32_t>& img = *std::get<DistributedImage<float32_t>*>(args[1]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
			else if (argIndex == 2)
			{
				DistributedImage<uint64_t>& img = *std::get<DistributedImage<uint64_t>*>(args[2]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
			else if (argIndex == 3)
			{
				DistributedImage<float32_t>* img = std::get<DistributedImage<float32_t>*>(args[3]);
				if (img)
				{
					readStart = Vec3c(0, 0, 0);
					readSize = img->dimensions();
				}
			}
			else if (argIndex == 4)
			{
				DistributedImage<int32_t>* img = std::get<DistributedImage<int32_t>*>(args[4]);
				if (img)
				{
					readStart = Vec3c(0, 0, 0);
					readSize = img->dimensions();
				}
			}
		}
	};

	template<typename pixel_t> class DrawGraph2Command : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DrawGraph2Command() : OneImageInPlaceCommand<pixel_t>("drawgraph", "Draws a graph into the image. Vertices are drawn as spheres and edges as lines.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "vertices", "Image where vertex coordinates are stored. The size of the image must be 3xN, where N is the number of vertices in the graph."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "edges", "Image where vertex indices corresponding to each edge will be set. The size of the image must be 2xM where M is the number of edges. Each row of the image consists of a pair of indices to the vertex array."),
				//CommandArgument<Image<float32_t> >(ParameterDirection::In, "measurements", "Image that stores properties of each edge as output from traceskeleton command."),
				// positions...
				CommandArgument<double>(ParameterDirection::In, "vertex radius", "Radius of the spheres corresponding to the vertices.", 2),
				CommandArgument<double>(ParameterDirection::In, "vertex color", "Value that is used to fill spheres corresponding to vertices.", (double)std::numeric_limits<pixel_t>::max()),
				CommandArgument<double>(ParameterDirection::In, "edge color", "Value that is used to draw lines corresponding to edges.", (double)std::numeric_limits<pixel_t>::max()),
				//CommandArgument<bool>(ParameterDirection::In, "use measured edge area", "Set to true to draw edges as capsules with cross-sectional area read from edge properties.", true),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image. This argument is used internally in distributed processing to shift the vertices to correct locations when only a part of the image is processed. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0)),
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>*>(args);
			//Image<float32_t>* meas = pop<Image<float32_t>*>(args);
			//Image<int32_t>* positions = pop<Image<int32_t>*>(args);
			double vertexRadius = pop<double>(args);
			double vertexColor = pop<double>(args);
			double edgeColor = pop<double>(args);
			//bool useMeasArea = pop<bool>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);


			// Insert null pointer to arguments list instead of pointer to original image.
			//args.insert(args.begin() + 2, (Image<float32_t>*)0);
			//args.insert(args.begin() + 2, (Image<float32_t>*)0);
			//args.insert(args.begin() + 6, true);
			//CommandList::get<DrawGraphCommand<pixel_t> >().run(in, args);
			
			CommandList::get<DrawGraphCommand<pixel_t> >().run(in, vertices, edges, nullptr, nullptr, vertexRadius, vertexColor, edgeColor, true, origin);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			Distributable* d = (Distributable*)&CommandList::get<DrawGraphCommand<pixel_t> >();
			return d->getRefIndex(args);
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			Distributable* d = (Distributable*)&CommandList::get<DrawGraphCommand<pixel_t> >();
			return d->getCorrespondingBlock(args, argIndex, readStart, readSize, writeFilePos, writeImPos, writeSize);
		}
	};
}
