
#include "fbp.h"
#include "pointprocess.h"
#include "fft.h"
#include "utilities.h"
#include "math/vectoroperations.h"
#include "interpolation.h"
#include "transform.h"
#include "filters.h"
#include "math/vec4.h"
#include "io/raw.h"
#include "generation.h"

#if defined(USE_OPENCL)
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#endif

#include <iostream>
#include <functional>
#include <array>
using namespace std;


namespace itl2
{
	


	ostream& operator<<(ostream& stream, const RecSettings& s)
	{
		// TODO: Do this through ImageMetaData object
		stream << "source_to_ra = " << s.sourceToRA << endl;
		stream << "rotation_direction = " << s.rotationDirection << endl;
		stream << "bhc = " << s.bhc << endl;
		stream << "rec_as_180_deg_scan = " << s.reconstructAs180degScan << endl;
		stream << "central_angle = " << s.centralAngleFor180degScan << endl;
		stream << "hswp = " << s.heuristicSinogramWindowingParameter << endl;
		stream << "rotation = " << s.rotation << endl;
		stream << "roi_center = " << s.roiCenter << endl;
		stream << "roi_size = " << s.roiSize << endl;
		stream << "crop_size = " << s.cropSize << endl;
		stream << "binning = " << s.binning << endl;
		stream << "remove_dead_pixels = " << s.removeDeadPixels << endl;
		stream << "dead_pixel_median_radius = " << s.deadPixelMedianRadius << endl;
		stream << "dead_pixel_std_dev_count = " << s.deadPixelStdDevCount << endl;
		stream << endl;
		stream << "center_shift = " << s.centerShift << endl;
		stream << "camera_z_shift = " << s.cameraZShift << endl;
		stream << "camera_rotation = " << s.cameraRotation << endl;
		stream << "cs_angle_slope = " << s.csAngleSlope << endl;
		stream << endl;
		stream << "pad_type = " << s.padType << endl;
		stream << "pad_size = " << s.padFraction << endl;
		stream << "filter_type = " << s.filterType << endl;
		stream << "filter_cut_off = " << s.filterCutOff << endl;
		stream << endl;
		stream << "phase_mode = " << s.phaseMode << endl;
		stream << "phase_pad_type = " << s.phasePadType << endl;
		stream << "phase_pad_size = " << s.phasePadFraction << endl;
		stream << "propagation_distance = " << s.objectCameraDistance << endl;
		stream << "delta = " << s.delta << endl;
		stream << "mu = " << s.mu << endl;
		stream << endl;
		stream << "range_min = " << s.dynMin << endl;
		stream << "range_max = " << s.dynMax << endl;
		stream << endl;
		stream << "shift_scale = " << s.shiftScaling << endl;
		stream << "use_shifts = " << s.useShifts << endl;
		stream << endl;
		stream << "angles" << endl;

		stream.precision(17);
		for (size_t n = 0; n < s.angles.size(); n++)
			stream << fixed << s.angles[n] << endl;
		stream << endl;

		stream << "sample_shifts" << endl;
		for (size_t n = 0; n < s.sampleShifts.size(); n++)
			stream << fixed << s.sampleShifts[n] << endl;
		stream << endl;

		stream << "source_shifts" << endl;
		for (size_t n = 0; n < s.sourceShifts.size(); n++)
			stream << fixed << s.sourceShifts[n] << endl;
		stream << endl;

		stream << "camera_shifts" << endl;
		for (size_t n = 0; n < s.cameraShifts.size(); n++)
			stream << fixed << s.cameraShifts[n] << endl;
		stream << endl;

		stream << "rotation_axis_shifts" << endl;
		for (size_t n = 0; n < s.rotationAxisShifts.size(); n++)
			stream << fixed << s.rotationAxisShifts[n] << endl;
		stream << endl;

		return stream;
	}

	/**
	Performs phase retrieval for single projection slice.
	*/
	void paganinSlice(Image<float32_t>& slice, PadType padType, float32_t padFraction, float32_t objectSourceDistance, float32_t objectCameraDistance, float32_t delta, float32_t mu)
	{
		clamp(padFraction, 0.0f, 1.0f);
		size_t mirrorSizeX = (size_t)round(padFraction * slice.width());
		size_t mirrorSizeY = (size_t)round(padFraction * slice.height());

		size_t paddedWidth = mirrorSizeX + slice.width() + mirrorSizeX;
		size_t paddedHeight = mirrorSizeY + slice.height() + mirrorSizeY;

		size_t width = (size_t)slice.width();
		size_t height = (size_t)slice.height();
		
		Image<float32_t> in(paddedWidth, paddedHeight);
		Image<complex32_t> out;

		for (size_t z0 = 0; z0 < paddedHeight; z0++)
		{
			int z = (int)z0 - (int)mirrorSizeY;
			if (padType == PadType::Mirror)
			{
				if (z < 0)
					z = -z;
				else if (z >= (int)height)
					z = (int)height - 1 - (z - (int)height);
			}
			else if (padType == PadType::Nearest)
			{
				if (z < 0)
					z = 0;
				else if (z >= (int)height)
					z = (int)height - 1;
			}
			else
				throw runtime_error("Pad type not supported.");

			// Copy data to input buffer and do padding
			for (size_t x = 0; x < mirrorSizeX; x++)
			{
				if (padType == PadType::Mirror)
					in(mirrorSizeX - x - 1, z0) = slice(x, z);
				else if (padType == PadType::Nearest)
					in(mirrorSizeX - x - 1, z0) = slice(0, z); //in[z0 * paddedWidth + mirrorSize - x - 1]
				else
					throw runtime_error("Pad type not supported.");
			}
			for (size_t x = 0; x < width; x++)
			{
				in(mirrorSizeX + x, z0) = slice(x, z);
			}
			for (size_t x = 0; x < mirrorSizeX; x++)
			{
				if (padType == PadType::Mirror)
					in(mirrorSizeX + width + x, z0) = slice(width - x - 1, z);
				else if (padType == PadType::Nearest)
					in(mirrorSizeX + width + x, z0) = slice(width - 1, z);
				else
					throw runtime_error("Pad type not supported.");
			}
		}

		float32_t M = (objectCameraDistance + objectSourceDistance) / objectSourceDistance;

		multiply(in, M*M);

		// FFT
		fft(in, out);

		// Multiply by (1 + 4pi^2 d delta/mu |w|^2)^-1
		// w = 0 are at (0, 0), and (0, outHeight)
		for (coord_t y = 0; y < out.height(); y++)
		{
			for (coord_t x = 0; x < out.width(); x++)
			{
				double dx = (double)(x - 0);
				double dy = (double)std::min(y - 0, out.height() - y);
				double w2 = dx * dx + dy * dy;

				double trans = 1 / (1 + 4 * PI * PI * objectCameraDistance * delta / mu * w2 / M);

				out(x, y) *= (float)trans;
			}
		}

		// Inverse FFT
		ifft(out, in);

		// Remove padding and replace original data
		for (size_t z = 0; z < height; z++)
		{
			for (size_t x = 0; x < width; x++)
			{
				slice(x, z) = in(x + mirrorSizeX, z + mirrorSizeY);
			}
		}

		// Take negative logarithm
		negLog(slice);
	}


	/**
	Performs Paganin phase retrieval.
	*/
	void paganin(Image<float32_t>& transmissionProjections, PadType padType, float32_t padFraction, float32_t objectSourceDistance, float32_t objectCameraDistance, float32_t delta, float32_t mu)
	{
		cout << "Single material Paganin phase-retrieval..." << endl;

		// Calculate Paganin (single material) filtering

		size_t counter = 0;
		#pragma omp parallel for
		for (coord_t anglei = 0; anglei < transmissionProjections.depth(); anglei++)
		{
			Image<float32_t> slice(transmissionProjections, anglei, anglei);

			paganinSlice(slice, padType, padFraction, objectSourceDistance, objectCameraDistance, delta, mu);

			showThreadProgress(counter, transmissionProjections.depth());
		}
	}

	/**
	Performs phase retrieval or -log operation on a single transmission projection slice.
	*/
	void phaseRetrievalSlice(Image<float32_t>& slice, PhaseMode phaseMode, PadType padType, float32_t padFraction, float32_t objectSourceDistance, float32_t objectCameraDistance, float32_t delta, float32_t mu)
	{
		switch (phaseMode)
		{
		case PhaseMode::Absorption:
			negLog(slice);
			break;
		case PhaseMode::Paganin:
			paganinSlice(slice, padType, padFraction, objectSourceDistance, objectCameraDistance, delta, mu);
			break;
		case PhaseMode::Direct:
			// The data is already in -ln(I/I0) format or equivalent, so do not do anything.
			break;
		default:
			throw ITLException("Unsupported phase retrieval mode.");
		}
	}

	/**
	Performs phase retrieval or -log operation on transmission projections.
	*/
	void phaseRetrieval(Image<float32_t>& transmissionProjections, PhaseMode phaseMode, PadType padType, float32_t padFraction, float32_t objectSourceDistance, float32_t objectCameraDistance, float32_t delta, float32_t mu)
	{
		switch (phaseMode)
		{
		case PhaseMode::Absorption:
			negLog(transmissionProjections);
			break;
		case PhaseMode::Paganin:
			paganin(transmissionProjections, padType, padFraction, objectSourceDistance, objectCameraDistance, delta, mu);
			break;
		case PhaseMode::Direct:
			// The data is already in -ln(I/I0) format or equivalent, so do not do anything.
			break;
		default:
			throw ITLException("Unsupported phase retrieval mode.");
		}
	}

	/*
	Performs beam hardening correction for data for which -log has been taken. Replaces original data.
	*/
	void beamHardeningCorrection(Image<float32_t>& transmissionProjections, float32_t bhc)
	{
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < transmissionProjections.pixelCount(); n++)
		{
			float32_t x = transmissionProjections(n);

			x = x + bhc * x * x;
			
			if (!isfinite(x))
				x = 0;

			transmissionProjections(n) = x;
		}
	}

	namespace internals
	{

		/**
		Calculate true central angle for scan from user-specified central angle.
		@param gammamax0 Maximum half cone angle on optical axis.
		*/
		float32_t calculateTrueCentralAngle(float32_t centralAngleFor180degScan, const vector<float32_t>& angles, float32_t gammamax0)
		{
			float32_t minAngle = min(angles) + 90 + gammamax0;
			float32_t maxAngle = max(angles) - 90 - gammamax0;

			clamp(centralAngleFor180degScan, minAngle, maxAngle);

			return centralAngleFor180degScan;
		}

		/**
		Calculates maximum half cone angle on optical axis.
		@param d Distance between object and cone tip.
		*/
		float32_t calculateGammaMax0(float32_t projectionWidth, float32_t d)
		{
			return atan2f(projectionWidth / 2.0f, d) / PIf * 180.0f;
		}


		float32_t csAnglePerturbation(size_t angleIndex, float32_t centralAngle, const vector<float32_t>& angles, float32_t csAngleSlope)
		{
			return (angles[angleIndex] - centralAngle) * csAngleSlope;
		}

		/**
		Checks that projection images and settings correspond to each other.
		Adjusts zero elements in roi size vector to full image dimension.
		*/
		void sanityCheck(const Image<float32_t>& transmissionProjections, RecSettings& settings, bool projectionsAreBinned)
		{
			if (transmissionProjections.depth() != settings.angles.size())
				throw ITLException("Count of projection images and count of angles do not match.");

			size_t projCount = settings.angles.size();

			if (settings.sampleShifts.size() > 0 && settings.sampleShifts.size() != projCount)
				throw ITLException("Count of sample shifts and count of angles do not match.");

			if (settings.rotationAxisShifts.size() > 0 && settings.rotationAxisShifts.size() != projCount)
				throw ITLException("Count of rotation axis shifts and count of angles do not match.");

			if (settings.cameraShifts.size() > 0 && settings.cameraShifts.size() != projCount)
				throw ITLException("Count of camera shifts and count of angles do not match.");

			if (settings.sourceShifts.size() > 0 && settings.sourceShifts.size() != projCount)
				throw ITLException("Count of source shifts and count of angles do not match.");

			if (NumberUtils<float32_t>::lessThanOrEqual(settings.sourceToRA, 0))
				throw ITLException("Non-positive source to rotation axis distance.");

			if (settings.binning < 1)
				throw ITLException("Binning must be at least 1.");

			// Roi size and position
			if (settings.roiSize.x <= 0)
				settings.roiSize.x = projectionsAreBinned ? transmissionProjections.width() * settings.binning : transmissionProjections.width();
			if (settings.roiSize.x <= 0)
				settings.roiSize.x = 1;
			if (settings.roiSize.y <= 0)
				settings.roiSize.y = projectionsAreBinned ? transmissionProjections.width() * settings.binning : transmissionProjections.width();
			if (settings.roiSize.y <= 0)
				settings.roiSize.y = 1;
			if (settings.roiSize.z <= 0)
				settings.roiSize.z = projectionsAreBinned ? transmissionProjections.height() * settings.binning : transmissionProjections.height();
			if (settings.roiSize.z <= 0)
				settings.roiSize.z = 1;

		}
	}

	/**
	Applies cone-beam related weighting of -log data (after beam hardening correction, before filtering) for single projection image slice.
	*/
	void fbpWeightingSlice(Image<float32_t>& slice, coord_t angleIndex, bool reconstructAs180DegScan, const vector<float32_t>& angles, float32_t baseCenterShift, float32_t csAngleSlope, float32_t sourceToRA, float32_t cameraZShift, float32_t centralAngle, float32_t gammaMax0, float32_t heuristicSinogramWindowingParameter)
	{
		// TODO: How to calculate these in the new geometry? These shifts are probably not very important...
		//float32_t dx = settings.objectShifts[anglei].x;
		//float32_t dz = settings.objectShifts[anglei].y;
		float32_t dx = 0;
		float32_t dz = 0;

		float proj_weight = 1;
		if (!reconstructAs180DegScan)
		{
			// Remove data outside of 360 deg range
			if (angleIndex > 0)
			{
				if (NumberUtils<float32_t>::greaterThan(abs(angles[angleIndex] - angles[0]), 360 - (angles[1] - angles[0]), 0.1f * abs(angles[1] - angles[0])))
				{
					cout << "Excess projection removed (index = " << angleIndex << ", angle = " << angles[angleIndex] << " deg)." << endl;
					proj_weight = 0;
				}
			}
		}

		float centerShift = baseCenterShift + internals::csAnglePerturbation(angleIndex, centralAngle, angles, csAngleSlope);

		//cout << "angle = " << angles[phi] << ", cs = " << centerShift << endl;

		float32_t d = sourceToRA;
		coord_t height = slice.height();
		coord_t width = slice.width();

		for (coord_t z = 0; z < height; z++)
		{
			float Z = (float)z - (height / 2.0f - dz) - cameraZShift;
			for (coord_t y = 0; y < width; y++)
			{
				float Y = (float)y - ((float)width / 2.0f - dx) - centerShift;

				// Weight inherent in FDK algorithm
				float w = d / sqrt(d * d + Y * Y + Z * Z);


				float sin_weight = 1;
				if (reconstructAs180DegScan)
				{
					// Smooth sinogram windowing for 180 deg scan


					// Original version with symmetric smoothing
					float f = abs(Z / (0.5f * height));
					float gammamax = (1 - f) * gammaMax0 + f * heuristicSinogramWindowingParameter * gammaMax0;

					// Angle from view direction to current ray
					float gamma = atan2f(Y, d) / PIf * 180.0f;

					// Angle of the current projection
					float beta = (float)(angles[angleIndex] - (centralAngle - 90));

					// Determine weight
					if (-gammamax <= beta && beta <= gammamax - 2 * gamma) // Left side shape
					{
						float sinterm = sinf((45 * ((beta + gammamax) / (gammamax - gamma))) / 180 * PIf);
						sin_weight = sinterm * sinterm;
					}
					else if (-gammamax <= beta - 360 && beta - 360 <= gammamax - 2 * gamma)	// Left side shape modulo
					{
						float sinterm = sinf((45 * ((beta - 360) + gammamax) / (gammamax - gamma)) / 180 * PIf);
						sin_weight = sinterm * sinterm;
					}
					else if (gammamax - 2 * gamma <= beta && beta <= 180 - gammamax - 2 * gamma ||					// Middle flat
						gammamax - 2 * gamma <= beta + 360 && beta + 360 <= 180 - gammamax - 2 * gamma ||		// Middle flat modulo
						gammamax - 2 * gamma <= beta - 360 && beta - 360 <= 180 - gammamax - 2 * gamma)			// Middle flat modulo
					{
						sin_weight = 1;
					}
					else if (180 - gammamax - 2 * gamma <= beta && beta <= 180 + gammamax) // Right side shape
					{
						float sinterm = sinf((45 * ((180 + gammamax - beta) / (gammamax + gamma))) / 180 * PIf);
						sin_weight = sinterm * sinterm;
					}
					else if (180 - gammamax - 2 * gamma <= beta + 360 && beta + 360 <= 180 + gammamax)	// Right side shape modulo
					{
						float sinterm = sinf((45 * (180 + gammamax - (beta + 360)) / (gammamax + gamma)) / 180 * PIf);
						sin_weight = sinterm * sinterm;
					}
					else
					{
						sin_weight = 0;
					}
				}

				slice(y, z) *= w * sin_weight * proj_weight;
			}
		}
	}

	/**
	Applies cone-beam related weighting of -log data (after beam hardening correction, before filtering).
	*/
	void fbpWeighting(Image<float32_t>& transmissionProjections, bool reconstructAs180DegScan, const vector<float32_t>& angles, float32_t centerShift, float32_t csAngleSlope, float32_t sourceToRA, float32_t cameraZShift, float32_t centralAngle, float32_t gammaMax0, float32_t heuristicSinogramWindowingParameter)
	{
		// Apply weighting
		#pragma omp parallel for
		for (int anglei = 0; anglei < transmissionProjections.depth(); anglei++)
		{
			Image<float32_t> slice(transmissionProjections, anglei, anglei);
			fbpWeightingSlice(slice, anglei, reconstructAs180DegScan, angles, centerShift, csAngleSlope, sourceToRA, cameraZShift, centralAngle, gammaMax0, heuristicSinogramWindowingParameter);
		}
	}

	/**
	Constructs a backprojection filter into given image. Image width must be set to correct value.
	@param cutoff Frequency cutoff, 1 corresponds to Nyquist frequency, 0 to DC.
	*/
	void createFilter(Image<float32_t>& filter, FilterType filterType, float32_t cutoff)
	{

		switch (filterType)
		{
			case FilterType::IdealRamp:
			{
				// Normal sharp |w| filter
				// This filter was used in the old versions.
				double k1 = (1.0 - 0.0) / (filter.width() - 1 - 0.0);
				for (coord_t x = 0; x < filter.width(); x++)
				{
					float r = (float)(k1*x);
					filter(x) = r;
				}

				break;
			}
			case FilterType::Ramp:
			{
				initFFTW();
				coord_t s = filter.width() - 1;
				coord_t pow2 = s << 1;
				Image<complex32_t> F(pow2);

				fftwf_plan plan;
				#pragma omp critical
				{
					plan = fftwf_plan_dft_1d((int)pow2, (fftwf_complex*)F.getData(), (fftwf_complex*)F.getData(), FFTW_FORWARD, FFTW_ESTIMATE);
				}

				F(0) = 0.25;
				for (coord_t i = 1; i < F.width(); i++)
				{
					if (i <= s)
					{
						if (i & 0x1)
							F(i) = -1.0f / (PIf * PIf * (float32_t)i * (float32_t)i);
						else
							F(i) = 0;
					}
					else
						F(i) = F(pow2 - i);
				}

				fftwf_execute(plan);

				for (coord_t n = 0; n < filter.width(); n++)
					filter(n) = 2 * F(n).real();

				#pragma omp critical
				{
					fftwf_destroy_plan(plan);
				}

				break;
			}
			case FilterType::SheppLogan:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double factor = (PI * filter(i)) / (2.0 * cutoff);
					filter(i) *= (float32_t)(sin(factor) / factor);
				}

				break;
			}
			case FilterType::Cosine:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double factor = (PI * filter(i)) / (2.0 * cutoff);
					filter(i) *= (float32_t)cos(factor);
				}

				break;
			}
			case FilterType::Hamming:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double factor = (PI * filter(i)) / cutoff;
					filter(i) *= (float32_t)(0.54 + 0.46 * cos(factor));
				}

				break;
			}
			case FilterType::Hann:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double factor = (PI * filter(i)) / cutoff;
					filter(i) *= (float32_t)(0.5 + 0.5 * cos(factor));
				}

				break;
			}
			case FilterType::Blackman:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double factor = (PI * filter(i)) / cutoff + PI;
					filter(i) *= (float32_t)(0.42 - 0.5 * cos(factor) + 0.08 * cos(2 * factor));
				}

				break;
			}
			case FilterType::Parzen:
			{
				createFilter(filter, FilterType::Ramp, cutoff);

				for (coord_t i = 1; i < filter.width(); i++)
				{
					double L = 2.0 * (double)filter.width() * cutoff;
					double n = (double)i;// filter(i);
					double q = n / (L / 2);

					double w;
					if (n <= L / 4)
						w = 1 - 6 * q * q * (1 - q);
					else
						w = 2 * (1 - q) * (1 - q) * (1 - q);
					filter(i) *= (float32_t)w;
				}

				break;
			}
			default:
			{
				throw ITLException("Unsupported backprojection filter.");
			}
		}
		
		// Reset everything above cutoff
		coord_t ind = pixelRound<coord_t>(cutoff * (filter.width() - 1));
		
		for (coord_t n = ind + 1; n < filter.width(); n++)
			filter(n) = 0;
	}

	struct FilterSettings
	{
		/**
		Input buffer.
		*/
		Image<float32_t> in;

		/**
		FFT buffer.
		*/
		Image<complex32_t> out;

		/**
		Filter that is to be applied.
		Currently this is created separately for each thread but the data is the same for all of them.
		*/
		Image<float32_t> H;

		/**
		FFT plans.
		*/
		fftwf_plan forward;
		fftwf_plan backward;

		/**
		Amount of padding on each size of the buffers in pixels.
		*/
		coord_t mirrorSize;

		FilterSettings(float32_t padFraction, coord_t projectionWidth, FilterType filterType, float32_t cutoff)
		{
			clamp(padFraction, 0.0f, 1.0f);
			mirrorSize = (coord_t)round(padFraction * projectionWidth);
			if (mirrorSize < 10)
				mirrorSize = 10;

			size_t paddedWidth = mirrorSize + projectionWidth + mirrorSize;

			in.ensureSize(paddedWidth);

			// Build filter
			size_t outSize = paddedWidth / 2 + 1;

			out.ensureSize(outSize);
			H.ensureSize(outSize);

			createFilter(H, filterType, cutoff);

			#pragma omp critical
			{
				forward = fftwf_plan_dft_r2c_1d((int)paddedWidth, in.getData(), (fftwf_complex*)out.getData(), FFTW_ESTIMATE);
				backward = fftwf_plan_dft_c2r_1d((int)paddedWidth, (fftwf_complex*)out.getData(), in.getData(), FFTW_ESTIMATE);
			}
		}

		~FilterSettings()
		{
			#pragma omp critical
			{
				fftwf_destroy_plan(backward);
				fftwf_destroy_plan(forward);
			}
		}
	};


	void filterSlice(Image<float32_t>& slice, FilterSettings& settings, PadType padType)
	{
		coord_t projectionWidth = slice.width();

		for (coord_t z = 0; z < slice.height(); z++)
		{
			// Copy data to input buffer
			for (coord_t x = 0; x < settings.mirrorSize; x++)
			{
				if (padType == PadType::Mirror)
					settings.in(settings.mirrorSize - x - 1) = slice(x, z);
				else if (padType == PadType::Nearest)
					settings.in(settings.mirrorSize - x - 1) = slice(0, z);
				else
					throw runtime_error("Pad type not supported.");
			}
			for (coord_t x = 0; x < projectionWidth; x++)
			{
				settings.in(settings.mirrorSize + x) = slice(x, z);
			}
			for (coord_t x = 0; x < settings.mirrorSize; x++)
			{
				if (padType == PadType::Mirror)
					settings.in(settings.mirrorSize + projectionWidth + x) = slice(projectionWidth - x - 1, z);
				else if (padType == PadType::Nearest)
					settings.in(settings.mirrorSize + projectionWidth + x) = slice(projectionWidth - 1, z);
				else
					throw runtime_error("Pad type not supported.");
			}

			// Transform
			// NOTE: Output stores only non-negative frequencies!
			fftwf_execute(settings.forward);

			// Apply filter
			for (coord_t x = 0; x < settings.out.width(); x++)
				settings.out(x) *= settings.H(x);

			// Inverse transform
			fftwf_execute(settings.backward);

			// Normalize and copy to output
			for (coord_t x = 0; x < projectionWidth; x++)
			{
				slice(x, z) = settings.in(settings.mirrorSize + x) / settings.in.width();
			}
		}
	}

	/**
	Performs sinogram filtering.
	*/
	void fbpFilter(Image<float32_t>& transmissionProjections, PadType padType, float32_t padFraction, FilterType filterType, float32_t cutoff)
	{
		initFFTW();

		#pragma omp parallel
		{
			FilterSettings settings(padFraction, transmissionProjections.width(), filterType, cutoff);

			#pragma omp for
			for (coord_t anglei = 0; anglei < transmissionProjections.depth(); anglei++)
			{
				Image<float32_t> slice(transmissionProjections, anglei, anglei);
				filterSlice(slice, settings, padType);
			}			
		}
	}

	/**
	Remove bad pixels from one slice of projection data.
	@param slice View of the slice.
	@param med, tmp Temporary images.
	@return Count of bad pixels in the slice.
	*/
	size_t deadPixelRemovalSlice(Image<float32_t>& slice, Image<float32_t>& med, Image<float32_t>& tmp, coord_t medianRadius = 2, float32_t stdDevCount = 30)
	{
		// Calculate median filtering of slice
		nanMedianFilter(slice, med, medianRadius, NeighbourhoodType::Rectangular, BoundaryCondition::Nearest);

		// Calculate abs(slice - median)
		setValue(tmp, slice);
		subtract(tmp, med);
		abs(tmp);

		// Calculate its mean and standard deviation
		Vec2d v = maskedMeanAndStdDev(tmp, numeric_limits<float32_t>::signaling_NaN());
		//float32_t meandifference = (float32_t)v.x;
		float32_t stddifference = (float32_t)v.y;


		// Perform filtering
		size_t badPixelCount = 0;
		for (coord_t y = 0; y < slice.height(); y++)
		{
			for (coord_t x = 0; x < slice.width(); x++)
			{
				float32_t p = slice(x, y);
				float32_t m = med(x, y);

				if (NumberUtils<float32_t>::isnan(p) || abs(m - p) > stdDevCount * stddifference)
				{
					slice(x, y) = m;
					badPixelCount++;
				}
			}
		}

		return badPixelCount;
	}

	void printBadPixelInfo(float averageBadPixels, size_t maxBadPixels)
	{
		cout << "Average number of bad pixels per slice: " << averageBadPixels << endl;
		cout << "Maximum number of bad pixels in slice: " << maxBadPixels << endl;

		if (maxBadPixels > 100)
			cout << "WARNING: Maximum number of bad pixels is high: " << maxBadPixels << ". Consider changing settings for bad pixel removal." << endl;
	}

	/**
	Removes dead pixels from each slice of the input image.
	*/
	void deadPixelRemoval(Image<float32_t>& img, coord_t medianRadius = 1, float32_t stdDevCount = 30)
	{
		if (medianRadius <= 0)
			throw ITLException("Invalid median radius.");
		if (stdDevCount <= 0)
			throw ITLException("Invalid standard deviation count.");

		float32_t averageBadPixels = 0;
		size_t maxBadPixels = 0;

		size_t counter = 0;
		#pragma omp parallel
		{

			Image<float32_t> med(img.width(), img.height());
			Image<float32_t> tmp(img.width(), img.height());
			
			#pragma omp for
			for (coord_t n = 0; n < img.depth(); n++)
			{
				Image<float32_t> slice(img, n, n);

				size_t badPixelCount = deadPixelRemovalSlice(img, med, tmp, medianRadius, stdDevCount);

				#pragma omp critical(badpixels)
				{
					averageBadPixels += badPixelCount;
					maxBadPixels = std::max(maxBadPixels, badPixelCount);
				}

				showThreadProgress(counter, img.depth());
			}
		}

		printBadPixelInfo(averageBadPixels / (float)img.depth(), maxBadPixels);
	}

	

	namespace internals
	{
		void applyBinningToParameters(RecSettings& settings)
		{
			if (settings.binning > 1)
			{
				// Apply binning to all parameters
				float32_t b = (float32_t)settings.binning;
				settings.cameraZShift /= b;
				settings.centerShift /= b;
				//settings.cropSize /= b;
				settings.csAngleSlope /= b;
				settings.objectCameraDistance /= b;
				for (size_t n = 0; n < settings.sampleShifts.size(); n++)
					settings.sampleShifts[n] /= b;
				settings.roiCenter /= (coord_t)settings.binning;
				settings.roiSize /= (coord_t)settings.binning;
				
				if (settings.roiSize.x < 1)
					settings.roiSize.x = 1;

				if (settings.roiSize.y < 1)
					settings.roiSize.y = 1;

				if (settings.roiSize.z < 1)
					settings.roiSize.z = 1;

				settings.sourceToRA /= b;
				settings.binning = 1;

				settings.delta *= b;
				settings.mu *= b;

				settings.deadPixelMedianRadius = std::max((coord_t)1, (coord_t)round(settings.deadPixelMedianRadius / b));
			}
		}
	}

	/**
	Caps values to range ]0, 1[.
	*/
	void replaceBadValues(Image<float32_t>& img)
	{
		float32_t eps = 0.00001f;
		for(coord_t n = 0; n < img.pixelCount(); n++)
		{
			float32_t& p = img(n);
			if(p < eps)
				p = eps;
			else if(p > 1 - eps)
				p = 1 - eps;
		}
	}
	
	void fbpPreprocess(const Image<float32_t>& transmissionProjections, Image<float32_t>& preprocessedProjections, RecSettings settings)
	{
		internals::sanityCheck(transmissionProjections, settings, false);

		size_t origBinning = settings.binning;

		// Calculate size of preprocessed projections
		Vec3c croppedSize(transmissionProjections.dimensions());
		if (settings.cropSize.max() > 0)
		{
			croppedSize.x = transmissionProjections.width() - 2 * settings.cropSize.x;
			croppedSize.y = transmissionProjections.height() - 2 * settings.cropSize.y;

			//if (croppedSize.x / settings.binning <= 0 || croppedSize.y / settings.binning <= 0)
			//	throw ITLException("Too large crop size. Cropped projection image size must be at least 1x1 pixels after cropping and binning.");
		}

		Vec3c outputSize = croppedSize;
		if (settings.binning > 1)
		{
			outputSize.x /= (coord_t)settings.binning;
			outputSize.y /= (coord_t)settings.binning;
		}

		if (outputSize.min() <= 0)
			throw ITLException("Too large crop size or binning. A projection image must have at least 1x1 pixels after cropping and binning.");
		
		// Adjust parameters for binning
		internals::applyBinningToParameters(settings);

		preprocessedProjections.ensureSize(outputSize);

		// Maximum angle for any ray
		float32_t gammamax0 = internals::calculateGammaMax0((float32_t)preprocessedProjections.width(), settings.sourceToRA);

		cout << "Maximum half cone angle on optical axis = " << gammamax0 << " deg" << endl;

		float32_t centralAngle = internals::calculateTrueCentralAngle(settings.centralAngleFor180degScan, settings.angles, gammamax0);

		if (settings.reconstructAs180degScan)
			cout << "Central angle for 180 deg reconstruction: " << centralAngle << " deg (available angular range = " << min(settings.angles) << " deg - " << max(settings.angles) << " deg)" << endl;

		cout << "Preprocessing..." << endl;

		// Process slice by slice to reduce disk I/O when the transmission projection image is memory-mapped.
		float32_t averageBadPixels = 0;
		size_t maxBadPixels = 0;
		size_t counter = 0;
		#pragma omp parallel
		{
			Image<float32_t> med;
			Image<float32_t> tmp;
			FilterSettings filterSettings(settings.padFraction, preprocessedProjections.width(), settings.filterType, settings.filterCutOff);

			Image<float32_t> cropTmp(croppedSize.x, croppedSize.y);

			#pragma omp for
			for (coord_t z = 0; z < preprocessedProjections.depth(); z++)
			{
				//cout << "Processing slice " << z << endl;

				Image<float32_t> origSlice(transmissionProjections, z, z);
				Image<float32_t> slice(preprocessedProjections, z, z);

				if (settings.cropSize.max() > 0 && origBinning <= 1)
				{
					// Cropping but no binning
					crop(origSlice, slice, Vec3c(settings.cropSize.x, settings.cropSize.y, 0));
				}
				else if (settings.cropSize.max() <= 0 && origBinning > 1)
				{
					// Binning but no cropping
					binning(origSlice, slice, Vec3c(origBinning, origBinning, 1), false);
				}
				else if (settings.cropSize.max() > 0 && origBinning > 1)
				{
					// Cropping and binning
					crop(origSlice, cropTmp, Vec3c(settings.cropSize.x, settings.cropSize.y, 0));
					binning(cropTmp, slice, Vec3c(origBinning, origBinning, 1), false);
				}
				else
				{
					// No binning, no cropping
					setValue(slice, origSlice);
				}
				
				replaceBadValues(slice);

				if(settings.removeDeadPixels)
				{
					size_t badPixelCount = deadPixelRemovalSlice(slice, med, tmp, settings.deadPixelMedianRadius, settings.deadPixelStdDevCount);
#pragma omp critical(badpixelsslice)
					{
						averageBadPixels += badPixelCount;
						maxBadPixels = std::max(maxBadPixels, badPixelCount);
					}
				}

				phaseRetrievalSlice(slice, settings.phaseMode, settings.phasePadType, settings.phasePadFraction, settings.sourceToRA, settings.objectCameraDistance, settings.delta, settings.mu);
				
				if(!NumberUtils<float32_t>::equals(settings.bhc, 0))
					beamHardeningCorrection(slice, settings.bhc);

				fbpWeightingSlice(slice, z, settings.reconstructAs180degScan, settings.angles, settings.centerShift, settings.csAngleSlope, settings.sourceToRA, settings.cameraZShift, centralAngle, gammamax0, settings.heuristicSinogramWindowingParameter);
				
				filterSlice(slice, filterSettings, settings.padType);
				
				showThreadProgress(counter, preprocessedProjections.depth());
			}
		}

		if (settings.removeDeadPixels)
			printBadPixelInfo(averageBadPixels / (float)preprocessedProjections.depth(), maxBadPixels);


		//setValue(preprocessedProjections, transmissionProjections);

		//// Perform dead pixel correction
		//if (settings.removeDeadPixels)
		//{
		//	cout << "Dead pixel removal..." << endl;
		//	deadPixelRemoval(preprocessedProjections);
		//}

		//// Calculate -log or phase retrieval
		//if (settings.phaseMode == PhaseMode::Absorption)
		//	cout << "Negative logarithm..." << endl;
		//else
		//	cout << "Phase retrieval..." << endl;
		//phaseRetrieval(preprocessedProjections, settings.phaseMode, settings.phasePadType, settings.phasePadFraction, settings.objectCameraDistance, settings.delta, settings.mu);

		//cout << "Beam hardening correction..." << endl;
		//beamHardeningCorrection(preprocessedProjections, settings.bhc);

		//
		//cout << "Weighting..." << endl;
		//fbpWeighting(preprocessedProjections, settings.reconstructAs180degScan, settings.angles, settings.centerShift, settings.csAngleSlope, settings.sourceToRA, settings.cameraZShift, settings.csZSlope, centralAngle, gammamax0, settings.heuristicSinogramWindowingParameter);

		//cout << "Filter..." << endl;
		//filter(preprocessedProjections, settings.padType, settings.padFraction, settings.filterType);
	}



#if defined(USE_OPENCL)

	namespace internals
	{

		struct CLEnv
		{
			size_t globalMemSize;
			size_t maxAllocSize;
			Vec3c max3DImageSize;
			cl::Context context;
			cl::Device device;
			cl::CommandQueue queue;
			cl::Kernel kernel;
			cl::Image3D output;
			cl::Image3D temp;
		};

		/**
		Fills OpenCL backprojection output image with zeroes.
		*/
		void backprojectOpenCLReset(CLEnv& env, const Image<float32_t>& output)
		{
			// Create output image and fill it with zeroes
			cl::ImageFormat format(CL_INTENSITY, CL_FLOAT);
			env.output = cl::Image3D(env.context, CL_MEM_READ_WRITE, format, output.width(), output.height(), output.depth());
			env.temp = cl::Image3D(env.context, CL_MEM_READ_WRITE, format, output.width(), output.height(), output.depth());
			
			cl_uint4 zero;
			zero.s[0] = 0;
			zero.s[1] = 0;
			zero.s[2] = 0;
			zero.s[3] = 0;

			cl::size_t<3> outputSize;
			outputSize[0] = output.width();
			outputSize[1] = output.height();
			outputSize[2] = output.depth();

			env.queue.enqueueFillImage(env.temp, zero, cl::size_t<3>(), outputSize);
		}

		/**
		Initializes OpenCL system for backprojection.
		*/
		CLEnv backprojectOpenCLPrepare(const RecSettings& settings)
		{

			const char* newBackprojectProgram = R"***(

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

float csZPerturbation(float projectionZ, float projectionHeight, float csZSlope)
{
	return (projectionZ - projectionHeight / 2) * csZSlope;
}

kernel void backproject(read_only image3d_t transmissionProjections,
						global read_only float3* pss,				// Source positions
						global read_only float3* pds,				// Detector center point positions
						global read_only float3* us,				// Detector right vectors
						global read_only float3* vs,				// Detector up vectors
						global read_only float3* ws,				// Detector normal vectors, w = u x w

						read_only float2 projectionShift,
						read_only float2 fullProjectionHalfSize,	// Projection size divided by 2

						read_only float3 center,					// Vec3f(settings.roiSize) / 2.0f - Vec3f(settings.roiCenter) - Vec3f(0.5, 0.5, 0.5)
						read_only float sourceToRA,

						read_only image3d_t initialValue,
						write_only image3d_t output)
{

	const int3 projSize = (int3)(get_image_width(transmissionProjections), get_image_height(transmissionProjections), get_image_depth(transmissionProjections));
	const int3 outSize = (int3)(get_image_width(output), get_image_height(output), get_image_depth(output));

	float projectionHalfWidth = fullProjectionHalfSize.x;
	float projectionHalfHeight = fullProjectionHalfSize.y;
	
	
	// Calculate position of current pixel
    const int3 pos = (int3)(get_global_id(0), get_global_id(1), get_global_id(2));
    const float3 posf = convert_float3(pos);
    if (pos.x < 0 || pos.y < 0 || pos.z < 0 || pos.x >= outSize.x || pos.y >= outSize.y || pos.z >= outSize.z)
        return;

	const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;

	// Read initial value of the pixel as this kernel may be called multiple times with different set of projections if
	// all of them do not fit into memory at once.
	//const sampler_t readSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;
	//float sum = read_imagef(initialValue, readSampler, (int4)(pos.x, pos.y, pos.z, 0)).x;
	float sum = read_imagef(initialValue, (int4)(pos.x, pos.y, pos.z, 0)).x;
	//float sum = 0;


	float3 p = posf - center;
	for (int anglei = 0; anglei < projSize.z; anglei++)
	{
		float3 ps = pss[anglei];
		float3 pd = pds[anglei];
		float3 uVec = us[anglei];
		float3 vVec = vs[anglei];
		float3 wHat = ws[anglei];

		float3 dVec = p - ps;
		float denom = dot(dVec, wHat);

		if(fabs(denom) > 1e-6)
		{
			float3 psmpd = ps - pd; // NOTE: In principle we could pre-calculate ps-pd as we don't need plain pd anywhere.
			float d = (-dot(psmpd, wHat)) / denom;

			// (ps + d * dVec) is projection of source point through reconstruction point to the detector plane,
			// and pd is detector position.
			float3 pDot = psmpd + d * dVec;

			// Convert to camera coordinates
			float u = dot(pDot, uVec);
			float v = dot(pDot, vVec);

			// NOTE: We have pre-divided uVec and vVec by M, so division is not necessary here.
			// If uVec and vVec would be unit vectors, then we have to divide by M here.
			//u /= M;
			//v /= M;
 
			u += projectionHalfWidth;
			v += projectionHalfHeight;

			// Apply projection shift (we load only part of projections)
            u -= projectionShift.x;
			v -= projectionShift.y;

			// This weight is needed in the FDK algorithm.
			float weight = sourceToRA / (sourceToRA + dot(p, wHat));
			weight *= weight;

			float imgVal = read_imagef(transmissionProjections, sampler, (float4)(u + 0.5f, v + 0.5f, anglei + 0.5f, 0)).s0;
			sum += weight * imgVal;
		}
	}

	write_imagef(output, (int4)(pos.x, pos.y, pos.z, 0), sum);
};
)***";

//			const char* backprojectProgram = R"***(
//
//#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable
//
//float csZPerturbation(float z, float projectionHeight, float csZSlope)
//{
//	return (z - projectionHeight / 2) * csZSlope;
//}
//
//kernel void backproject(read_only image3d_t transmissionProjections,
//						global read_only float3* xHatArray,
//						global read_only float3* nHatArray,
//						read_only float centerShift, read_only float csZSlope,
//						global read_only float* csAnglePerturbations,
//						read_only float3 center,
//						read_only float d,
//						global read_only float2* objectShifts,
//						read_only float cameraZShift,
//						//read_only float dynMin, read_only float dynMax, read_only float scale,
//						read_only float2 projectionShift,
//						read_only float2 fullProjectionSize,
//						read_only image3d_t initialValue,
//						write_only image3d_t output)
//{
//
//	const int3 projSize = (int3)(get_image_width(transmissionProjections), get_image_height(transmissionProjections), get_image_depth(transmissionProjections));
//	const int3 outSize = (int3)(get_image_width(output), get_image_height(output), get_image_depth(output));
//
//	//float projectionWidth = (float)projSize.x;
//	//float projectionHeight = (float)projSize.y;
//	float projectionWidth = fullProjectionSize.x;
//	float projectionHeight = fullProjectionSize.y;
//	
//	
//	// Calculate position of current pixel
//    const int3 pos = (int3)(get_global_id(0), get_global_id(1), get_global_id(2));
//    const float3 posf = convert_float3(pos);
//    if (pos.x < 0 || pos.y < 0 || pos.z < 0 || pos.x >= outSize.x || pos.y >= outSize.y || pos.z >= outSize.z)
//        return;
//
//	float3 zHat = (float3)(0, 0, 1);
//
//	const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;
//
//	float currentCS = centerShift + csZPerturbation(posf.z, projectionHeight, csZSlope);
//
//	// Read initial value of the pixel as this kernel may be called multiple times with different set of projections if
//	// all of them do not fit into memory at once.
//	//const sampler_t readSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;
//	//float sum = read_imagef(initialValue, readSampler, (int4)(pos.x, pos.y, pos.z, 0)).x;
//	float sum = read_imagef(initialValue, (int4)(pos.x, pos.y, pos.z, 0)).x;
//	//float sum = 0;
//
//	for (int anglei = 0; anglei < projSize.z; anglei++)
//	{
//		float3 rho = posf - center;
//		float3 xHat = xHatArray[anglei];
//		float3 nHat = nHatArray[anglei];
//
//		float dprhox = d + dot(rho, xHat);
//		float Y = (d * dot(rho, nHat)) / dprhox;
//		float Z = (d * dot(rho, zHat)) / dprhox;
//
//		float w = (d * d) / (dprhox * dprhox);
//
//		// Get data from ideal detector position (Y, Z)
//		// First account for detector shifts and rotation
//		// TODO: Actually we have object shifts so this is only approximation that is correct for parallel beam case.
//		float sdx = objectShifts[anglei].x;
//		float sdz = objectShifts[anglei].y;
//		float angleCS = currentCS + csAnglePerturbations[anglei];
//		float ix = Y + projectionWidth / 2.0f + angleCS - sdx;
//		float iy = Z + projectionHeight / 2.0f + cameraZShift - sdz;
//
//		// TODO: Handle camera rotation here
//
//		// Apply projection shift (we load only part of projections)
//		ix -= projectionShift.x;
//		iy -= projectionShift.y;
//
//		float imgVal = read_imagef(transmissionProjections, sampler, (float4)(ix + 0.5f, iy + 0.5f, anglei + 0.5f, 0)).s0;
//
//		sum += w * imgVal;
//	}
//
//	// This is done later as this kernel may be called multiple times with different set of projections if
//	// all the projections do not fit into device memory at once.
//	//sum *= normFactor;
//	//sum = (sum - dynMin) / (dynMax - dynMin) * scale;
//
//	write_imagef(output, (int4)(pos.x, pos.y, pos.z, 0), sum);
//
//
//
//
//	//// Scaling
//	//sum = (sum - settings.dynMin) / (settings.dynMax - settings.dynMin) * NumberUtils<out_t>::scale();
//	////output(x, y, zi) = pixelRound<out_t>(sum);
//	//output(x, y, z) = pixelRound<out_t>(sum);
//
//};
//)***";

			try
			{

				std::vector<cl::Platform> platforms;
				cl::Platform::get(&platforms);
				if (platforms.size() == 0)
					throw ITLException("No OpenCL platforms available.");

				cout << platforms[0].getInfo<CL_PLATFORM_NAME>() << endl;
				cout << platforms[0].getInfo<CL_PLATFORM_VERSION>() << endl;

				cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0 };
				cl::Context context(CL_DEVICE_TYPE_GPU, properties);

				std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
				cl::Device device = devices[0];
				devices.clear();
				devices.push_back(device);

				string devName = device.getInfo<CL_DEVICE_NAME>();
				cout << "Using " << devName << endl;
				size_t globalMemSize = device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
				cout << "Global memory size " << bytesToString((double)globalMemSize) << endl;

				size_t maxAllocSize = device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
				cout << "Maximum size of image " << bytesToString((double)maxAllocSize) << endl;

				Vec2c max2DSize(device.getInfo<CL_DEVICE_IMAGE2D_MAX_WIDTH>(), device.getInfo<CL_DEVICE_IMAGE2D_MAX_HEIGHT>());
				cout << "Maximum linear 2D image size " << max2DSize << endl;
				Vec3c max3DSize(device.getInfo<CL_DEVICE_IMAGE3D_MAX_WIDTH>(), device.getInfo<CL_DEVICE_IMAGE3D_MAX_HEIGHT>(), device.getInfo<CL_DEVICE_IMAGE3D_MAX_DEPTH>());
				cout << "Maximum linear 3D image size " << max3DSize << endl;

				cout << "Maximum work group size " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << endl;
				cout << "Maximum work item dimensions " << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << endl;
				auto v = device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
				cout << "Maximum work item sizes ";
				for (auto x : v)
					cout << x << " ";
				cout << endl;

				string extensions = device.getInfo<CL_DEVICE_EXTENSIONS>();
				if (!contains(extensions, "cl_khr_3d_image_writes"))
                    cout << string("Warning: The OpenCL device ") + devName + string(" does not support cl_khr_3d_image_writes extension required by this program. On some devices the functionality is available even if the extension is not reported to be supported.") << endl;
					
                    //throw ITLException(string("The OpenCL device ") + devName + string(" does not support cl_khr_3d_image_writes extension required by this program. List of all supported extensions: ") + extensions);

				cl::Program::Sources source(1, std::make_pair(newBackprojectProgram, strlen(newBackprojectProgram)));
				cl::Program program = cl::Program(context, source);
				try
				{
					program.build(devices);
				}
				catch (cl::Error err)
				{
					string msg = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
					throw ITLException("OpenCL program build failed.\n" + msg);
				}

				cl_int err;
				cl::Kernel kernel(program, "backproject", &err);

				cl::CommandQueue queue(context, device, 0, &err);

				CLEnv env = { globalMemSize, maxAllocSize, max3DSize, context, device, queue, kernel, cl::Image3D() };

				return env;
			}
			catch (cl::Error err)
			{
				throw ITLException("Unable to initialize OpenCL. Error at " + string(err.what()) + ", error code = " + toString(err.err()));
			}
		}

		/**
		Internal OpenCL backprojection routine.
		@param firstAngleIndex Index of angle in settings.angles corresponding to the first image in transmissionProjections.
		*/
		void backprojectOpenCL(const Image<float32_t>& transmissionProjections, const Vec3c& fullProjectionSize, const Vec3c& projBlockStart,
			//const RecSettings& settings,
			std::vector<Vec4f> pss,			// Source positions
			std::vector<Vec4f> pds,			// Detector center point positions
			std::vector<Vec4f> us,			// Detector right vectors
			std::vector<Vec4f> vs,			// Detector up vectors
			std::vector<Vec4f> ws,			// Detector normal vectors, w = u x w
			Vec3f roiSize, Vec3f roiCenter, float32_t sourceToRA,
			Image<float32_t>& output, CLEnv& clEnv)
		{

			

			Vec3f center = roiSize / 2.0f - roiCenter - Vec3f(0.5, 0.5, 0.5);

			Vec2f fullProjectionHalfSizeFloat((float32_t)fullProjectionSize.x / 2.0f, (float32_t)fullProjectionSize.y / 2.0f);
			Vec2f projectionShift((float32_t)projBlockStart.x, (float32_t)projBlockStart.y);

			Timer timer;

			try
			{
				cl::ImageFormat format(CL_INTENSITY, CL_FLOAT);
				cl::Image3D transmissionProjectionsCL = cl::Image3D(clEnv.context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, format, transmissionProjections.width(), transmissionProjections.height(), transmissionProjections.depth());

				cl::Buffer pssCL = cl::Buffer(clEnv.context, pss.begin(), pss.end(), true, true);
				cl::Buffer pdsCL = cl::Buffer(clEnv.context, pds.begin(), pds.end(), true, true);
				cl::Buffer usCL = cl::Buffer(clEnv.context, us.begin(), us.end(), true, true);
				cl::Buffer vsCL = cl::Buffer(clEnv.context, vs.begin(), vs.end(), true, true);
				cl::Buffer wsCL = cl::Buffer(clEnv.context, ws.begin(), ws.end(), true, true);

				Vec4f centerCL(center.x, center.y, center.z, 0); // CL defines float3 as float4, so we need extra zero after center position.

				// Set all arguments
				clEnv.kernel.setArg(0, transmissionProjectionsCL);
				clEnv.kernel.setArg(1, pssCL);
				clEnv.kernel.setArg(2, pdsCL);
				clEnv.kernel.setArg(3, usCL);
				clEnv.kernel.setArg(4, vsCL);
				clEnv.kernel.setArg(5, wsCL);
				clEnv.kernel.setArg(6, projectionShift);
				clEnv.kernel.setArg(7, fullProjectionHalfSizeFloat);
				clEnv.kernel.setArg(8, centerCL);
				clEnv.kernel.setArg(9, sourceToRA);
				clEnv.kernel.setArg(10, clEnv.temp);
				clEnv.kernel.setArg(11, clEnv.output);


				cl::size_t<3> projSizeCL;
				projSizeCL[0] = transmissionProjections.width();
				projSizeCL[1] = transmissionProjections.height();
				projSizeCL[2] = transmissionProjections.depth();
				clEnv.queue.enqueueWriteImage(transmissionProjectionsCL, true, cl::size_t<3>(), projSizeCL, 0, 0, (void*)transmissionProjections.getData());


				cl::size_t<3> outSizeCL;
				outSizeCL[0] = output.width();
				outSizeCL[1] = output.height();
				outSizeCL[2] = output.depth();
				clEnv.queue.enqueueCopyImage(clEnv.output, clEnv.temp, cl::size_t<3>(), cl::size_t<3>(), outSizeCL);

				clEnv.queue.finish();

				// Enqueue in batches so that watchdog timer issues are avoided
				size_t blockSize = 1;
				coord_t zStart = 0;
				do
				{
					coord_t zEnd = zStart + blockSize;
					if (zEnd > output.depth())
						zEnd = output.depth();

					//cout << "enqueueNDRangeKernel z = " << zStart << " - " << zEnd << endl;

					timer.start();
					clEnv.queue.enqueueNDRangeKernel(clEnv.kernel,
						cl::NDRange(0, 0, zStart),
						cl::NDRange(output.width(), output.height(), zEnd - zStart),
						cl::NullRange);

					clEnv.queue.finish();
					timer.stop();

					zStart = zEnd;

					// Auto-tune block size
					if (timer.getSeconds() < 0.9)
						blockSize++;
					else if (timer.getSeconds() > 1.1)
						blockSize--;

					if (blockSize < 1)
						blockSize = 1;


				} while (zStart < output.depth());

			}
			catch (cl::Error err)
			{
				timer.stop();
				//cout << timer.getSeconds() << endl;
				throw ITLException("OpenCL error at " + string(err.what()) + ", error code = " + toString(err.err()) + ", OpenCL kernel run time = " + toString(timer.getSeconds()) + " s.");
			}

		}

		///**
		//Internal OpenCL backprojection routine.
		//@param firstAngleIndex Index of angle in settings.angles corresponding to the first image in transmissionProjections.
		//*/
		//void backprojectOpenCL(const Image<float32_t>& transmissionProjections, const Vec3c& fullProjectionSize, const Vec3c& projBlockStart, const RecSettings& settings, Image<float32_t>& output, CLEnv& clEnv)
		//{

		//	float32_t gammamax0 = internals::calculateGammaMax0((float32_t)transmissionProjections.width(), settings.sourceToRA);
		//	float32_t centralAngle = internals::calculateTrueCentralAngle(settings.centralAngleFor180degScan, settings.angles, gammamax0);


		//	// Pre-calculate direction vectors. OpenCL requires 4-component vectors!
		//	std::vector<Vec4f> xHatArray;
		//	std::vector<Vec4f> nHatArray;
		//	xHatArray.reserve(transmissionProjections.depth());
		//	nHatArray.reserve(transmissionProjections.depth());

		//	// Pre-calculate center shift angle perturbations
		//	vector<float32_t> csAnglePerturbations;
		//	csAnglePerturbations.reserve(transmissionProjections.depth());

		//	// Make arrays containing correct shifts
		//	vector<Vec2f> objectShifts;
		//	objectShifts.reserve(transmissionProjections.depth());

		//	double rotMul = 1.0;
		//	if (settings.rotationDirection == RotationDirection::Counterclockwise)
		//		rotMul = -1.0;

		//	for (size_t anglei = projBlockStart.z; anglei < projBlockStart.z + (size_t)transmissionProjections.depth(); anglei++)
		//	{
		//		double angle = rotMul * (settings.angles[anglei] - 90 + settings.rotation) / 180.0 * PI;

		//		float32_t c = (float32_t)cos(angle);
		//		float32_t s = (float32_t)sin(angle);
		//		xHatArray.push_back(Vec4f(c, s, 0, 0));
		//		nHatArray.push_back(Vec4f(-s, c, 0, 0));

		//		csAnglePerturbations.push_back(internals::csAnglePerturbation(anglei, centralAngle, settings.angles, settings.csAngleSlope));

		//		objectShifts.push_back(settings.objectShifts[anglei] * settings.shiftScaling * (settings.useShifts ? 1.0f : 0.0f));
		//	}
		//	Vec3f zHat(0, 0, 1);

		//	// Backproject
		//	float32_t d = settings.sourceToRA;

		//	
		//	// NOTE: -0.5 ensures that if roiSize.z == 1, roiCenter.z = 1 / 2 - roiCenter.z - 0.5 = roiCenter.z
		//	// TODO: Should subtraction be done for all the components?
		//	Vec3f center = Vec3f(settings.roiSize) / 2.0f - Vec3f(settings.roiCenter) - Vec3f(0, 0, 0.5);

		//	Vec2f fullProjectionSizeFloat((float32_t)fullProjectionSize.x, (float32_t)fullProjectionSize.y);
		//	Vec2f projectionShift((float32_t)projBlockStart.x, (float32_t)projBlockStart.y);

		//	Timer timer;

		//	try
		//	{
		//		cl::ImageFormat format(CL_INTENSITY, CL_FLOAT);
		//		cl::Image3D transmissionProjectionsCL = cl::Image3D(clEnv.context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, format, transmissionProjections.width(), transmissionProjections.height(), transmissionProjections.depth());

		//		cl::Buffer xHatArrayCL = cl::Buffer(clEnv.context, xHatArray.begin(), xHatArray.end(), true, true);
		//		cl::Buffer nHatArrayCL = cl::Buffer(clEnv.context, nHatArray.begin(), nHatArray.end(), true, true);
		//		cl::Buffer csAnglePerturbationsCL = cl::Buffer(clEnv.context, csAnglePerturbations.begin(), csAnglePerturbations.end(), true, true);

		//		Vec4f centerCL(center.x, center.y, center.z, 0); // CL defines float3 as float4, so we need extra zero after center position.

		//		cl::Buffer objectShiftsCL = cl::Buffer(clEnv.context, objectShifts.begin(), objectShifts.end(), true, true);

		//		// Set all arguments
		//		clEnv.kernel.setArg(0, transmissionProjectionsCL);
		//		clEnv.kernel.setArg(1, xHatArrayCL);
		//		clEnv.kernel.setArg(2, nHatArrayCL);
		//		clEnv.kernel.setArg(3, settings.centerShift);
		//		clEnv.kernel.setArg(4, settings.csZSlope);
		//		clEnv.kernel.setArg(5, csAnglePerturbationsCL);
		//		clEnv.kernel.setArg(6, centerCL);
		//		clEnv.kernel.setArg(7, d);
		//		clEnv.kernel.setArg(8, objectShiftsCL);
		//		clEnv.kernel.setArg(9, settings.cameraZShift);
		//		//clEnv.kernel.setArg(10, settings.dynMin);
		//		//clEnv.kernel.setArg(11, settings.dynMax);
		//		//clEnv.kernel.setArg(12, NumberUtils<float32_t>::scale());
		//		//clEnv.kernel.setArg(13, outputCL);
		//		clEnv.kernel.setArg(10, projectionShift); // Shift that must be applied to projection coordinates
		//		clEnv.kernel.setArg(11, fullProjectionSizeFloat); // Width and height of full projection images
		//		clEnv.kernel.setArg(12, clEnv.temp);
		//		clEnv.kernel.setArg(13, clEnv.output);


		//		cl::size_t<3> projSizeCL;
		//		projSizeCL[0] = transmissionProjections.width();
		//		projSizeCL[1] = transmissionProjections.height();
		//		projSizeCL[2] = transmissionProjections.depth();
		//		clEnv.queue.enqueueWriteImage(transmissionProjectionsCL, true, cl::size_t<3>(), projSizeCL, 0, 0, (void*)transmissionProjections.getData());
		//		

		//		cl::size_t<3> outSizeCL;
		//		outSizeCL[0] = output.width();
		//		outSizeCL[1] = output.height();
		//		outSizeCL[2] = output.depth();
		//		clEnv.queue.enqueueCopyImage(clEnv.output, clEnv.temp, cl::size_t<3>(), cl::size_t<3>(), outSizeCL);

		//		clEnv.queue.finish();

		//		// Enqueue in batches so that watchdog timer issues are avoided
		//		size_t blockSize = 1; // TODO: Optimize block size somehow.
		//		coord_t zStart = 0;
		//		do
		//		{
		//			coord_t zEnd = zStart + blockSize;
		//			if (zEnd > output.depth())
		//				zEnd = output.depth();
		//			
		//			//cout << "enqueueNDRangeKernel z = " << zStart << " - " << zEnd << endl;
		//			
		//			timer.start();
		//			clEnv.queue.enqueueNDRangeKernel(clEnv.kernel,
		//				cl::NDRange(0, 0, zStart),
		//				cl::NDRange(output.width(), output.height(), zEnd - zStart),
		//				cl::NullRange);

		//			clEnv.queue.finish();

		//			zStart = zEnd;
		//		} while (zStart < output.depth());

		//	}
		//	catch (cl::Error err)
		//	{
		//		timer.stop();
		//		//cout << timer.getSeconds() << endl;
		//		throw ITLException("OpenCL error at " + string(err.what()) + ", error code = " + toString(err.err()) + ", OpenCL kernel run time = " + toString(timer.getSeconds()) + " s.");
		//	}

		//}

		/**
		Performs scaling of output pixel values to desired range.
		*/
		void backprojectOpenCLFinalize(const RecSettings& settings, Image<float32_t>& output, CLEnv& clEnv)
		{
			cl::size_t<3> fullSize;
			fullSize[0] = output.width();
			fullSize[1] = output.height();
			fullSize[2] = output.depth();
			clEnv.queue.enqueueReadImage(clEnv.output, true, cl::size_t<3>(), fullSize, 0, 0, (void*)output.getData());

			float32_t normFact = normFactor(settings);

			// TODO: This normalization process can be done in OpenCL, but it is also quite fast if done on CPU.
			#pragma omp parallel for
			for (coord_t n = 0; n < output.pixelCount(); n++)
				output(n) = ((output(n) * normFact) - settings.dynMin) / (settings.dynMax - settings.dynMin) * NumberUtils<float32_t>::scale();
		}
	}

	/**
	OpenCL backprojection without any special treatment for large input and output arrays.
	*/
	void backprojectOpenCLSimple(Image<float32_t>& transmissionProjections, RecSettings settings, Image<float32_t>& output)
	{
		internals::sanityCheck(transmissionProjections, settings, true);
		output.mustNotBe(transmissionProjections);
		
		internals::applyBinningToParameters(settings);

		output.ensureSize(settings.roiSize);

		// Pre-calculate direction vectors. OpenCL requires 4-component vectors!
		std::vector<Vec4f> pss;			// Source positions
		std::vector<Vec4f> pds;			// Detector center point positions
		std::vector<Vec4f> us;			// Detector right vectors
		std::vector<Vec4f> vs;			// Detector up vectors
		std::vector<Vec4f> ws;			// Detector normal vectors, w = u x w
		internals::determineBackprojectionGeometry(settings, transmissionProjections.width(), pss, pds, us, vs, ws);

		internals::CLEnv env = internals::backprojectOpenCLPrepare(settings);
		internals::backprojectOpenCLReset(env, output);
		internals::backprojectOpenCL(transmissionProjections, transmissionProjections.dimensions(), Vec3c(0, 0, 0),
			pss, pds, us, vs, ws,
			Vec3f(settings.roiSize), Vec3f(settings.roiCenter), settings.sourceToRA,
			output, env);
		internals::backprojectOpenCLFinalize(settings, output, env);
	}

	namespace internals
	{
		void backprojectOpenCLProjectionBlocks(const Image<float32_t>& transmissionProjections, const RecSettings& settings, Image<float32_t>& output, CLEnv& env, coord_t maxProjBlockSizeZ = numeric_limits<coord_t>::max())
		{
			// Memory problems:
			// Projection data size > device memory size
			// Output data size > device memory size
			// Divide output to ROIs and process each ROI separately
			// In this function reconstruct one ROI



			// Pre-calculate direction vectors. OpenCL requires 4-component vectors!
			std::vector<Vec4f> pss;			// Source positions
			std::vector<Vec4f> pds;			// Detector center point positions
			std::vector<Vec4f> us;			// Detector right vectors
			std::vector<Vec4f> vs;			// Detector up vectors
			std::vector<Vec4f> ws;			// Detector normal vectors, w = u x w
			internals::determineBackprojectionGeometry(settings, transmissionProjections.width(), pss, pds, us, vs, ws);



			// Calculate positions of all 8 corners of the ROI
			// Project all 8 ROI corners to the projection plane in all angles,
			// select minimum and maximum coordinates in Y and Z.
			// In the angle direction always load all the projections as
			// those are needed anyway.
			float32_t w = (float32_t)settings.roiSize.x;
			float32_t h = (float32_t)settings.roiSize.y;
			float32_t d = (float32_t)settings.roiSize.z;
			array<Vec3f, 8> corners =
				{
					Vec3f(0, 0, 0),
					Vec3f(w, 0, 0),
					Vec3f(0, h, 0),
					Vec3f(w, h, 0),
					Vec3f(0, 0, d),
					Vec3f(w, 0, d),
					Vec3f(0, h, d),
					Vec3f(w, h, d),
				};


			Vec2f projMinCoords(numeric_limits<float32_t>::infinity(), numeric_limits<float32_t>::infinity());
			Vec2f projMaxCoords = -projMinCoords;

			float32_t projectionHalfWidth = (float32_t)transmissionProjections.width() / 2.0f;
			float32_t projectionHalfHeight = (float32_t)transmissionProjections.height() / 2.0f;

			for (size_t m = 0; m < corners.size(); m++)
			{
				Vec3f p = corners[m] - Vec3f(settings.roiSize) / 2.0f + Vec3f(settings.roiCenter)
					//+ Vec3f(0, 0, 0.5); // Old version was like this
					+ Vec3f(0.5, 0.5, 0.5);

				float32_t sum = 0;
				for (coord_t anglei = 0; anglei < transmissionProjections.depth(); anglei++)
				{
					Vec3f ps = Vec3f(pss[anglei]);
					Vec3f pd = Vec3f(pds[anglei]);
					Vec3f uVec = Vec3f(us[anglei]);
					Vec3f vVec = Vec3f(vs[anglei]);
					Vec3f wHat = Vec3f(ws[anglei]);

					Vec3f dVec = p - ps;
					float32_t denom = dVec.dot(wHat);
					if (!NumberUtils<float32_t>::equals(denom, 0))
					{
						Vec3f psmpd = ps - pd; // NOTE: In principle we could pre-calculate ps-pd as we don't need plain pd anywhere.
						float32_t d = (-psmpd.dot(wHat)) / denom;

						// (ps + d * dVec) is projection of source point through reconstruction point to the detector plane,
						// and pd is detector position.
						Vec3f pDot = psmpd + d * dVec;

						// Convert to camera coordinates
						float32_t u = pDot.dot(uVec);
						float32_t v = pDot.dot(vVec);

						// NOTE: We have pre-divided uVec and vVec by M, so division is not necessary here.
						// If uVec and vVec would be unit vectors, then we have to divide by M here.
						//u /= M;
						//v /= M;

						u += projectionHalfWidth;
						v += projectionHalfHeight;

						projMinCoords = min(projMinCoords, Vec2f(u, v));
						projMaxCoords = max(projMaxCoords, Vec2f(u, v));
					}
				}
			}



			//float32_t gammamax0 = internals::calculateGammaMax0((float32_t)transmissionProjections.width(), settings.sourceToRA);
			//float32_t centralAngle = internals::calculateTrueCentralAngle(settings.centralAngleFor180degScan, settings.angles, gammamax0);

			//double rotMul = 1.0;
			//if (settings.rotationDirection == RotationDirection::Counterclockwise)
			//	rotMul = -1.0;

			//float32_t projectionWidth = (float32_t)transmissionProjections.width();
			//float32_t projectionHeight = (float32_t)transmissionProjections.height();

			//Vec3f center = Vec3f(settings.roiSize) / 2.0f - Vec3f(settings.roiCenter);

			//Vec2f projMinCoords(numeric_limits<float32_t>::infinity(), numeric_limits<float32_t>::infinity());
			//Vec2f projMaxCoords = -projMinCoords;
			//for (size_t anglei = 0; anglei < settings.angles.size(); anglei++)
			//{
			//	double angle = rotMul * (settings.angles[anglei] - 90 + settings.rotation) / 180.0 * PI;
			//	
			//	float32_t c = (float32_t)cos(angle);
			//	float32_t s = (float32_t)sin(angle);

			//	Vec3f xHat(c, s, 0);
			//	Vec3f nHat(-s, c, 0);
			//	Vec3f zHat(0, 0, 1);
			//	float32_t d = settings.sourceToRA;
			//	
			//	for (size_t m = 0; m < corners.size(); m++)
			//	{
			//		Vec3f xyz = corners[m];

			//		float32_t currentCS = settings.centerShift + internals::csZPerturbation(xyz.z, (float32_t)settings.roiCenter.z, (float32_t)settings.roiSize.z, projectionHeight, settings.csZSlope);

			//		Vec3f rho = xyz - center;

			//		float32_t dprhox = d + rho.dot(xHat);
			//		float32_t Y = (d * rho.dot(nHat)) / dprhox;
			//		float32_t Z = (d * rho.dot(zHat)) / dprhox;

			//		float32_t sdx = settings.objectShifts[anglei].x * settings.shiftScaling * (settings.useShifts ? 1 : 0);
			//		float32_t sdz = settings.objectShifts[anglei].y * settings.shiftScaling * (settings.useShifts ? 1 : 0);
			//		float32_t angleCS = currentCS + internals::csAnglePerturbation(anglei, centralAngle, settings.angles, settings.csAngleSlope);
			//		float32_t ix = Y + projectionWidth / 2.0f + angleCS - sdx;
			//		float32_t iy = Z + projectionHeight / 2.0f + settings.cameraZShift - sdz;

			//		projMinCoords = min(projMinCoords, Vec2f(ix, iy));
			//		projMaxCoords = max(projMaxCoords, Vec2f(ix, iy));
			//	}
			//}

			Vec2c projMin = floor(projMinCoords - Vec2f(2, 2));
			Vec2c projMax = ceil(projMaxCoords + Vec2f(2, 2) + Vec2f(1, 1));

			projMin = max(projMin, Vec2c(0, 0));
			projMax = min(projMax, Vec2c(transmissionProjections.width(), transmissionProjections.height()));

			Vec2c projBlockSize2D = projMax - projMin;

			Vec3c projBlockSize = env.max3DImageSize;

			if (projBlockSize.x < projBlockSize2D.x ||
				projBlockSize.y < projBlockSize2D.y)
				throw ITLException("OpenCL device does not support large enough 3D images to hold the required block of projection image.");

			projBlockSize.x = projBlockSize2D.x;
			projBlockSize.y = projBlockSize2D.y;

			projBlockSize.z = std::min(projBlockSize.z, transmissionProjections.depth());
			projBlockSize.z = std::min(projBlockSize.z, maxProjBlockSizeZ);

			size_t outputSize = output.pixelCount() * sizeof(float32_t);
			size_t allocSize = projBlockSize.x * projBlockSize.y * projBlockSize.z * sizeof(float32_t);
			while (outputSize + allocSize >= 0.95 * env.globalMemSize ||
				allocSize >= env.maxAllocSize)
			{
				projBlockSize.z--;
				allocSize = projBlockSize.x * projBlockSize.y * projBlockSize.z * sizeof(float32_t);
			}

			cout << "Dividing projection images into blocks of size " << projBlockSize << endl;


			// Copy data to new image (TODO: Remove this step somehow, perhaps improving the 'view' capabilities. Right now that
			// is not done as a simple addition of row and slice stride breaks capability for simple linear loop over all pixels.)
			Image<float32_t> projectionBlockData(Vec3c(projBlockSize.x, projBlockSize.y, transmissionProjections.depth()));
			crop(transmissionProjections, projectionBlockData, Vec3c(projMin.x, projMin.y, 0));


			coord_t zStart = 0;
			do
			{
				coord_t zEnd = zStart + projBlockSize.z - 1;
				if (zEnd >= projectionBlockData.depth())
					zEnd = projectionBlockData.depth() - 1;
				Image<float32_t> projBlock(projectionBlockData, zStart, zEnd);

				cout << "Backprojecting projections " << zStart << " - " << zEnd << endl;

				// Backproject this block of projections
				internals::backprojectOpenCL(projBlock, transmissionProjections.dimensions(), Vec3c(projMin.x, projMin.y, zStart), 
					pss, pds, us, vs, ws,
					Vec3f(settings.roiSize), Vec3f(settings.roiCenter), settings.sourceToRA,
					output, env);

				zStart = zEnd + 1;
			} while (zStart < projectionBlockData.depth());


			internals::backprojectOpenCLFinalize(settings, output, env);
		}
	}

	/**
	OpenCL backprojection that divides transmission projections into blocks that fit into the device memory.
	Algorithm:
	Take as many projections as fits into the memory in the device,
	backproject them, take next set of projections, backproject again, etc.
	Continue until all projections have been backprojected.
	Run division by 2pi, scaling, and conversion to output pixel format after backprojections.
	*/
	void backprojectOpenCLProjectionBlocks(Image<float32_t>& transmissionProjections, RecSettings settings, Image<float32_t>& output, coord_t maxProjBlockSizeZ = numeric_limits<coord_t>::max())
	{
		// Memory problems:
		// Projection data size > device memory size
		// Output data size > device memory size
		// Divide output to ROIs and process each ROI separately
		// In this function reconstruct one ROI

		internals::sanityCheck(transmissionProjections, settings, true);
		output.mustNotBe(transmissionProjections);
		
		internals::applyBinningToParameters(settings);
		
		output.ensureSize(settings.roiSize);

		try
		{
			internals::CLEnv env = internals::backprojectOpenCLPrepare(settings);
			internals::backprojectOpenCLReset(env, output);
			internals::backprojectOpenCLProjectionBlocks(transmissionProjections, settings, output, env, maxProjBlockSizeZ);
		}
		catch (cl::Error err)
		{
			throw ITLException("OpenCL error at " + string(err.what()) + ", error code = " + toString(err.err()));
		}
	}

	/**
	OpenCL backprojection that divides both projection and output data to blocks.
	*/
	void backprojectOpenCL(const Image<float32_t>& transmissionProjections, RecSettings settings, Image<float32_t>& output, coord_t maxOutputBlockSizeZ, coord_t maxProjectionBlockSizeZ)
	{
		internals::sanityCheck(transmissionProjections, settings, true);
		output.mustNotBe(transmissionProjections);
		
		internals::applyBinningToParameters(settings);
		
		output.ensureSize(settings.roiSize);

		try
		{
			// Divide output into blocks that fit well into device memory with projection data
			internals::CLEnv env = internals::backprojectOpenCLPrepare(settings);

			Vec3c blockSize = env.max3DImageSize;

			if (blockSize.x < output.width() ||
				blockSize.y < output.height())
				throw ITLException("OpenCL device does not support large enough 3D images to hold a single slice of the result image.");

			blockSize.x = output.width();
			blockSize.y = output.height();

			blockSize.z = std::min(blockSize.z, output.depth());
			blockSize.z = std::min(blockSize.z, maxOutputBlockSizeZ);

			while (blockSize.x * blockSize.y * blockSize.z * sizeof(float32_t) > std::min((size_t)round(0.45 * env.globalMemSize), env.maxAllocSize)) // TODO: The block size condition could probably be perfected...
				blockSize.z--;
			
			cout << "Dividing output to blocks of size " << blockSize << endl;

			Vec3c origCenter = settings.roiCenter;

			// Call backprojectOpenCLBlocks for each of the blocks
			coord_t zStart = 0;
			do
			{
				coord_t zEnd = zStart + blockSize.z - 1;
				if (zEnd >= output.depth())
					zEnd = output.depth() - 1;

				cout << "Reconstructing block of output image in z range " << zStart << " - " << zEnd << endl;

				Image<float32_t> outputBlock(output, zStart, zEnd);

				// Backproject this block of projections. Override ROI in the settings to reconstruct
				// the correct location.
				settings.roiCenter.z = origCenter.z - output.depth() / 2 + zStart + outputBlock.depth() / 2;
				settings.roiSize = outputBlock.dimensions();

				internals::backprojectOpenCLReset(env, outputBlock);
				internals::backprojectOpenCLProjectionBlocks(transmissionProjections, settings, outputBlock, env, maxProjectionBlockSizeZ);

				zStart = zEnd + 1;
			} while (zStart < output.depth());

		}
		catch (cl::Error err)
		{
			throw ITLException("OpenCL error at " + string(err.what()) + ", error code = " + toString(err.err()));
		}
	}


#endif



	namespace tests
	{
		void recSettings()
		{
		
			RecSettings settings;
        
			settings.angles.push_back(7);
			settings.sampleShifts.push_back(Vec3f(1, 2, 3));
			settings.sampleShifts.push_back(Vec3f(4, 5, 6));
			
			stringstream s;
			s << settings;

			cout << settings << endl;
			
			RecSettings settings2 = fromString<RecSettings>(s.str());
		
			
			stringstream s2;
			s2 << settings;

			testAssert(s.str() == s2.str(), "rec settings save/load");
		}


		void paganin()
		{

			Image<float32_t> slice(329, 33, 1);
			draw(slice, AABox<coord_t>(Vec3c(50, 10, 0), Vec3c(250, 25, 1)), 0.95f);

			raw::writed(slice, "./fbp/before_paganin");
			paganinSlice(slice, PadType::Mirror, 0, 0, 1000, 1, 1);
			raw::writed(slice, "./fbp/after_paganin");
		}

		void fbp()
		{
			Image<float32_t> original;
			raw::read(original, "projections/plates");

			// Generate weirdly shifted projections

			//raw::read(projections, "projections/projected120");
			coord_t projectionCount = 120;
			coord_t projWidth = 100;
			coord_t projHeight = 100;
			coord_t projPixelSize = 1;
			Vec2f angularRange(-180, 180);
			float32_t sourceDistance = 300;
			float32_t cameraDistance = 100;
			float32_t cameraPixelSize = 1 + cameraDistance / sourceDistance;

			Vec3f rotAxis(0, 0, 1);

			Image<float32_t> projections(projWidth, projHeight, projectionCount);
			vector<float32_t> angles;
			vector<Vec3f> sampleShifts;
			vector<Vec3f> raShifts;
			vector<Vec3f> sourceShifts;
			vector<Vec3f> cameraShifts;
			for (coord_t anglei = 0; anglei < projectionCount; anglei++)
			{
				float32_t rotationAngle = projectionCount > 1 ? angularRange.x + (angularRange.y - angularRange.x) / (projectionCount - 1) * anglei : angularRange.x;

				Vec3f srcPos(-sourceDistance, 0, 0);
				Vec3f cameraPos(cameraDistance, 0, 0);
				Vec3f u(0, 1, 0);
				Vec3f v(0, 0, 1);

				// rotation axis alignment stage position (stage below rotation axis, moves rotation axis and sample)
				float32_t tmp1 = 3 * PIf * (float32_t)anglei / projectionCount;
				float32_t A1 = 6;
				Vec3f raPos(A1 * sin(tmp1), A1 * sin(tmp1), A1 * sin(tmp1));
				srcPos -= raPos;
				cameraPos -= raPos;

				// fine-alignment stage position (stage on top of rotation axis, moves sample after rotation)
				float32_t tmp = 2 * PIf * (float32_t)anglei / projectionCount;
				float32_t A = 5;
				Vec3f sampleShift(A*cos(tmp), A*sin(tmp), A * cos(tmp));

				// Camera shift
				float32_t tmp3 = 5 * PIf * (float32_t)anglei / projectionCount;
				float32_t A3 = 9;
				Vec3f cameraShift(A3 * -sin(tmp3), A3 * cos(tmp3), A3 * -sin(tmp3));
				cameraPos += cameraShift;

				// Source shift
				float32_t tmp4 = 6 * PIf * (float32_t)anglei / projectionCount;
				float32_t A4 = 12;
				Vec3f srcShift(10 * A4 * -cos(tmp4), A4 * cos(tmp4), A4 * -cos(tmp4));
				srcPos += srcShift;



				float32_t angleRad = rotationAngle / 180.0f * PIf;
				itl2::rotate<float32_t>(srcPos, angleRad, 0.0f);
				itl2::rotate<float32_t>(cameraPos, angleRad, 0.0f);
				itl2::rotate<float32_t>(u, angleRad, 0.0f);
				itl2::rotate<float32_t>(v, angleRad, 0.0f);
				
				u *= cameraPixelSize;
				v *= cameraPixelSize;

				Image<float32_t> proj(projWidth, projHeight);
				itl2::siddonProject(original, srcPos, cameraPos, sampleShift, u, v, proj);
				copyValues(projections, proj, Vec3c(0, 0, anglei));

				angles.push_back(rotationAngle);
				sampleShifts.push_back(sampleShift);
				raShifts.push_back(raPos);
				cameraShifts.push_back(cameraShift);
				sourceShifts.push_back(srcShift);
			}

			raw::writed(projections, "fbp/projections");


			// Here projections = integrals of mu where mu=1 for objects
			// Convert it to I/I0, with mu = 0.01 for objects
			multiply(projections, 1e-2);
			negate(projections);
			exponentiate(projections);

			// Rotate projections -> camera rotation
			float32_t cameraRotation = 1.0f;
			Image<float32_t> temp(projections.dimensions());
			itl2::rotate(projections, temp, cameraRotation / 180.0 * 3.14, Vec3d(0, 0, 1), LinearInterpolator<float32_t, float32_t>(BoundaryCondition::Nearest));
			itl2::setValue(projections, temp);

			// Translate in z and x -> camera z shift and center shift
			float32_t centerShift = 5;
			float32_t zShift = 8;
			itl2::translate(projections, temp, Vec3d(centerShift, zShift, 0), LinearInterpolator<float32_t, float32_t>(BoundaryCondition::Nearest));
			itl2::setValue(projections, temp);

			raw::writed(projections, "fbp/transmission");



			RecSettings settings;
			settings.angles = angles;
			settings.sampleShifts = sampleShifts;
			settings.rotationAxisShifts = raShifts;
			settings.cameraShifts = cameraShifts;
			settings.sourceShifts = sourceShifts;
			settings.rotationDirection = RotationDirection::Clockwise;
			settings.filterType = FilterType::Ramp;
			settings.sourceToRA = sourceDistance;
			settings.objectCameraDistance = cameraDistance;
			settings.reconstructAs180degScan = false;

			// Above we shift the images, so that means the camera will shift to the inverse direction.
			settings.centerShift = centerShift;
			settings.cameraZShift = zShift;

			// Here we need minus sign in front of angle as the detector normal in the reconstruction is inverse of what
			// we used above to rotate the images.
			settings.cameraRotation = -cameraRotation;

			Image<float32_t> preProcProjections;
			itl2::fbpPreprocess(projections, preProcProjections, settings);
			

#if defined(USE_OPENCL)
			Image<float32_t> outputCL;
			itl2::backprojectOpenCL(preProcProjections, settings, outputCL);
			raw::writed(outputCL, "./fbp/reconstructed_cl");
#endif

			Image<float32_t> outputNew;
			itl2::backproject(preProcProjections, settings, outputNew);
			raw::writed(outputNew, "./fbp/reconstructed_new");

			//Image<float32_t> output;
			//itl2::backproject(preProcProjections, settings, output);
			//raw::writed(output, "./fbp/reconstructed");




			//Image<float32_t> output2;
			//reslice(output, output2, ResliceDirection::Top);
			//raw::writed(output2, "./fbp/reconstructed_resliced");

			//Image<float32_t> output2New;
			//reslice(outputNew, output2New, ResliceDirection::Top);
			//raw::writed(output2New, "./fbp/reconstructed_new_resliced");
		}

		void openCLBackProjection()
		{
#if defined(USE_OPENCL)
			// TODO: Test data is not available.

			string logFile = readText("./rec/fibre_log.txt", true);
			RecSettings settings = fromString<RecSettings>(logFile);

			Image<float32_t> projections;
			Image<float32_t> preProcProjections;
			raw::read(projections, "./rec/fibre_projections");

			Timer timer;

			Image<float32_t> output;
			Image<float32_t> output2;

			timer.start();
			itl2::fbpPreprocess(projections, preProcProjections, settings);
			timer.stop();
			cout << "Preprocessing took " << timer.getSeconds() << " s." << endl;

			cout << "OpenCL backprojection..." << endl;
			timer.start();
			itl2::backprojectOpenCLSimple(preProcProjections, settings, output);
			timer.stop();
			cout << "Backprojection took " << timer.getSeconds() << " s." << endl;

			raw::writed(output, "./rec/fibre_reconstruction_mid_gpu");



			cout << "OpenCL backprojection, projections in blocks" << endl;

			timer.start();
			itl2::backprojectOpenCLProjectionBlocks(preProcProjections, settings, output2, 100);
			timer.stop();
			cout << "Backprojection took " << timer.getSeconds() << " s." << endl;

			raw::writed(output2, "./rec/fibre_reconstruction_mid_gpu_projection_blocks");



			cout << "OpenCL backprojection, output and projections in blocks" << endl;

			timer.start();
			Image<float32_t> output3;
			itl2::backprojectOpenCL(preProcProjections, settings, output3, 100, 50);
			timer.stop();
			cout << "Backprojection took " << timer.getSeconds() << " s." << endl;

			raw::writed(output3, "./rec/fibre_reconstruction_mid_gpu_output_projection_blocks");


			testAssert(equals(output, output2, 1e-4f), "images reconstructed in whole and in projection blocks do not match.");
			testAssert(equals(output, output3, 1e-4f), "images reconstructed in whole and in output and projection blocks do not match.");
#endif
		}

		void openCLBackProjectionRealBin2()
		{

#if defined(USE_OPENCL)
			string logFile = readText("C:\\mytemp\\cfrp\\CRPE_small.log");
			RecSettings settings = fromString<RecSettings>(logFile);

			cout << settings.filterType << endl;

			//Image<float32_t> projections(988, 988, 1441);
			//raw::read(projections, "C:\\mytemp\\cfrp\\CRPE_small_988x988x1441.raw");
			//cout << "Copy input to output..." << endl;
			//copyFile("C:\\mytemp\\cfrp\\CRPE_small_988x988x1441.raw", "C:\\mytemp\\cfrp\\CRPE_small_988x988x1441-preprocessed.raw");
			Image<float32_t> projections("C:\\mytemp\\cfrp\\CRPE_small", true, 988, 988, 1441);
			Image<float32_t> preProcessedProjections("C:\\mytemp\\cfrp\\CRPE_small_preprocessed", false, 988, 988, 1441);

			itl2::internals::sanityCheck(projections, settings, false);

			Timer timer;

			timer.start();
			itl2::fbpPreprocess(projections, preProcessedProjections, settings);
			timer.stop();
			cout << "Preprocessing took " << timer.getSeconds() << " s" << endl;

			timer.start();
			Image<float32_t> output;
			itl2::backprojectOpenCL(preProcessedProjections, settings, output, 8);
			timer.stop();
			cout << "Backprojection took " << timer.getSeconds() << " s." << endl;

			raw::writed(output, "C:\\mytemp\\cfrp\\CRPE_small_rec");
#endif
		}


		void openCLBackProjectionRealBin1()
		{
#if defined(USE_OPENCL)
			string logFile = readText("C:\\mytemp\\cfrp\\crpe_log_big.txt");
			RecSettings settings = fromString<RecSettings>(logFile);

			//cout << "Copy input to output..." << endl;
			//copyFile("C:\\mytemp\\cfrp\\CRPE_test_scans_20X_Tomo_AreaB_1976x1976x1441.raw", "C:\\mytemp\\cfrp\\CRPE_test_scans_20X_Tomo_AreaB_1976x1976x1441-preprocessed.raw");

			// Memory map
			Image<float32_t> preProcessedProjections("C:\\mytemp\\cfrp\\CRPE_test_scans_20X_Tomo_AreaB_preprocessed", false, 1976, 1976, 1441);
			
			Timer timer;
			{
				Image<float32_t> projections("C:\\mytemp\\cfrp\\CRPE_test_scans_20X_Tomo_AreaB", true, 1976, 1976, 1441);
				itl2::internals::sanityCheck(projections, settings, false);

				timer.start();
				itl2::fbpPreprocess(projections, preProcessedProjections, settings);
				timer.stop();
				cout << "Preprocessing took " << timer.getSeconds() << " s" << endl;
			}

			timer.start();
			Image<float32_t> output("C:\\mytemp\\cfrp\\CRPE_big_rec", false, 1976, 1976, 1976);
			//Image<float32_t> output;
			itl2::backprojectOpenCL(preProcessedProjections, settings, output, 8);
			timer.stop();
			cout << "Backprojection took " << timer.getSeconds() << " s." << endl;

			//raw::writed(output, "C:\\mytemp\\cfrp\\CRPE_big_rec");
#endif
		}



	}
}
