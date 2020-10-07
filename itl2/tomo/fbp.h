#pragma once

//#if defined(_WIN32)
#if !defined(NO_OPENCL)
	#define USE_OPENCL
#else
	#pragma message("NO_OPENCL defined, skipping OpenCL compilation.")
#endif

#include <vector>

#include "image.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "stringutils.h"
#include "interpolation.h"
#include "math/vectoroperations.h"
#include "imagemetadata.h"

namespace itl2
{
	enum class RotationDirection
	{
		Clockwise = 0,
		Counterclockwise
	};

	inline std::ostream& operator<<(std::ostream& stream, const RotationDirection& x)
	{
		switch (x)
		{
		case RotationDirection::Clockwise: stream << "Clockwise"; return stream;
		case RotationDirection::Counterclockwise: stream << "Counterclockwise"; return stream;
		}
		throw ITLException("Invalid rotation direction.");
	}

	template<>
	inline RotationDirection fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "clockwise")
			return RotationDirection::Clockwise;
		if (str == "counterclockwise")
			return RotationDirection::Counterclockwise;

		throw ITLException("Invalid rotation direction: " + str0);
	}



	enum class PadType
	{
		Mirror = 0,
		Nearest
	};

	inline std::ostream& operator<<(std::ostream& stream, const PadType& x)
	{
		switch (x)
		{
		case PadType::Mirror: stream << "Mirror"; return stream;
		case PadType::Nearest: stream << "Nearest"; return stream;
		}
		throw ITLException("Invalid pad type.");
	}

	template<>
	inline PadType fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "mirror")
			return PadType::Mirror;
		if(str == "nearest")
			return PadType::Nearest;

		throw ITLException("Invalid pad type: " + str0);
	}



	
	enum class PhaseMode
	{
		Absorption = 0,
		Paganin
	};

	inline std::ostream& operator<<(std::ostream& stream, const PhaseMode& x)
	{
		switch (x)
		{
		case PhaseMode::Absorption: stream << "Absorption"; return stream;
		case PhaseMode::Paganin: stream << "Paganin"; return stream;
		}
		throw ITLException("Invalid phase mode.");
	}

	template<>
	inline PhaseMode fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "absorption")
			return PhaseMode::Absorption;
		if (str == "paganin")
			return PhaseMode::Paganin;

		throw ITLException("Invalid phase mode: " + str0);
	}




	enum class FilterType
	{
		/**
		Standard |w| filter.
		*/
		IdealRamp = 0,
		/**
		Standard |w| filter adjusted to keep image average as it should be.
		*/
		Ramp,
		/**
		Shepp-Logan filter.
		*/
		SheppLogan,
		/**
		Cosine filter.
		*/
		Cosine,
		/**
		Hamming window.
		*/
		Hamming,
		/**
		Hann filter.
		*/
		Hann,
		/**
		Blackman filter.
		*/
		Blackman,
		/**
		Parzen filter.
		*/
		Parzen
	};

	inline std::ostream& operator<<(std::ostream& stream, const FilterType& x)
	{
		switch (x)
		{
		case FilterType::IdealRamp: stream << "Ideal ramp"; return stream;
		case FilterType::Ramp: stream << "Ramp"; return stream;
		case FilterType::SheppLogan: stream << "Shepp-Logan"; return stream;
		case FilterType::Cosine: stream << "Cosine"; return stream;
		case FilterType::Hamming: stream << "Hamming"; return stream;
		case FilterType::Hann: stream << "Hann"; return stream;
		case FilterType::Blackman: stream << "Blackman"; return stream;
		case FilterType::Parzen: stream << "Parzen"; return stream;
		}
		throw ITLException("Invalid filter type.");
	}

	template<>
	inline FilterType fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "ideal ramp")
			return FilterType::IdealRamp;
		if (str == "ramp")
			return FilterType::Ramp;
		if (str == "shepp-logan")
			return FilterType::SheppLogan;
		if (str == "cosine")
			return FilterType::Cosine;
		if (str == "hamming")
			return FilterType::Hamming;
		if (str == "hann")
			return FilterType::Hann;
		if (str == "blackman")
			return FilterType::Blackman;
		if (str == "parzen")
			return FilterType::Parzen;

		throw ITLException("Invalid filter type: " + str0);
	}

	/**
	Constructs a backprojection filter into given image. Image width must be set to desired value.
	*/
	void createFilter(Image<float32_t>& filter, FilterType filterType, float32_t cutoff);


	struct RecSettings
	{

		/**
		Distance between rotation axis and X-ray source, or more generally distance between radiation cone tip and rotation axis.
		Given in pixels.
		Set to zero to assume parallel beam geometry.
		*/
		float32_t sourceToRA = 0;

		/**
		Rotation direction, clockwise or counterclockwise.
		*/
		RotationDirection rotationDirection = RotationDirection::Clockwise;



		/**
		Beam hardening correction parameter.
		This is multiplier of the second-degree term in correction equation
		x' = x + bhc * x^2,
		where x = -log(I/I_0) and
		x' is the corrected value.
		*/
		float32_t bhc = 0.0f;

		/**
		Flag that indicates whether the data represents 180 deg scan or not.
		*/
		bool reconstructAs180degScan = true;

		/**
		If the scan is 180 degree scan, 180 deg angular range around this angle is used for the reconstruction.
		*/
		float32_t centralAngleFor180degScan = 0.0f;

		/**
		Set to 1 to disable this heuristic parameter. Set to positive value (e.g. 6) to smooth edges of sinogram more further from optical axis.
		*/
		float32_t heuristicSinogramWindowingParameter = 1.0f;

		/**
		Use to rotate the output image around rotation axis without interpolation.
		Give value in degrees.
		*/
		float32_t rotation = 0.0f;

		/**
		Reconstruction ROI start and size.
		Zero position is at the intersection of optical and rotation axes.
		Specify zeroes to ROI size to reconstruct full image.
		*/
		Vec3c roiCenter = Vec3c(0, 0, 0);
		Vec3c roiSize = Vec3c(0, 0, 0);

		/**
		Crop this many pixels away from edges of projections in horizontal and vertical direction.
		Cropping is done before binning.
		*/
		Vec2c cropSize = Vec2c(0, 0);

		/**
		Use to bin projection data before reconstruction (for faster reconstruction, better signal-to-noise ratio and lower resolution).
		*/
		size_t binning = 1;

		/**
		Should dead/stuck pixels be corrected?
		*/
		bool removeDeadPixels = false;


		/**
		Shift between shadow of rotation axis at camera and centerline of the camera.
		Center shift is assume to be shift in camera location (whereas xy-shifts are assumed to be shifts in sample location).
		Specify value in pixels.
		*/
		float32_t centerShift = 0.0f;

		/**
		Shift between optical axis and camera centerline in z-direction (up-down).
		Specify value in pixels.
		*/
		float32_t cameraZShift = 0.0f;

		/**
		Rate of change of center shift as a function of angle from the image taken at central angle.
		*/
		float32_t csAngleSlope = 0.0f;

		/**
		Rate of change of center shift as a function of z-directional distance from optical axis.
		*/
		float32_t csZSlope = 0.0f;






		/**
		Type of padding used in backprojection filtering.
		*/
		PadType padType = PadType::Mirror;

		/**
		Size of padding used in backprojection filtering. Specify as a fraction of image size.
		*/
		float32_t padFraction = 0.5f;

		/**
		Type of filter to use.
		*/
		FilterType filterType = FilterType::Ramp;

		/**
		Cut-off frequency, 0 corresponds to DC and 1 to Nyquist.
		*/
		float32_t filterCutOff = 1.0f;





		/**
		Type of phase retrieval used.
		*/
		PhaseMode phaseMode = PhaseMode::Absorption;

		/**
		Size of padding used in phase retrieval. Specify as a fraction of image size.
		*/
		float32_t phasePadFraction = 0.5f;

		/**
		Type of padding used in backprojection filtering.
		*/
		PadType phasePadType = PadType::Mirror;

		/**
		Distance between object and camera.
		Specify value in pixels.
		*/
		float32_t objectCameraDistance = 1000.0f;

		/**
		"delta" parameter for Paganin phase retrieval.
		Refractive index is 1-delta + beta * i, and
		attenuation coefficient mu = 4 * pi * beta / lambda,
		where lambda is wavelength.
		*/
		float32_t delta = 5e-16f;

		/**
		Effective attenuation coefficient for Paganin phase retrieval.
		Refractive index is 1-delta + beta * i, and
		attenuation coefficient mu = 4 * pi * beta / lambda,
		where lambda is wavelength.
		*/
		float32_t mu = 1e-7f;


		


		/**
		Values that are mapped to 0 and maximum possible value in output pixel data type.
		If output type is float32, maximum is assumed to be 1, and no clipping is performed.
		Therefore, when using float32 data type, set dynMin to 0 and dynMax to 1 to retain
		original values of the reconstructed data.
		*/
		float32_t dynMin = 0.0f, dynMax = 1.0f;

		/**
		Rotation angle for each projection image in degrees.
		*/
		std::vector<float32_t> angles = std::vector<float32_t>();

		/**
		Shift of object from its initial location for each projection image.
		*/
		std::vector<Vec2f> objectShifts = std::vector<Vec2f>();

		/**
		Scaling factor that will be applied to the object shift values during reconstruction.
		Use to, e.g. correct for inaccurately calibrated pixel size that was used for determination of object shifts.
		*/
		float32_t shiftScaling = 1.0f;

		/**
		Indicates if XY shifts should be accounted for.
		*/
		bool useShifts = true;


		/**
		Convert reconstruction settings to string.
		*/
		friend std::ostream& operator<<(std::ostream& stream, const RecSettings& s);
	};

    

	/**
	Create reconstruction settings from string similar to what << operator outputs.
	*/
	template<>
	inline RecSettings fromString(const string& strOrig)
	{
		// This makes default settings
		RecSettings s;

		ImageMetadata id;
		id.readFromString(strOrig);

		s.sourceToRA = id.get("source_to_ra", s.sourceToRA);
		s.rotationDirection = id.get("rotation_direction", s.rotationDirection);
		s.bhc = id.get("bhc", s.bhc);
		s.reconstructAs180degScan = id.get("rec_as_180_deg_scan", s.reconstructAs180degScan);
		s.centralAngleFor180degScan = id.get("central_angle", s.centralAngleFor180degScan);
		s.heuristicSinogramWindowingParameter = id.get("hswp", s.heuristicSinogramWindowingParameter);
		s.rotation = id.get("rotation", s.rotation);
		s.roiCenter = id.get("roi_center", s.roiCenter);
		s.roiSize = id.get("roi_size", s.roiSize);
		s.cropSize = id.get("crop_size", s.cropSize);
		s.binning = id.get("binning", s.binning);
		s.removeDeadPixels = id.get("remove_dead_pixels", s.removeDeadPixels);
		s.centerShift = id.get("center_shift", s.centerShift);
		s.cameraZShift = id.get("camera_z_shift", s.cameraZShift);
		s.csAngleSlope = id.get("cs_angle_slope", s.csAngleSlope);
		s.csZSlope = id.get("cs_z_slope", s.csZSlope);
		s.padType = id.get("pad_type", s.padType);
		s.padFraction = id.get("pad_size", s.padFraction);
		s.filterType = id.get("filter_type", s.filterType);
		s.filterCutOff = id.get("filter_cut_off", s.filterCutOff);
		s.phaseMode = id.get("phase_mode", s.phaseMode);
		s.phasePadType = id.get("phase_pad_type", s.phasePadType);
		s.phasePadFraction = id.get("phase_pad_size", s.phasePadFraction);
		s.objectCameraDistance = id.get("propagation_distance", s.objectCameraDistance);
		s.delta = id.get("delta", s.delta);
		s.mu = id.get("mu", s.mu);
		s.dynMin = id.get("range_min", s.dynMin);
		s.dynMax = id.get("range_max", s.dynMax);
		s.shiftScaling = id.get("shift_scale", s.shiftScaling);
		s.useShifts = id.get("use_shifts", s.useShifts);

        std::vector<float32_t> emptyV1;
		s.angles = id.getList<float32_t>("angles", emptyV1);
        std::vector<Vec2f> emptyV2;
		s.objectShifts = id.getList<Vec2f>("sample_shifts", emptyV2);

		// If there are no shifts supplied, set all shifts to zero.
		if (s.objectShifts.size() <= 0)
		{
			while (s.objectShifts.size() < s.angles.size())
				s.objectShifts.push_back(Vec2f(0, 0));
		}

		return s;
	}
	



	/**
	Pre-processing of data for filtered backprojection.
	Call this before calling backproject method.
	*/
	void fbpPreprocess(const Image<float32_t>& transmissionProjections, Image<float32_t>& preprocessedProjections, RecSettings settings);

	namespace internals
	{
		void sanityCheck(const Image<float32_t>& transmissionProjections, RecSettings& settings, bool projectionsAreBinned);

		float32_t calculateTrueCentralAngle(float32_t centralAngleFor180degScan, const std::vector<float32_t>& angles, float32_t gammamax0);

		float32_t calculateGammaMax0(float32_t projectionWidth, float32_t d);

		float32_t csAnglePerturbation(size_t angleIndex, float32_t centralAngle, const std::vector<float32_t>& angles, float32_t csAngleSlope);

		float32_t  csZPerturbation(float32_t z, float32_t roiCenterZ, float32_t roiSizeZ, float32_t projectionHeight, float32_t csZSlope);

		/**
		Adjusts reconstruction settings for binning.
		*/
		void applyBinningToParameters(RecSettings& settings);
	}

	/**
	Calculates normalization factor for FBP. Output pixel values must be multiplied by the returned value.
	*/
	inline float normFactor(const RecSettings& settings)
	{
		float32_t normFactor;
		if (settings.reconstructAs180degScan)
		{
			// 180 deg scan, calculate fraction of angles that form the 180 deg rotation.
			float angleFraction = 180.0f / (itl2::max(settings.angles) - itl2::min(settings.angles));
			float N = settings.angles.size() * angleFraction;
			normFactor = 1.0f / N * PIf;
		}
		else
		{
			// 360 deg scan
			float N = (float)settings.angles.size();
			normFactor = 1.0f / N * PIf / 2.0f;
		}

		return normFactor;
	}

	/**
	Backprojection step of filtered backprojection algorithm.
	*/
	template<typename out_t> void backproject(const Image<float32_t>& transmissionProjections, RecSettings settings, Image<out_t>& output)
	{
		internals::sanityCheck(transmissionProjections, settings, true);
		output.mustNotBe(transmissionProjections);
		
		internals::applyBinningToParameters(settings);

		output.ensureSize(settings.roiSize);

		float32_t gammamax0 = internals::calculateGammaMax0((float32_t)transmissionProjections.width(), settings.sourceToRA);
		float32_t centralAngle = internals::calculateTrueCentralAngle(settings.centralAngleFor180degScan, settings.angles, gammamax0);

		float32_t normFact = normFactor(settings);
		


		// Pre-calculate direction vectors.
		std::vector<Vec3f> xHatArray;
		std::vector<Vec3f> nHatArray;
		xHatArray.reserve(settings.angles.size());
		nHatArray.reserve(settings.angles.size());

		double rotMul = 1.0;
		if (settings.rotationDirection == RotationDirection::Counterclockwise)
			rotMul = -1.0;

		for (size_t anglei = 0; anglei < settings.angles.size(); anglei++)
		{
			double angle = rotMul * ((double)settings.angles[anglei] - 90 + (double)settings.rotation) / 180.0 * PI;

			float32_t c = (float32_t)cos(angle);
			float32_t s = (float32_t)sin(angle);
			xHatArray.push_back(Vec3f(c, s, 0));
			nHatArray.push_back(Vec3f(-s, c, 0));
		}
		Vec3f zHat(0, 0, 1);

		// Pre-calculate center shift angle perturbations
		std::vector<float32_t> csAnglePerturbations;
		csAnglePerturbations.reserve(settings.angles.size());
		for (size_t anglei = 0; anglei < settings.angles.size(); anglei++)
		{
			csAnglePerturbations.push_back(internals::csAnglePerturbation(anglei, centralAngle, settings.angles, settings.csAngleSlope));
		}

		// Backproject
		float32_t projectionWidth = (float32_t)transmissionProjections.width();
		float32_t projectionHeight = (float32_t)transmissionProjections.height();
		float32_t d = settings.sourceToRA;

		LinearInterpolator<float32_t, float32_t> interpolator(BoundaryCondition::Zero);

		//output.ensureSize(settings.roiSize);

		// NOTE: -0.5 ensures that if roiSize.z == 1, roiCenter.z = 1 / 2 - roiCenter.z - 0.5 = roiCenter.z
		// TODO: Should subtraction be done for all the components?
		Vec3f center = Vec3f(settings.roiSize) / 2.0f - Vec3f(settings.roiCenter) - Vec3f(0, 0, 0.5);

		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel() && settings.roiSize.z > 1)
		for(coord_t z = 0; z < settings.roiSize.z; z++)
		{
			#pragma omp parallel for if(!omp_in_parallel() && settings.roiSize.y > 1)
			for(coord_t y = 0; y < settings.roiSize.y; y++)
			{
				#pragma omp parallel for if(!omp_in_parallel() && settings.roiSize.x > 1)
				for(coord_t x = 0; x < settings.roiSize.x; x++)
				{
					float32_t currentCS = settings.centerShift + internals::csZPerturbation((float32_t)z, (float32_t)settings.roiCenter.z, (float32_t)settings.roiSize.z, projectionHeight, settings.csZSlope);

					// Sum contributions from all projections
					float32_t sum = 0;
					for (coord_t anglei = 0; anglei < transmissionProjections.depth(); anglei++)
					{
						Vec3f rho = Vec3f((float32_t)x, (float32_t)y, (float32_t)z) - center;


						Vec3f xHat = xHatArray[anglei];
						Vec3f nHat = nHatArray[anglei];
						float32_t dprhox = d + rho.dot(xHat);
						float32_t Y = (d * rho.dot(nHat)) / dprhox;
						float32_t Z = (d * rho.dot(zHat)) / dprhox;

						float32_t w = (d * d) / (dprhox * dprhox);

						// Get data from ideal detector position (Y, Z)
						// First account for detector shifts and rotation
						// TODO: Actually we have object shifts so this is only approximation that is correct for parallel beam case.
						float32_t sdx = settings.objectShifts[anglei].x * settings.shiftScaling * (settings.useShifts ? 1 : 0);
						float32_t sdz = settings.objectShifts[anglei].y * settings.shiftScaling * (settings.useShifts ? 1 : 0);
						float32_t angleCS = currentCS + csAnglePerturbations[anglei];
						float32_t ix = Y + projectionWidth / 2.0f + angleCS -sdx;
						float32_t iy = Z + projectionHeight / 2.0f + settings.cameraZShift -sdz;

						// TODO: Handle camera rotation here

						sum += w * interpolator(transmissionProjections, ix, iy, (float32_t)anglei);
					}

					sum *= normFact;

					// Scaling
					sum = (sum - settings.dynMin) / (settings.dynMax - settings.dynMin) * NumberUtils<out_t>::scale();
					//output(x, y, zi) = pixelRound<out_t>(sum);
					output(x, y, z) = pixelRound<out_t>(sum);
				}
			}

			showThreadProgress(counter, settings.roiSize.z);
		}
	}

#if defined(USE_OPENCL)
	void backprojectOpenCLProjectionOutputBlocks(const Image<float32_t>& transmissionProjections, RecSettings settings, Image<float32_t>& output, coord_t maxOutputBlockSizeZ = 8, coord_t maxProjectionBlockSizeZ = std::numeric_limits<coord_t>::max());
#endif


	namespace tests
	{
		void recSettings();
		void fbp();
		void paganin();

		void openCLBackProjection();
		void openCLBackProjectionRealBin2();
		void openCLBackProjectionRealBin1();
	}
}
