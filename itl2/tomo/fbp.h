#pragma once

//#if defined(_WIN32)
#if !defined(NO_OPENCL)
	#define USE_OPENCL
#else
	#pragma message("NO_OPENCL defined, skipping OpenCL compilation.")
#endif

#include <vector>

#include "image.h"
#include "iteration.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "stringutils.h"
#include "interpolation.h"
#include "math/vectoroperations.h"
#include "imagemetadata.h"
#include "tomo/siddonprojections.h"
#include "neighbourhoodtype.h"

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
		/**
		Input data is in I/I0 format and absorption reconstruction should be made.
		*/
		Absorption = 0,
		/**
		Input data is in I/I0 format and phase reconstruction should me made using Paganin method.
		*/
		Paganin,
		/**
		Input data is in -ln(I/I0) format and no negative logarithm or phase extraction is made before reconstruction.
		*/
		Direct
	};

	inline std::ostream& operator<<(std::ostream& stream, const PhaseMode& x)
	{
		switch (x)
		{
		case PhaseMode::Absorption: stream << "Absorption"; return stream;
		case PhaseMode::Paganin: stream << "Paganin"; return stream;
		case PhaseMode::Direct: stream << "Direct"; return stream;
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
		if (str == "direct")
			return PhaseMode::Direct;

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
		Filtering radius (defect radius) parameter for dead pixel removal.
		*/
		coord_t deadPixelMedianRadius = 1;

		/**
		Standard deviation (defect intensity magnitude) parameter for dead pixel removal.
		*/
		float32_t deadPixelStdDevCount = 20.0f;


		/**
		Shift between shadow of rotation axis at the camera and centerline of the camera.
		Center shift is assume to be shift in camera location (whereas xy-shifts are assumed to be shifts in sample location).
		By convention, this is shift of image, not shift of camera.
		Specify value in pixels.
		*/
		float32_t centerShift = 0.0f;

		/**
		Shift between optical axis and camera centerline in z-direction (up-down).
		By convention, this is shift of image, not shift of camera.
		Specify value in pixels.
		*/
		float32_t cameraZShift = 0.0f;

		/**
		Specify rotation angle of camera around its normal, in degrees.
		Camera rotation is not accounted for in the filtering step of FDK algorithm, but it is accounted for in
		the backprojection phase.
		*/
		float32_t cameraRotation = 0.0f;

		/**
		Rate of change of center shift as a function of angle from the image taken at central angle.
		*/
		float32_t csAngleSlope = 0.0f;

		/**
		This value is used to increase or decrease the angular range by the given amount of degrees.
		*/
		float32_t angleTweak = 0.0f;






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
		Type of reconstruction to make.
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
		Shift of sample applied after rotation.
		This shift is typically caused by stages installed on top of rotation axis.
		When rotation angle is 0 deg,
		X corresponds to the optical axis direction (from the source to the camera)
		Y corresponds to camera right.
		Z corresponds to camera up.
		*/
		std::vector<Vec3f> sampleShifts = std::vector<Vec3f>();

		/**
		Shifts of the rotation axis.
		This shift is typically caused by stages installed below the rotation axis.
		X corresponds to the direction from the source to the camera
		Y corresponds to the camera right.
		Z corresponds to the camera up.
		*/
		std::vector<Vec3f> rotationAxisShifts = std::vector<Vec3f>();

		/**
		Shifts applied to the camera.
		X corresponds to the direction from the source to the camera.
		Y corresponds to the camera right.
		Z corresponds to the camera up.
		*/
		std::vector<Vec3f> cameraShifts = std::vector<Vec3f>();

		/**
		Shifts applied to the source.
		X corresponds to the direction from the source to the camera.
		Y corresponds to the camera right.
		Z corresponds to the camera up.
		*/
		std::vector<Vec3f> sourceShifts = std::vector<Vec3f>();

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
		s.deadPixelMedianRadius = id.get("dead_pixel_median_radius", s.deadPixelMedianRadius);
		s.deadPixelStdDevCount = id.get("dead_pixel_std_dev_count", s.deadPixelStdDevCount);
		s.centerShift = id.get("center_shift", s.centerShift);
		s.cameraZShift = id.get("camera_z_shift", s.cameraZShift);
		s.cameraRotation = id.get("camera_rotation", s.cameraRotation);
		s.csAngleSlope = id.get("cs_angle_slope", s.csAngleSlope);
		s.angleTweak = id.get("angle_tweak", s.angleTweak);
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
		std::vector<Vec3f> emptyV3;
		s.sampleShifts = id.getList<Vec3f>("sample_shifts", emptyV3);
		s.sourceShifts = id.getList<Vec3f>("source_shifts", emptyV3);
		s.cameraShifts = id.getList<Vec3f>("camera_shifts", emptyV3);
		s.rotationAxisShifts = id.getList<Vec3f>("rotation_axis_shifts", emptyV3);

		// If there are no shifts supplied, set all shifts to zero.
		if (s.sampleShifts.size() <= 0)
		{
			while (s.sampleShifts.size() < s.angles.size())
				s.sampleShifts.push_back(Vec3f(0, 0, 0));
		}

		return s;
	}
	



	/**
	Pre-processing of data for filtered backprojection.
	Call this before calling backproject method.
	*/
	void fbpPreprocess(const Image<float32_t>& transmissionProjections, Image<float32_t>& preprocessedProjections, RecSettings settings);


	/**
	Remove bad pixels from one slice of projection data.
	@param slice View of the slice.
	@param med, tmp Temporary images.
	@return Count of bad pixels in the slice.
	*/
	template<typename pixel_t, typename intermediate_t = typename math_intermediate_type<pixel_t, pixel_t>::type > size_t deadPixelRemovalSlice(Image<pixel_t>& slice, Image<pixel_t>& med, Image<intermediate_t>& tmp, coord_t medianRadius = 2, float32_t stdDevCount = 30)
	{
		// Calculate median filtering of the slice
		nanMedianFilter(slice, med, medianRadius, NeighbourhoodType::Rectangular, BoundaryCondition::Nearest);

		// Calculate abs(slice - median)
		setValue(tmp, slice);
		subtract(tmp, med);
		abs(tmp);

		// Calculate its mean and standard deviation
		Vec2d v = meanAndStdDev(tmp);
		float32_t meandifference = (float32_t)v.x;
		float32_t stddifference = (float32_t)v.y;


		// Perform filtering
		size_t badPixelCount = 0;
		for (coord_t y = 0; y < slice.height(); y++)
		{
			for (coord_t x = 0; x < slice.width(); x++)
			{
				pixel_t p = slice(x, y);
				pixel_t m = med(x, y);

				if (NumberUtils<pixel_t>::isnan(p) || abs((float32_t)m - (float32_t)p) > stdDevCount * stddifference)
				{
					slice(x, y) = m;
					badPixelCount++;
				}
			}
		}

		return badPixelCount;
	}

	/**
	Removes dead pixels from each slice of the input image.
	*/
	template<typename pixel_t> void deadPixelRemoval(Image<pixel_t>& img, coord_t medianRadius, float32_t stdDevCount)
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

			Image<pixel_t> med(img.width(), img.height());
			using intermediate_t = typename math_intermediate_type<pixel_t, pixel_t>::type;
			Image<intermediate_t> tmp(img.width(), img.height());

#pragma omp for
			for (coord_t n = 0; n < img.depth(); n++)
			{
				Image<pixel_t> slice(img, n, n);

				size_t badPixelCount = deadPixelRemovalSlice(img, med, tmp, medianRadius, stdDevCount);

#pragma omp critical(badpixels)
				{
					averageBadPixels += badPixelCount;
					maxBadPixels = std::max(maxBadPixels, badPixelCount);
				}

				showThreadProgress(counter, img.depth());
			}
		}

		std::cout << "Average number of bad pixels per slice: " << averageBadPixels / (float)img.depth() << std::endl;
		std::cout << "Maximum number of bad pixels in slice: " << maxBadPixels << std::endl;

		if (maxBadPixels > 100)
			std::cout << "WARNING: Maximum number of bad pixels is high: " << maxBadPixels << ". Consider changing settings for bad pixel removal." << std::endl;
	}


	namespace internals
	{
		void sanityCheck(const Image<float32_t>& transmissionProjections, RecSettings& settings, bool projectionsAreBinned);

		float32_t calculateTrueCentralAngle(float32_t centralAngleFor180degScan, const std::vector<float32_t>& angles, float32_t gammamax0);

		float32_t calculateGammaMax0(float32_t projectionWidth, float32_t d);

		float32_t csAnglePerturbation(size_t angleIndex, float32_t centralAngle, const std::vector<float32_t>& angles, float32_t csAngleSlope);

		/**
		Calculates perturbation to be applied to CS values, given projection height and z-coordinate in the projection in range [0, projectionHeight - 1[.
		*/
		inline float32_t csZPerturbation(float32_t projectionZ, float32_t projectionHeight, float32_t csZSlope)
		{
			return (projectionZ - projectionHeight / 2) * csZSlope;
		}

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
			//normFactor = 1.0f / N * PIf;
			normFactor = 1.0f / N * PIf / 2.0f; // NOTE: This seems to be the correct factor according to tests where 180 deg and 360 deg scans are compared to each other
		}
		else
		{
			// 360 deg scan
			float N = (float)settings.angles.size();
			normFactor = 1.0f / N * PIf / 2.0f;
		}

		return normFactor;
	}

	namespace internals
	{
		inline void fillTo(std::vector<Vec3f>& v, size_t count)
		{
			while (v.size() < count)
				v.push_back(Vec3f());
		}

		/**
		Calculates backprojection geometry (source position, detector position, detector right, detector up, detector normal)
		from reconstruction settings.
		The output vector type is templated as OpenCL processing needs the output in Vec4f, but CPU processing is
		fine with Vec3f.
		*/
		template<typename VEC> void determineBackprojectionGeometry(RecSettings settings, coord_t projectionWidth,
			std::vector<VEC>& pss,		// Source positions
			std::vector<VEC>& pds,		// Detector center point positions
			std::vector<VEC>& us,			// Detector right vectors
			std::vector<VEC>& vs,			// Detector up vectors
			std::vector<VEC>& ws			// Detector normal vectors, w = u x w
			)
		{
			pss.clear();
			pds.clear();
			us.clear();
			vs.clear();
			ws.clear();

			size_t projCount = settings.angles.size();
			pss.reserve(projCount);
			pds.reserve(projCount);
			us.reserve(projCount);
			vs.reserve(projCount);
			ws.reserve(projCount);

			// Ensure that we have all shifts available.
			fillTo(settings.cameraShifts, projCount);
			fillTo(settings.sourceShifts, projCount);
			fillTo(settings.rotationAxisShifts, projCount);
			fillTo(settings.sampleShifts, projCount);


			float32_t gammamax0 = internals::calculateGammaMax0((float32_t)projectionWidth, settings.sourceToRA);
			float32_t centralAngle = internals::calculateTrueCentralAngle(settings.centralAngleFor180degScan, settings.angles, gammamax0);

			float32_t rotMul = 1 + (settings.angleTweak / max(settings.angles));
			if (settings.rotationDirection == RotationDirection::Counterclockwise)
				rotMul *= -1;

			if (!settings.useShifts)
				settings.shiftScaling = 0;

			float32_t M = 1 + settings.objectCameraDistance / settings.sourceToRA;
			for (size_t anglei = 0; anglei < settings.angles.size(); anglei++)
			{
				//double angleRad = rotMul * ((double)settings.angles[anglei] - 90 + (double)settings.rotation) / 180.0 * PI;

				//float32_t c = (float32_t)cos(angle);
				//float32_t s = (float32_t)sin(angle);

				//Vec3f u(-s, c, 0); // At 0 deg the camera right vector is +y
				//Vec3f v(0, 0, 1);
				//Vec3f w = u.cross(v);
				//Vec3f ps = Vec3f(0, 0, 0) - settings.sourceToRA * w; // TODO: Add source shift
				//Vec3f pd = Vec3f(0, 0, 0) + settings.objectCameraDistance * w;


				Vec3f u(0, 1, 0);	// Camera right at 0 deg is +y
				Vec3f v(0, 0, 1);	// Camera up at 0 deg is +z
				Vec3f ps(-settings.sourceToRA, 0, 0);
				Vec3f pd(settings.objectCameraDistance, 0, 0);

				ps -= settings.rotationAxisShifts[anglei] * settings.shiftScaling;
				pd -= settings.rotationAxisShifts[anglei] * settings.shiftScaling;

				ps += settings.sourceShifts[anglei] * settings.shiftScaling;
				pd += settings.cameraShifts[anglei] * settings.shiftScaling;
				

				float32_t angleRad = (rotMul * settings.angles[anglei] + settings.rotation) / 180.0f * PIf;
				itl2::rotate<float32_t>(ps, angleRad, 0.0f);
				itl2::rotate<float32_t>(pd, angleRad, 0.0f);
				itl2::rotate<float32_t>(u, angleRad, 0.0f);
				itl2::rotate<float32_t>(v, angleRad, 0.0f);
				Vec3f w = u.cross(v);

				//float32_t sdx = settings.objectShifts[anglei].x * settings.shiftScaling * (settings.useShifts ? 1 : 0);
				//float32_t sdz = settings.objectShifts[anglei].y * settings.shiftScaling * (settings.useShifts ? 1 : 0);
				//Vec3f objShift = sdx * u + sdz * v;
				//ps -= objShift;
				//pd -= objShift;

				// Add sample fine alignment shifts by shifting both source and camera to the inverse direction
				ps -= settings.sampleShifts[anglei] * settings.shiftScaling;
				pd -= settings.sampleShifts[anglei] * settings.shiftScaling;


				// Add camera calibration shifts
				// Camera shift in Z-direction and center shift (in optical plane) are given in image pixels, so they must
				// be multiplied by magnification to get the correct effect on camera position.
				float32_t csAnglePert = internals::csAnglePerturbation(anglei, centralAngle, settings.angles, settings.csAngleSlope);
				pd += Vec3f(0, 0, -settings.cameraZShift * M) + u * (-settings.centerShift + -csAnglePert) * M;

				// Add camera rotation around its normal
				float32_t cameraRotRad = settings.cameraRotation / 180.0f * PIf;
				u = u.rotate(w, cameraRotRad);
				v = v.rotate(w, cameraRotRad);


				pss.push_back(VEC(ps));
				pds.push_back(VEC(pd));
				us.push_back(VEC(u / M)); // NOTE: We store u / M and v / M instead of plain u and v as this way we save two vector divisions in the inner backprojection loop.
				vs.push_back(VEC(v / M));
				ws.push_back(VEC(w));
			}
		}
	}

	template<typename out_t> void backproject(const Image<float32_t>& transmissionProjections, RecSettings settings, Image<out_t>& output)
	{
		internals::sanityCheck(transmissionProjections, settings, true);
		output.mustNotBe(transmissionProjections);

		internals::applyBinningToParameters(settings);

		output.ensureSize(settings.roiSize);

		float32_t normFact = normFactor(settings);

		// Pre-calculate direction and position vectors that define the backprojection geometry
		// ------------------------------------------------------------------------------------
		std::vector<Vec3f> pss;			// Source positions
		std::vector<Vec3f> pds;			// Detector center point positions
		std::vector<Vec3f> us;			// Detector right vectors
		std::vector<Vec3f> vs;			// Detector up vectors
		std::vector<Vec3f> ws;			// Detector normal vectors, w = u x w
		internals::determineBackprojectionGeometry(settings, transmissionProjections.width(), pss, pds, us, vs, ws);

		// Backproject
		// -----------
		float32_t projectionHalfWidth = (float32_t)transmissionProjections.width() / 2.0f;
		float32_t projectionHalfHeight = (float32_t)transmissionProjections.height() / 2.0f;

		LinearInterpolator<float32_t, float32_t> interpolator(BoundaryCondition::Zero);

		size_t counter = 0;
		forAllPixels(output, [&](coord_t x, coord_t y, coord_t z)
			{
				// NOTE: Adding 0.5 ensures that if x = 0, roiSize = 1, and roiCenter = 0,
				// p.x = 0 - 1/2.0f + c + 0.5 = 0, as expected.
				// NOTE: In this new version we use the more correct +0.5 in all coordinate directions. This is different from the old version
				// where +0.5 was used only in the z direction
				Vec3f p = Vec3f((float)x, (float)y, (float)z) - Vec3f(settings.roiSize) / 2.0f + Vec3f(settings.roiCenter)
					//+ Vec3f(0, 0, 0.5); // Old version was like this
					+ Vec3f(0.5, 0.5, 0.5);

				// Sum contributions from all projections
				float32_t sum = 0;
				for (coord_t anglei = 0; anglei < transmissionProjections.depth(); anglei++)
				{
					Vec3f ps = pss[anglei];
					Vec3f pd = pds[anglei];
					Vec3f uVec = us[anglei];
					Vec3f vVec = vs[anglei];
					Vec3f wHat = ws[anglei];

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

						// This weight is needed in the FDK algorithm.
						float32_t weight = settings.sourceToRA / (settings.sourceToRA + p.dot(wHat));
						weight *= weight;

						sum += weight * interpolator(transmissionProjections, u, v, (float32_t)anglei);
					}
				}

				sum *= normFact;

				// Scaling
				sum = (sum - settings.dynMin) / (settings.dynMax - settings.dynMin) * NumberUtils<out_t>::scale();
				output(x, y, z) = pixelRound<out_t>(sum);
			},
			true);
	}




#if defined(USE_OPENCL)
	void backprojectOpenCL(const Image<float32_t>& transmissionProjections, RecSettings settings, Image<float32_t>& output, coord_t maxOutputBlockSizeZ = 8, coord_t maxProjectionBlockSizeZ = std::numeric_limits<coord_t>::max());
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
