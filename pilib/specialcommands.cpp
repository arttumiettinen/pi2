
#include "specialcommands.h"
#include "pisystem.h"
#include "distributedimage.h"
#include "utilities.h"
#include "stringutils.h"
#include <omp.h>
#include "io/vol.h"
#include "io/io.h"

#include "commandmacros.h"

namespace pilib
{
	void addSpecialCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			new ClearCommand(),
			new DistributeCommand(),
			new EchoCommandsCommand(),
			new HelloCommand(),
			new HelpCommand(),
			new InfoCommand(),
			new ListCommand(),
			new LicenseCommand(),
			new ReadCommand(),
			new MapRawCommand(),
			new NewImageCommand(),
			new ReadBlockCommand(),
			new ReadRawBlockCommand(),
			new ReadRawCommand(),
			new ReadSequenceBlockCommand(),
			new ReadSequenceCommand(),
			new ReadVolCommand(),
			new WaitReturnCommand()
			}
		);
	}




	void InfoCommand::run(vector<ParamVariant>& args) const
	{
		cout << "Process Image = pi2 Copyright (C) 2018 Arttu Miettinen" << endl;
		cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
		cout << "This is free software, and you are welcome to redistribute it" << endl;
		cout << "under certain conditions; run `license()' command for more details." << endl;
		cout << endl;
		cout << "Based on work done at" << endl;
		cout << "X-ray tomography research group, TOMCAT beamline, Swiss Light Source, Paul Scherrer Institute, Switzerland" << endl;
		cout << "Centre d'Imagerie BioMedicale (CIBM), Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland" << endl;
		cout << "Complex materials research group, Department of Physics, University of Jyvaskyla, Finland" << endl;
		cout << "Contact: arttu.miettinen@psi.ch" << endl;
		cout << endl;

		if (sizeof(void*) == 4)
			cout << "32-bit version" << endl;
		else if (sizeof(void*) == 8)
			cout << "64-bit version" << endl;
		else
			cout << "Warning: Using neither 32-bit nor 64-bit version." << endl;

		cout << "Number of threads: " << omp_get_max_threads() << endl;
		cout << "Available RAM: " << bytesToString((double)memorySize()) << endl;
	}

	void LicenseCommand::run(vector<ParamVariant>& args) const
	{
		cout << "This program is licensed under GNU General Public License Version 3." << endl;
		cout << "Please see LICENSE.txt bundled with this software or" << endl;
		cout << "https://www.gnu.org/licenses/gpl-3.0.html for more information." << endl;




		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "Windows versions of this program contain code from libtiff library licensed under the libtiff license:" << endl;
		cout << R"END(
Copyright(c) 1988 - 1997 Sam Leffler
Copyright(c) 1991 - 1997 Silicon Graphics, Inc.

Permission to use, copy, modify, distribute, and sell this software and
its documentation for any purpose is hereby granted without fee, provided
that(i) the above copyright notices and this permission notice appear in
all copies of the software and related documentation, and (ii)the names of
Sam Leffler and Silicon Graphics may not be used in any advertising or
publicity relating to the software without the specific, prior written
permission of Sam Leffler and Silicon Graphics.

THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
R ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
OF THIS SOFTWARE.
)END" << endl;
		cout << endl;



		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << R"END(This program contains code owned by Joachim Kopp, used for some matrix operations and licensed
under the terms of the GNU Lesser General Public License(LGPL), see
https://www.gnu.org/licenses/lgpl.txt for details. The original source code can be
found at https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/.
The functionality of the code is further described in article
Joachim Kopp - Efficient numerical diagonalization of hermitian 3x3 matrices
Int. J. Mod. Phys. C 19 (2008) 523-548
arXiv.org: physics/0610206)END" << endl;
		cout << endl;


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from https://github.com/ITHare/util with the following license:" << endl;
		cout << R"END(
Copyright(c) 2018, ITHare.com
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :

*Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.)END" << endl;
		cout << endl;

		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from http://www.davekoelle.com/files/alphanum.hpp licensed under the MIT license:" << endl;
		cout << R"END(
Released under the MIT License - https://opensource.org/licenses/MIT

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.)END" << endl;
		cout << endl;

		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "Windows versions of this program contain code from libpng library licensed under the libpng license:" << endl;
		cout << R"END(
COPYRIGHT NOTICE, DISCLAIMER, and LICENSE:

If you modify libpng you may insert additional notices immediately following
this sentence.

This code is released under the libpng license.

libpng versions 1.0.7, July 1, 2000 through 1.6.34, September 29, 2017 are
Copyright (c) 2000-2002, 2004, 2006-2017 Glenn Randers-Pehrson, are
derived from libpng-1.0.6, and are distributed according to the same
disclaimer and license as libpng-1.0.6 with the following individuals
added to the list of Contributing Authors:

   Simon-Pierre Cadieux
   Eric S. Raymond
   Mans Rullgard
   Cosmin Truta
   Gilles Vollant
   James Yu
   Mandar Sahastrabuddhe
   Google Inc.
   Vadim Barkov

and with the following additions to the disclaimer:

   There is no warranty against interference with your enjoyment of the
   library or against infringement.  There is no warranty that our
   efforts or the library will fulfill any of your particular purposes
   or needs.  This library is provided with all faults, and the entire
   risk of satisfactory quality, performance, accuracy, and effort is with
   the user.

Some files in the "contrib" directory and some configure-generated
files that are distributed with libpng have other copyright owners and
are released under other open source licenses.

libpng versions 0.97, January 1998, through 1.0.6, March 20, 2000, are
Copyright (c) 1998-2000 Glenn Randers-Pehrson, are derived from
libpng-0.96, and are distributed according to the same disclaimer and
license as libpng-0.96, with the following individuals added to the list
of Contributing Authors:

   Tom Lane
   Glenn Randers-Pehrson
   Willem van Schaik

libpng versions 0.89, June 1996, through 0.96, May 1997, are
Copyright (c) 1996-1997 Andreas Dilger, are derived from libpng-0.88,
and are distributed according to the same disclaimer and license as
libpng-0.88, with the following individuals added to the list of
Contributing Authors:

   John Bowler
   Kevin Bracey
   Sam Bushell
   Magnus Holmgren
   Greg Roelofs
   Tom Tanner

Some files in the "scripts" directory have other copyright owners
but are released under this license.

libpng versions 0.5, May 1995, through 0.88, January 1996, are
Copyright (c) 1995-1996 Guy Eric Schalnat, Group 42, Inc.

For the purposes of this copyright and license, "Contributing Authors"
is defined as the following set of individuals:

   Andreas Dilger
   Dave Martindale
   Guy Eric Schalnat
   Paul Schmidt
   Tim Wegner

The PNG Reference Library is supplied "AS IS".  The Contributing Authors
and Group 42, Inc. disclaim all warranties, expressed or implied,
including, without limitation, the warranties of merchantability and of
fitness for any purpose.  The Contributing Authors and Group 42, Inc.
assume no liability for direct, indirect, incidental, special, exemplary,
or consequential damages, which may result from the use of the PNG
Reference Library, even if advised of the possibility of such damage.

Permission is hereby granted to use, copy, modify, and distribute this
source code, or portions hereof, for any purpose, without fee, subject
to the following restrictions:

  1. The origin of this source code must not be misrepresented.

  2. Altered versions must be plainly marked as such and must not
     be misrepresented as being the original source.

  3. This Copyright notice may not be removed or altered from any
     source or altered source distribution.

The Contributing Authors and Group 42, Inc. specifically permit, without
fee, and encourage the use of this source code as a component to
supporting the PNG file format in commercial products.  If you use this
source code in a product, acknowledgment is not required but would be
appreciated.

END OF COPYRIGHT NOTICE, DISCLAIMER, and LICENSE.

TRADEMARK:

The name "libpng" has not been registered by the Copyright owner
as a trademark in any jurisdiction.  However, because libpng has
been distributed and maintained world-wide, continually since 1995,
the Copyright owner claims "common-law trademark protection" in any
jurisdiction where common-law trademark is recognized.

OSI CERTIFICATION:

Libpng is OSI Certified Open Source Software.  OSI Certified Open Source is
a certification mark of the Open Source Initiative. OSI has not addressed
the additional disclaimers inserted at version 1.0.7.

EXPORT CONTROL:

The Copyright owner believes that the Export Control Classification
Number (ECCN) for libpng is EAR99, which means not subject to export
controls or International Traffic in Arms Regulations (ITAR) because
it is open source, publicly available software, that does not contain
any encryption software.  See the EAR, paragraphs 734.3(b)(3) and
734.7(b).

Glenn Randers-Pehrson
glennrp at users.sourceforge.net
September 29, 2017)END" << endl;
		cout << endl;


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "Windows versions of this program contain code from zlib library licensed with the following notice:" << endl;
		cout << R"END(
(C)1995 - 2017 Jean - loup Gailly and Mark Adler

This software is provided 'as-is', without any express or implied
warranty.In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions :

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software.If you use this software
in a product, an acknowledgment in the product documentation would be
appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Jean - loup Gailly        Mark Adler
jloup@gzip.org          madler@alumni.caltech.edu

If you use the zlib library in a product, we would appreciate * not* receiving
lengthy legal documents to sign.The sources are provided for free but without
warranty of any kind.The library has been entirely written by Jean - loup
Gailly and Mark Adler; it does not include third - party code.

If you redistribute modified sources, we would appreciate that you include in
the file ChangeLog history information documenting your changes.Please read
the FAQ for more information on the distribution of modified source versions.)END" << endl;
		cout << endl;


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "Windows versions of this program are shipped with FFTW3 library licensed" << endl;
		cout << "under the GNU General Public License Version 2." << endl;
		cout << "Please see https://www.gnu.org/licenses/gpl-2.0.html for full license text." << endl;


	}

	void ListCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		double totalSize = 0;
		for (const string& name : system->getImageNames())
		{
			ImageBase* img = system->getImage(name);
			auto img8 = dynamic_cast<Image<uint8_t>*>(img);
			auto img16 = dynamic_cast<Image<uint16_t>*>(img);
			auto img32 = dynamic_cast<Image<uint32_t>*>(img);
			auto img64 = dynamic_cast<Image<uint64_t>*>(img);
			auto imgf32 = dynamic_cast<Image<float32_t>*>(img);
			auto imgC32 = dynamic_cast<Image<complex32_t>*>(img);

			Vec3c dimensions = img->dimensions();
			ImageDataType dt;
			double pixelSize = 0;

			if (img8)
			{
				dt = ImageDataType::UInt8;
				pixelSize = 1;
			}
			else if (img16)
			{
				dt = ImageDataType::UInt16;
				pixelSize = 2;
			}
			else if (img32)
			{
				dt = ImageDataType::UInt32;
				pixelSize = 4;
			}
			else if (img64)
			{
				dt = ImageDataType::UInt64;
				pixelSize = 8;
			}
			else if (imgf32)
			{
				dt = ImageDataType::Float32;
				pixelSize = 4;
			}
			else if (imgC32)
			{
				dt = ImageDataType::Complex32;
				pixelSize = 8;
			}
			else
			{
				dimensions = Vec3c(0, 0, 0);
				dt = ImageDataType::Unknown;
				pixelSize = 0;
			}
			
			double dataSize = dimensions.x * dimensions.y * dimensions.z * pixelSize;
			totalSize += dataSize;

			cout << name << ", " << dimensions << ", " << itl2::toString(dt) << ", " << bytesToString(dataSize) << endl;
		}
		cout << "Total " << bytesToString(totalSize) << endl;
	}

	vector<string> ListCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		double totalSize = 0;
		PISystem* system = distributor.getSystem();
		for (const string& name : system->getDistributedImageNames())
		{
			DistributedImageBase* img = system->getDistributedImage(name);
			auto img8 = dynamic_cast<DistributedImage<uint8_t>*>(img);
			auto img16 = dynamic_cast<DistributedImage<uint16_t>*>(img);
			auto img32 = dynamic_cast<DistributedImage<uint32_t>*>(img);
			auto img64 = dynamic_cast<DistributedImage<uint64_t>*>(img);
			auto imgf32 = dynamic_cast<DistributedImage<float32_t>*>(img);
			auto imgC32 = dynamic_cast<DistributedImage<complex32_t>*>(img);

			Vec3c dimensions = img->dimensions();
			ImageDataType dt;
			double pixelSize = 0;

			if (img8)
			{
				dt = ImageDataType::UInt8;
				pixelSize = 1;
			}
			else if (img16)
			{
				dt = ImageDataType::UInt16;
				pixelSize = 2;
			}
			else if (img32)
			{
				dt = ImageDataType::UInt32;
				pixelSize = 4;
			}
			else if (img64)
			{
				dt = ImageDataType::UInt64;
				pixelSize = 8;
			}
			else if (imgf32)
			{
				dt = ImageDataType::Float32;
				pixelSize = 4;
			}
			else if (imgC32)
			{
				dt = ImageDataType::Complex32;
				pixelSize = 8;
			}
			else
			{
				dimensions = Vec3c(0, 0, 0);
				dt = ImageDataType::Unknown;
				pixelSize = 0;
			}

			double dataSize = dimensions.x * dimensions.y * dimensions.z * pixelSize;
			totalSize += dataSize;

			cout << name << ", " << dimensions << ", " << itl2::toString(dt) << ", " << bytesToString(dataSize) << endl;
		}

		return vector<string>();
	}

	void NewImageCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		if (dt == ImageDataType::UInt8)
			system->replaceImage(name, new itl2::Image<uint8_t>(w, h, d));
		else if (dt == ImageDataType::UInt16)
			system->replaceImage(name, new itl2::Image<uint16_t>(w, h, d));
		else if (dt == ImageDataType::UInt32)
			system->replaceImage(name, new itl2::Image<uint32_t>(w, h, d));
		else if (dt == ImageDataType::UInt64)
			system->replaceImage(name, new itl2::Image<uint64_t>(w, h, d));
		else if (dt == ImageDataType::Float32)
			system->replaceImage(name, new itl2::Image<float32_t>(w, h, d));
		else if(dt == ImageDataType::Complex32)
			system->replaceImage(name, new itl2::Image<complex32_t>(w, h, d));
		else
			throw ParseException(string("Invalid data type: ") + dts);
	}

	vector<string> NewImageCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		PISystem* system = distributor.getSystem();

		if (dt == ImageDataType::UInt8)
			system->replaceDistributedImage(name, new DistributedImage<uint8_t>(name, w, h, d));
		else if (dt == ImageDataType::UInt16)
			system->replaceDistributedImage(name, new DistributedImage<uint16_t>(name, w, h, d));
		else if (dt == ImageDataType::UInt32)
			system->replaceDistributedImage(name, new DistributedImage<uint32_t>(name, w, h, d));
		else if (dt == ImageDataType::UInt64)
			system->replaceDistributedImage(name, new DistributedImage<uint64_t>(name, w, h, d));
		else if (dt == ImageDataType::Float32)
			system->replaceDistributedImage(name, new DistributedImage<float32_t>(name, w, h, d));
		else if (dt == ImageDataType::Complex32)
			system->replaceDistributedImage(name, new DistributedImage<complex32_t>(name, w, h, d));
		else
			throw ParseException(string("Invalid data type: ") + dts);

		return vector<string>();
	}

	/**
	Parse arguments from parameter vector.
	Throws exception if arguments cannot be parsed or if file dimension parsing is requested but cannot be done.
	*/
	void ReadRawCommand::parseArgs(vector<ParamVariant>& args, string& name, string& fname, coord_t& w, coord_t& h, coord_t& d, ImageDataType& dt)
	{
		name = pop<string>(args);
		fname = pop<string>(args);
		string dts = pop<string>(args);
		w = pop<coord_t>(args);
		h = pop<coord_t>(args);
		d = pop<coord_t>(args);

		dt = fromString<ImageDataType>(dts);

		// Parse dimensions from file name if no dimensions are provided
		if (w <= 0 || h <= 0 || d <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			if (!raw::getInfo(fname, dims, dt2))
				throw ParseException(string("Unable to find dimensions from file name: ") + fname);

			if (dt == ImageDataType::Unknown)
				dt = dt2;

			w = dims.x;
			h = dims.y;
			d = dims.z;
		}
	}


	void ReadCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);

		Vec3c dimensions;
		ImageDataType dt;
		if (!io::getInfo(filename, dimensions, dt))
			throw ITLException(string("File type cannot be automatically recognized: ") + filename);

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(dimensions);
			system->replaceImage(name, img);
			io::read(*img, filename);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	vector<string> ReadCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);

		Vec3c dimensions;
		ImageDataType dt;
		if (!io::getInfo(filename, dimensions, dt))
			throw ITLException(string("File type cannot be automatically recognized: ") + filename);

		PISystem* system = distributor.getSystem();

		if (dt == ImageDataType::UInt8)
		{
			DistributedImage<uint8_t>* img = new DistributedImage<uint8_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt16)
		{
			DistributedImage<uint16_t>* img = new DistributedImage<uint16_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt32)
		{
			DistributedImage<uint32_t>* img = new DistributedImage<uint32_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt64)
		{
			DistributedImage<uint64_t>* img = new DistributedImage<uint64_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Float32)
		{
			DistributedImage<float32_t>* img = new DistributedImage<float32_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Complex32)
		{
			DistributedImage<complex32_t>* img = new DistributedImage<complex32_t>(name, dimensions, filename);
			system->replaceDistributedImage(name, img);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));

		return vector<string>();
	}


	void ReadRawCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name;
		string filename;
		coord_t w;
		coord_t h;
		coord_t d;
		ImageDataType dt;
		parseArgs(args, name, filename, w, h, d, dt);

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(w, h, d);
			system->replaceImage(name, img);
			raw::read(*img, filename);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	vector<string> ReadRawCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name;
		string filename;
		coord_t w;
		coord_t h;
		coord_t d;
		ImageDataType dt;
		parseArgs(args, name, filename, w, h, d, dt);

		PISystem* system = distributor.getSystem();

		if (dt == ImageDataType::UInt8)
		{
			DistributedImage<uint8_t>* img = new DistributedImage<uint8_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt16)
		{
			DistributedImage<uint16_t>* img = new DistributedImage<uint16_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt32)
		{
			DistributedImage<uint32_t>* img = new DistributedImage<uint32_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt64)
		{
			DistributedImage<uint64_t>* img = new DistributedImage<uint64_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Float32)
		{
			DistributedImage<float32_t>* img = new DistributedImage<float32_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Complex32)
		{
			DistributedImage<complex32_t>* img = new DistributedImage<complex32_t>(name, w, h, d, filename);
			system->replaceDistributedImage(name, img);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));

		return vector<string>();
	}

	void ReadSequenceCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		//coord_t z0 = pop<coord_t>(args);
		//coord_t z1 = pop<coord_t>(args);

		//if (z1 < 0)
		//	z1 = numeric_limits<coord_t>::max();

		coord_t w, h, d;
		ImageDataType dt;
		sequence::getInfo(fname, w, h, d, dt);

		//math::clamp<coord_t>(z0, 0, d - 1);
		//math::clamp<coord_t>(z1, 0, d - 1);

		//d = z1 - z0 + 1;

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(w, h, d);
			system->replaceImage(name, img);
			sequence::read(*img, fname);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	vector<string> ReadSequenceCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		
		coord_t w, h, d;
		ImageDataType dt;
		sequence::getInfo(fname, w, h, d, dt);

		PISystem* system = distributor.getSystem();
		
		if (dt == ImageDataType::UInt8)
		{
			DistributedImage<uint8_t>* img = new DistributedImage<uint8_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt16)
		{
			DistributedImage<uint16_t>* img = new DistributedImage<uint16_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt32)
		{
			DistributedImage<uint32_t>* img = new DistributedImage<uint32_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::UInt64)
		{
			DistributedImage<uint64_t>* img = new DistributedImage<uint64_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Float32)
		{
			DistributedImage<float32_t>* img = new DistributedImage<float32_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else if (dt == ImageDataType::Complex32)
		{
			DistributedImage<complex32_t>* img = new DistributedImage<complex32_t>(name, w, h, d, fname);
			system->replaceDistributedImage(name, img);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));

		return vector<string>();
	}

	void ReadVolCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);

		Vec3c dimensions;
		ImageDataType dt;
		string endianness;
		size_t headerSize;
		if (!vol::getInfo(fname, dimensions, dt, endianness, headerSize))
			throw ITLException(string("Not a .vol file: ") + fname);

		coord_t w = dimensions.x;
		coord_t h = dimensions.y;
		coord_t d = dimensions.z;

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(w, h, d);
			system->replaceImage(name, img);
			vol::read(*img, fname);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	void ReadBlockCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		coord_t x = pop<coord_t>(args);
		coord_t y = pop<coord_t>(args);
		coord_t z = pop<coord_t>(args);
		coord_t bw = pop<coord_t>(args);
		coord_t bh = pop<coord_t>(args);
		coord_t bd = pop<coord_t>(args);

		Vec3c dims;
		ImageDataType dt;
		if (!io::getInfo(fname, dims, dt))
			throw ITLException(string("Unable to read file dimensions: ") + fname);

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			io::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	void ReadRawBlockCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		coord_t x = pop<coord_t>(args);
		coord_t y = pop<coord_t>(args);
		coord_t z = pop<coord_t>(args);
		coord_t bw = pop<coord_t>(args);
		coord_t bh = pop<coord_t>(args);
		coord_t bd = pop<coord_t>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		// Parse dimensions from file name if no dimensions are provided
		if (w <= 0 || h <= 0 || d <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			if (!raw::getInfo(fname, dims, dt2))
				throw ParseException(string("Unable to find dimensions from file name: ") + fname);

			if (dt == ImageDataType::Unknown)
				dt = dt2;

			w = dims.x;
			h = dims.y;
			d = dims.z;
		}

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			raw::readBlockNoParse(*img, fname, Vec3c(w, h, d), Vec3c(x, y, z), true);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	void ReadSequenceBlockCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		coord_t x = pop<coord_t>(args);
		coord_t y = pop<coord_t>(args);
		coord_t z = pop<coord_t>(args);
		coord_t bw = pop<coord_t>(args);
		coord_t bh = pop<coord_t>(args);
		coord_t bd = pop<coord_t>(args);

		coord_t w, h, d;
		ImageDataType dt;
		sequence::getInfo(fname, w, h, d, dt);

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(bw, bh, bd);
			system->replaceImage(name, img);
			sequence::readBlock(*img, fname, Vec3c(x, y, z), true);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	void MapRawCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		// Parse dimensions from file name if no dimensions are provided
		if (w <= 0 || h <= 0 || d <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			if (!raw::getInfo(fname, dims, dt2))
				throw ParseException(string("Unable to find dimensions from file name: ") + fname);

			if (dt == ImageDataType::Unknown)
				dt = dt2;

			w = dims.x;
			h = dims.y;
			d = dims.z;
		}

		if (dt == ImageDataType::UInt8)
		{
			itl2::Image<uint8_t>* img = new itl2::Image<uint8_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else if (dt == ImageDataType::UInt16)
		{
			itl2::Image<uint16_t>* img = new itl2::Image<uint16_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else if (dt == ImageDataType::UInt32)
		{
			itl2::Image<uint32_t>* img = new itl2::Image<uint32_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else if (dt == ImageDataType::UInt64)
		{
			itl2::Image<uint64_t>* img = new itl2::Image<uint64_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else if (dt == ImageDataType::Float32)
		{
			itl2::Image<float32_t>* img = new itl2::Image<float32_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else if (dt == ImageDataType::Complex32)
		{
			itl2::Image<complex32_t>* img = new itl2::Image<complex32_t>(fname, w, h, d);
			system->replaceImage(name, img);
		}
		else
			throw ParseException(string("Invalid data type: ") + itl2::toString(dt));
	}

	void ClearCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		if (name != "")
		{
			system->replaceImage(name, 0);
			system->replaceDistributedImage(name, 0);
		}
		else
		{
			for (auto name : system->getImageNames())
			{
				system->replaceImage(name, 0);
				system->replaceDistributedImage(name, 0);
			}
		}
	}

	vector<string> ClearCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();
		runInternal(system, args);
		//string name = pop<string>(args);
		//if (name != "")
		//{
		//	system->replaceDistributedImage(name, 0);
		//}
		//else
		//{
		//	for (auto name : system->getDistributedImageNames())
		//		system->replaceDistributedImage(name, 0);
		//}

		return vector<string>();
	}

	void HelpCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string term = pop<string>(args);

		if (name.length() > 0)
		{
			vector<string> topics = system->getHelp(name);

			if (topics.size() <= 0)
			{
				cout << "Command " << name << " not found." << endl;
				return;
			}

			int printedCount = 0;
			for (size_t n = 0; n < topics.size(); n++)
			{
				if (term.length() <= 0 || topics[n].find(term, 0) < topics[n].length())
				{
					if (n > 0)
						cout << endl;
					//	cout << "--------------------" << endl << endl;
					cout << topics[n];
					printedCount++;
				}
			}

			if (printedCount <= 0)
			{
				cout << "No topics contained term " << term << "." << endl;
			}
		}
		else
		{
			cout << "List of commands follows. Use help(command_name) to retrieve help for specific command." << endl << endl;

			cout << system->commandList(true) << endl;
		}

	}

	void EchoCommandsCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		bool echo = pop<bool>(args);
		bool timing = pop<bool>(args);
		system->showCommands(echo, timing);
	}

	void DistributeCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string provider = pop<string>(args);
		system->distributedProcessing(provider);
	}
}
