
#include "specialcommands.h"
#include "pisystem.h"
#include "distributedimage.h"
#include "utilities.h"
#include "stringutils.h"
#include <omp.h>
#include "io/vol.h"
#include "io/io.h"
#include "pilibutilities.h"
#include "whereamicpp.h"
#include "commandmacros.h"
#include "timing.h"

using namespace std;

namespace pilib
{
	void addSpecialCommands()
	{
		CommandList::add<ClearCommand>();
		CommandList::add<DistributeCommand>();
		CommandList::add<MaxMemoryCommand>();
		CommandList::add<GetMaxMemoryCommand>();
		CommandList::add<MaxJobsCommand>();
		CommandList::add<ChunkSizeCommand>();
		CommandList::add<DelayingCommand>();
		CommandList::add<PrintTaskScriptsCommand>();
		CommandList::add<EchoCommandsCommand>();
		CommandList::add<HelloCommand>();
		CommandList::add<PrintCommand>();
		CommandList::add<HelpCommand>();
		CommandList::add<CommandReferenceCommand>();
		CommandList::add<ListCommand>();
		CommandList::add<LicenseCommand>();
		CommandList::add<ReadCommand>();
		CommandList::add<ReadSequenceCommand>();
		CommandList::add<ReadVolCommand>();
		CommandList::add<WaitReturnCommand>();

		CommandList::add<TimingCommand>();
		CommandList::add<SaveTimingCommand>();
		CommandList::add<ResetTimingCommand>();

		CommandList::add<MapRawCommand>();
		CommandList::add<MapRaw2Command>();

		CommandList::add<NewValueCommand>();
		CommandList::add<SetStringCommand>();
		CommandList::add<SetIntCommand>();
		CommandList::add<SetRealCommand>();
		CommandList::add<SetBoolCommand>();

		CommandList::add<NewImageCommand>();
		CommandList::add<NewImage2Command>();

		CommandList::add<ReadRawCommand>();
		CommandList::add<ReadRaw2Command>();

		CommandList::add<ReadBlockCommand>();
		CommandList::add<ReadBlock2Command>();

		CommandList::add<ReadRawBlockCommand>();
		CommandList::add<ReadRawBlock2Command>();

		CommandList::add<ReadSequenceBlockCommand>();
		CommandList::add<ReadSequenceBlock2Command>();

		CommandList::add<ReadNN5BlockCommand>();
		CommandList::add<ReadNN5Block2Command>();
		
		ADD_ALL(EnsureSizeCommand);
		ADD_ALL(EnsureSize2Command);

		ADD_ALL(GetMapFileCommand);
	}




	void LicenseCommand::run(vector<ParamVariant>& args) const
	{
		cout << "This program is licensed under GNU General Public License Version 3:" << endl;

		fs::path path = getModulePath();
		string gpl3 = readText(path.replace_filename("LICENSE.txt").string(), false);
		if (gpl3.length() > 0)
		{
			cout << endl << gpl3 << endl;
		}
		else
		{
			cout << "Please see LICENSE.txt bundled with this software or" << endl;
			cout << "https://www.gnu.org/licenses/gpl-3.0.html for full license text." << endl;
		}


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This software is based in part on the work of the Independent JPEG Group." << endl;
		cout << endl;


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from the LZ4 library, licensed" << endl;
		cout << "under the license shown below." << endl;
		cout << "See also https://github.com/lz4/lz4" << endl;
		cout << R"END(
LZ4 Library
Copyright (c) 2011-2020, Yann Collet
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
)END" << endl;
		cout << endl;




		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from the JSON for Modern C++ library, licensed" << endl;
		cout << "under the MIT license available at http://www.opensource.org/licenses/MIT." << endl;
		cout << "See also https://github.com/nlohmann/json" << endl;
		cout << R"END(
JSON for Modern C++
version 3.10.5
https://github.com/nlohmann/json

Licensed under the MIT License < http ://opensource.org/licenses/MIT>.
SPDX - License - Identifier : MIT
Copyright(c) 2013 - 2022 Niels Lohmann < http ://nlohmann.me>.

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this softwareand associated  documentation files(the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy, modify, merge, publish, distribute, sublicense, and /or  sell
copies  of  the Software, and to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions :

The above copyright noticeand this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT.IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM, DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
)END" << endl;
		cout << endl;


		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from the C++ Mathematical Expression Toolkit Library, licensed" << endl;
		cout << "under the MIT license available at http://www.opensource.org/licenses/MIT." << endl;
		cout << "See also http://www.partow.net/programming/exprtk/index.html" << endl;
		cout << R"END(
Copyright 2021 Arash Partow

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files(the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
)END" << endl;
		cout << endl;







		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code from PStreams library, licensed under the following license:" << endl;
		cout << R"END(
PStreams - POSIX Process I/O for C++

Copyright (C) 2001 - 2017 Jonathan Wakely
Distributed under the Boost Software License, Version 1.0.
	
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license(the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third - parties to whom the Software is furnished to
do so, all subject to the following :

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine - executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON - INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
)END" << endl;
		cout << endl;






		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code ported from public-domain reference implementation of JAMA : A Java Matrix Package." << endl;
		cout << "The original code is available at https://math.nist.gov/javanumerics/jama/" << endl;
		cout << endl;




		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains code ported from a public-domain library TNT : Template Numerical Toolkit." << endl;
		cout << "The original code is available at https://math.nist.gov/tnt/" << endl;
		cout << endl;




		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains optimized median calculation code by N. Devillard." << endl;
		cout << "The original code is downloadable from http://ndevilla.free.fr/median/ and is licensed under the following statement:" << endl;
		cout << R"END(
The following routines have been built from knowledge gathered
around the Web. I am not aware of any copyright problem with
them, so use it as you want.
N. Devillard - 1998
)END" << endl;
		cout << endl;




		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "This program contains 'whereami' library by Gregory Pakosz, licened under the MIT license:" << endl;
		cout << R"END(
Copyright Gregory Pakosz

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files(the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
)END" << endl;
		cout << endl;




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

* Redistributions of source code must retain the above copyright notice, this
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
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
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
warranty. In no event will the authors be held liable for any damages
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

Jean-loup Gailly        Mark Adler
jloup@gzip.org          madler@alumni.caltech.edu

If you use the zlib library in a product, we would appreciate *not* receiving
lengthy legal documents to sign. The sources are provided for free but without
warranty of any kind. The library has been entirely written by Jean-loup
Gailly and Mark Adler; it does not include third-party code.

If you redistribute modified sources, we would appreciate that you include in
the file ChangeLog history information documenting your changes. Please read
the FAQ for more information on the distribution of modified source versions.)END" << endl;
		cout << endl;





		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl << endl;
		cout << "Windows versions of this program are shipped with FFTW3 library licensed" << endl;
		cout << "under the GNU General Public License Version 2." << endl;
		cout << "Please see https://www.gnu.org/licenses/gpl-2.0.html for full license text." << endl;


	}

	void TimingCommand::run(vector<ParamVariant>& args) const
	{
		cout << GlobalTimer::results().toString() << endl;
	}

	void SaveTimingCommand::run(vector<ParamVariant>& args) const
	{
		string fname = pop<string>(args);
		GlobalTimer::results().toFile(fname);
	}

	void ResetTimingCommand::run(vector<ParamVariant>& args) const
	{
		GlobalTimer::reset();
	}

	void ListCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		cout << "Images:" << endl;
		cout << "-------" << endl;

		double totalSize = 0;
		for (const string& name : system->getImageNames())
		{
			ImageBase* img = system->getImage(name);
			
			Vec3c dimensions = img->dimensions();
			size_t pixelSize = img->pixelSize();
			size_t dataSize = img->pixelCount() * pixelSize;
			totalSize += dataSize;

			cout << name << ", " << dimensions << ", " << itl2::toString(img->dataType()) << ", " << bytesToString((double)dataSize) << endl;
		}
		if (totalSize > 0)
			cout << "Total " << bytesToString(totalSize) << endl;
		else
			cout << "-- none --" << endl;

		cout << "Variables:" << endl;
		cout << "----------" << endl;
		size_t count = 0;
		auto v = system->getValueNames();
		for (const string& name : v)
		{
			count++;
			Value* val = system->getValue(name);
			cout << name << " (" << pilib::toString(val->getType()) << ")" << endl;
		}

		if (count > 0)
			cout << "Total " << count << " variables." << endl;
		else
			cout << "-- none --" << endl;

	}

	vector<string> ListCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		double totalSize = 0;
		PISystem* system = distributor.getSystem();
		for (const string& name : system->getDistributedImageNames())
		{
			DistributedImageBase* img = system->getDistributedImage(name);

			Vec3c dimensions = img->dimensions();
			size_t pixelSize = img->pixelSize();
			size_t dataSize = img->pixelCount() * pixelSize;
			totalSize += dataSize;

			cout << name << ", " << dimensions << ", " << itl2::toString(img->dataType()) << ", " << bytesToString((double)dataSize) << endl;
		}

		return vector<string>();
	}

	

	void NewValueCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		string value = pop<string>(args);

		ValueType type;
		if (dts == "string")
			type = ValueType::String;
		else if (dts == "int" || dts == "integer")
			type = ValueType::Int;
		else if (dts == "real" || dts == "float" || dts == "double")
			type = ValueType::Real;
		else if (dts == "bool" || dts == "boolean")
			type = ValueType::Bool;
		else
			throw ITLException(string("Unsupported variable type: ") + dts);

		shared_ptr<Value> pival = make_shared<Value>(type);
		
		if (value != "")
		{
			if (dts == "string")
				pival->stringValue = value;
			else if (dts == "int" || dts == "integer")
				pival->intValue = fromString<coord_t>(value);
			else if (dts == "real" || dts == "float" || dts == "double")
				pival->realValue = fromString<double>(value);
			else if (dts == "bool" || dts == "boolean")
				pival->boolValue = fromString<bool>(value);
			else
				throw ITLException(string("Unsupported variable type: ") + dts);
		}

		system->replaceValue(name, pival);
	}

	vector<string> NewValueCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();

		runInternal(system, args);

		return vector<string>();
	}


	void SetStringCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string* name = pop<string*>(args);
		string value = pop<string>(args);

		*name = value;
	}

	vector<string> SetStringCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();

		runInternal(system, args);

		return vector<string>();
	}



	void SetIntCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		coord_t* name = pop<coord_t*>(args);
		coord_t value = pop<coord_t>(args);

		*name = value;
	}

	vector<string> SetIntCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();

		runInternal(system, args);

		return vector<string>();
	}



	void SetRealCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		double* name = pop<double*>(args);
		double value = pop<double>(args);

		*name = value;
	}

	vector<string> SetRealCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();

		runInternal(system, args);

		return vector<string>();
	}



	void SetBoolCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		bool* name = pop<bool*>(args);
		bool value = pop<bool>(args);

		*name = value;
	}

	vector<string> SetBoolCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();

		runInternal(system, args);

		return vector<string>();
	}



	void NewImageCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		Vec3c dimensions(w, h, d);

		pick<CreateImage>(dts, dimensions, name, system);
	}

	vector<string> NewImageCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		Vec3c dimensions(w, h, d);
		PISystem* system = distributor.getSystem();

		pick<CreateEmptyDistributedImage>(dts, name, dimensions, system);

		return vector<string>();
	}

	void NewImage2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dimensions = pop<Vec3c>(args);

		pick<CreateImage>(dts, dimensions, name, system);
	}

	vector<string> NewImage2Command::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dimensions = pop<Vec3c>(args);

		PISystem* system = distributor.getSystem();

		pick<CreateEmptyDistributedImage>(dts, name, dimensions, system);

		return vector<string>();
	}




	template<typename pixel_t> struct CreateDistributedImageFromFile
	{
		static void run(const string& imgName, const Vec3c& dimensions, const string& filename, PISystem* system)
		{
			shared_ptr<DistributedImageBase> ptr = make_shared<DistributedImage<pixel_t> >(*system->getDistributor(), imgName, dimensions, filename);
			system->replaceDistributedImage(imgName, ptr);
		}
	};
	


	template<typename pixel_t> struct CreateImageAndRead
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			io::read<pixel_t>(img, filename);
		}
	};

	void ReadCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		Vec3c dimensions;
		ImageDataType dt2;
		string reason;
		if (!io::getInfo(filename, dimensions, dt2, reason))
			throw ITLException(string("File type cannot be automatically recognized, the given template does not uniquely identify any file, or the file is not found: ") + filename + ". \n" + reason);

		if (dt == ImageDataType::Unknown)
			dt = dt2;
		
		pick<CreateImageAndRead>(dt, dimensions, name, system, filename);
	}


	vector<string> ReadCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		Vec3c dimensions;
		ImageDataType dt2;
		string reason;
		if (!io::getInfo(filename, dimensions, dt2, reason))
			throw ITLException(string("File type cannot be automatically recognized, the given template does not uniquely identify any file, or the file is not found: ") + filename + ". \n" + reason);

		if (dt == ImageDataType::Unknown)
			dt = dt2;

		pick<CreateDistributedImageFromFile>(dt, name, dimensions, filename, distributor.getSystem());

		return vector<string>();
	}

	


	///**
	//Parse arguments from parameter vector.
	//Throws exception if arguments cannot be parsed or if file dimension parsing is requested but cannot be done.
	//*/
	//void ReadRawCommand::parseArgs(vector<ParamVariant>& args, string& name, string& fname, coord_t& w, coord_t& h, coord_t& d, ImageDataType& dt)
	//{
	//	name = pop<string>(args);
	//	fname = pop<string>(args);
	//	string dts = pop<string>(args);
	//	w = pop<coord_t>(args);
	//	h = pop<coord_t>(args);
	//	d = pop<coord_t>(args);

	//	dt = fromString<ImageDataType>(dts);

	//	// Parse dimensions from file name if no dimensions are provided
	//	if (w <= 0 || h <= 0 || d <= 0 || dt == ImageDataType::Unknown)
	//	{
	//		Vec3c dims;
	//		ImageDataType dt2;
	//		string reason;
	//		if (!raw::getInfo(fname, dims, dt2, reason))
	//			throw ParseException(string("Unable to find dimensions from file name or file not found: ") + fname + ". " + reason);

	//		// Not required in this command. This is done later in raw::read.
	//		//internals::expandRawFilename(filename);

	//		if (dt == ImageDataType::Unknown)
	//			dt = dt2;

	//		w = dims.x;
	//		h = dims.y;
	//		d = dims.z;
	//	}
	//}

	/**
	Parse arguments from parameter vector.
	Throws exception if arguments cannot be parsed or if file dimension parsing is requested but cannot be done.
	*/
	void parseArgs(string& name, string& fname, Vec3c& dimensions, ImageDataType& dt)
	{
		raw::internals::expandRawFilename(fname);

		// Parse dimensions from file name if no dimensions are provided
		if (dimensions.x <= 0 || dimensions.y <= 0 || dimensions.z <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			string reason;
			if (!raw::getInfo(fname, dims, dt2, reason))
				throw ParseException(string("Unable to find dimensions from file name or file not found: ") + fname + ". " + reason);

			if (dt == ImageDataType::Unknown)
				dt = dt2;

			dimensions = dims;
		}
	}

	template<typename pixel_t> struct CreateImageAndReadRaw
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			raw::readNoParse<pixel_t>(img, filename);
		}
	};

	void ReadRawCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		Vec3c dimensions(w, h, d);
		ImageDataType dt = fromString<ImageDataType>(dts);

		parseArgs(name, filename, dimensions, dt);

		pick<CreateImageAndReadRaw>(dt, dimensions, name, system, filename);
	}

	vector<string> ReadRawCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		Vec3c dimensions(w, h, d);
		ImageDataType dt = fromString<ImageDataType>(dts);

		parseArgs(name, filename, dimensions, dt);

		pick<CreateDistributedImageFromFile>(dt, name, dimensions, filename, distributor.getSystem());

		return vector<string>();
	}




	void ReadRaw2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dimensions = pop<Vec3c>(args);
		
		ImageDataType dt = fromString<ImageDataType>(dts);

		parseArgs(name, filename, dimensions, dt);

		pick<CreateImageAndReadRaw>(dt, dimensions, name, system, filename);
	}

	vector<string> ReadRaw2Command::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dimensions = pop<Vec3c>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		parseArgs(name, filename, dimensions, dt);

		pick<CreateDistributedImageFromFile>(dt, name, dimensions, filename, distributor.getSystem());

		return vector<string>();
	}










	template<typename pixel_t> struct CreateImageAndReadSequence
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			sequence::read<pixel_t>(img, filename);
		}
	};

	void ReadSequenceCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);

		coord_t w, h, d;
		ImageDataType dt;
		string reason;
		if (!sequence::getInfo(filename, w, h, d, dt, reason))
			throw ITLException(string("Unable to recognize sequence: ") + filename + ". " + reason);

		Vec3c dimensions(w, h, d);

		pick<CreateImageAndReadSequence>(dt, dimensions, name, system, filename);
	}

	vector<string> ReadSequenceCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		
		coord_t w, h, d;
		ImageDataType dt;
		string reason;
		if (!sequence::getInfo(filename, w, h, d, dt, reason))
			throw ITLException(string("Unable to recognize sequence: ") + filename + ". " + reason);

		Vec3c dimensions(w, h, d);

		pick<CreateDistributedImageFromFile>(dt, name, dimensions, filename, distributor.getSystem());

		return vector<string>();
	}



	template<typename pixel_t> struct CreateImageAndReadVol
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			vol::read<pixel_t>(img, filename);
		}
	};

	void ReadVolCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);

		Vec3c dimensions;
		ImageDataType dt;
		string endianness;
		size_t headerSize;
		string reason;
		if (!vol::getInfo(filename, dimensions, dt, endianness, headerSize, reason))
			throw ITLException(string("Not a .vol file: ") + filename + ". " + reason);

		pick<CreateImageAndReadVol>(dt, dimensions, name, system, filename);
	}




	template<typename pixel_t> struct CreateImageAndReadBlock
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename, const Vec3c& blockStart)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			io::readBlock<pixel_t>(img, filename, blockStart, true);
		}
	};

	void readBlock(string name, string fname, Vec3c blockStart, Vec3c blockDims, string dts, PISystem* system)
	{
		ImageDataType dt = fromString<ImageDataType>(dts);

		Vec3c dims;
		ImageDataType dt2;
		string reason;
		if (!io::getInfo(fname, dims, dt2, reason))
			throw ITLException(string("File type cannot be automatically recognized, the given template does not uniquely identify any file, or the file is not found: ") + fname + ". \n" + reason);

		if (dt == ImageDataType::Unknown)
			dt = dt2;

		pick<CreateImageAndReadBlock>(dt, blockDims, name, system, fname, blockStart);
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
		string dts = pop<string>(args);

		Vec3c blockStart(x, y, z);
		Vec3c blockDims(bw, bh, bd);

		readBlock(name, fname, blockStart, blockDims, dts, system);
	}

	void ReadBlock2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		Vec3c blockStart = pop<Vec3c>(args);
		Vec3c blockDims = pop<Vec3c>(args);
		string dts = pop<string>(args);

		readBlock(name, fname, blockStart, blockDims, dts, system);
	}




	template<typename pixel_t> struct CreateImageAndReadRawBlock
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename, const Vec3c& fileDims, const Vec3c& blockStart)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			raw::readBlockNoParse<pixel_t>(img, filename, fileDims, blockStart, true);
		}
	};

	void readRawBlock(string name, string fname, const Vec3c& blockStart, const Vec3c& blockDims, Vec3c fileDims, ImageDataType dt, PISystem* system)
	{
		// Parse dimensions from file name if no dimensions are provided
		if (fileDims.x <= 0 || fileDims.y <= 0 || fileDims.z <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			string reason;
			if (!raw::getInfo(fname, dims, dt2, reason))
				throw ParseException(string("Unable to read file: ") + fname + ". " + reason);

			raw::internals::expandRawFilename(fname);

			if (dt == ImageDataType::Unknown)
				dt = dt2;

			fileDims.x = dims.x;
			fileDims.y = dims.y;
			fileDims.z = dims.z;
		}

		pick<CreateImageAndReadRawBlock>(dt, blockDims, name, system, fname, fileDims, blockStart);
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

		readRawBlock(name, fname, Vec3c(x, y, z), Vec3c(bw, bh, bd), Vec3c(w, h, d), dt, system);
	}

	void ReadRawBlock2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		Vec3c blockStart = pop<Vec3c>(args);
		Vec3c blockDims = pop<Vec3c>(args);
		string dts = pop<string>(args);
		Vec3c fileDims = pop<Vec3c>(args);

		ImageDataType dt = fromString<ImageDataType>(dts);

		readRawBlock(name, fname, blockStart, blockDims, fileDims, dt, system);
	}





	template<typename pixel_t> struct CreateImageAndReadSequenceBlock
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename, const Vec3c& blockStart)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			sequence::readBlock<pixel_t>(img, filename, blockStart, true);
		}
	};

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
		string reason;
		if(!sequence::getInfo(fname, w, h, d, dt, reason))
			throw ParseException(string("Unable to read sequence: ") + fname + ". " + reason);

		Vec3c blockDims(bw, bh, bd);

		Vec3c blockStart(x, y, z);

		pick<CreateImageAndReadSequenceBlock>(dt, blockDims, name, system, fname, blockStart);
	}

	void ReadSequenceBlock2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		Vec3c blockStart = pop<Vec3c>(args);
		Vec3c blockDims = pop<Vec3c>(args);
		
		coord_t w, h, d;
		ImageDataType dt;
		string reason;
		if(!sequence::getInfo(fname, w, h, d, dt, reason))
			throw ParseException(string("Unable to read sequence: ") + fname + ". " + reason);

		pick<CreateImageAndReadSequenceBlock>(dt, blockDims, name, system, fname, blockStart);
	}





	template<typename pixel_t> struct CreateImageAndReadNN5Block
	{
		static void run(const Vec3c& dimensions, const string& imgName, PISystem* system, const string& filename, const Vec3c& blockStart)
		{
			Image<pixel_t>& img = *CreateImage<pixel_t>::run(dimensions, imgName, system);
			nn5::readBlock<pixel_t>(img, filename, blockStart);
		}
	};

	void ReadNN5BlockCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		coord_t x = pop<coord_t>(args);
		coord_t y = pop<coord_t>(args);
		coord_t z = pop<coord_t>(args);
		coord_t bw = pop<coord_t>(args);
		coord_t bh = pop<coord_t>(args);
		coord_t bd = pop<coord_t>(args);

		Vec3c dimensions;
		ImageDataType dt;
		string reason;
		if (!nn5::getInfo(fname, dimensions, dt, reason))
			throw ParseException(string("Unable to read NN5 dataset: ") + fname + ". " + reason);

		Vec3c blockDims(bw, bh, bd);

		Vec3c blockStart(x, y, z);

		pick<CreateImageAndReadNN5Block>(dt, blockDims, name, system, fname, blockStart);
	}

	void ReadNN5Block2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		Vec3c blockStart = pop<Vec3c>(args);
		Vec3c blockDims = pop<Vec3c>(args);

		Vec3c dimensions;
		ImageDataType dt;
		string reason;
		if (!nn5::getInfo(fname, dimensions, dt, reason))
			throw ParseException(string("Unable to read NN5 dataset: ") + fname + ". " + reason);

		pick<CreateImageAndReadNN5Block>(dt, blockDims, name, system, fname, blockStart);
	}



	template<typename pixel_t> struct CreateImageAndMapRaw
	{
		static void run(const Vec3c& dimensions, const string& imgName, const string& filename, bool readOnly, PISystem* system)
		{
			//shared_ptr<ParamVariant> img = make_shared<ParamVariant>(new Image<pixel_t>(filename, readOnly, dimensions));
			//system->replaceImage(imgName, img);
			shared_ptr<ImageBase> img = make_shared<Image<pixel_t>>(filename, readOnly, dimensions);
			system->replaceImage(imgName, img);
		}
	};


	void mapRaw(const string& name, string fname, const string& dts, coord_t w, coord_t h, coord_t d, bool readOnly, PISystem* system)
	{
		ImageDataType dt = fromString<ImageDataType>(dts);

		// Parse dimensions from file name if no dimensions are provided
		if (w <= 0 || h <= 0 || d <= 0 || dt == ImageDataType::Unknown)
		{
			Vec3c dims;
			ImageDataType dt2;
			string reason;
			if (raw::getInfo(fname, dims, dt2, reason))
			{
				// File is found, use that.
				raw::internals::expandRawFilename(fname);

				if (dt == ImageDataType::Unknown)
					dt = dt2;

				w = dims.x;
				h = dims.y;
				d = dims.z;

				// Make sure that the fname is just the prefix.
				fname = itl2::getPrefix(fname);
			}
			else
			{
				// No file found, create new image with (1, 1, 1) dimensions and fname as prefix.
				w = 1;
				h = 1;
				d = 1;
				//throw ParseException(string("Unable to find dimensions from file name or file not found: ") + fname);
			}
		}
		else
		{
			// Dimensions are given, the file name is assumed to be a prefix.
			// Everything should be ok.
		}

		Vec3c dimensions(w, h, d);

		pick<CreateImageAndMapRaw>(dt, dimensions, name, fname, readOnly, system);
	}


	void MapRawCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);
		bool readOnly = pop<bool>(args);

		mapRaw(name, fname, dts, w, h, d, readOnly, system);
	}

	void MapRaw2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string fname = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dims = pop<Vec3c>(args);
		bool readOnly = pop<bool>(args);

		mapRaw(name, fname, dts, dims.x, dims.y, dims.z, readOnly, system);
	}



	void ClearCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		if (name != "")
		{
			system->replaceImage(name, nullptr);
			system->replaceDistributedImage(name, nullptr);
			system->replaceValue(name, nullptr);
		}
		else
		{
			for (auto name : system->getImageNames())
			{
				system->replaceImage(name, nullptr);
				system->replaceDistributedImage(name, nullptr);
			}
			for (auto name : system->getValueNames())
				system->replaceValue(name, nullptr);
		}
	}

	vector<string> ClearCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		PISystem* system = distributor.getSystem();
		runInternal(system, args);

		return vector<string>();
	}

	string HelpCommand::run(const string& name, HelpFormat format) const
	{
		vector<string> topics = CommandList::help(name, format);

		stringstream out;

		if (format == HelpFormat::Rst)
		{
			out << ".. _" << name << ":" << endl << endl;
			printTitle(out, name, 0);

			if (topics.size() > 1)
				out << "There are " << topics.size() << " forms of this command." << endl << endl;
		}

		for (size_t n = 0; n < topics.size(); n++)
		{
			if (n > 0)
				out << endl;
			out << topics[n];
		}

		return out.str();
	}

	void HelpCommand::run(vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		//string term = pop<string>(args);
		string sformat = pop<string>(args);

		// TODO: Do this with HelpFormat fromString(const string& str)
		HelpFormat format = HelpFormat::Text;
		toLower(sformat);
		if (sformat == "rst" || sformat == "restructuredtext")
			format = HelpFormat::Rst;


		if (name.length() > 0)
		{
			string help = run(name, format);

			if (help.length() > 0)
			{
				cout << help << endl;
			}
			else
			{
				cout << "Command " << name << " not found." << endl;

				vector<string> apro = CommandList::apropos(name);
				if (apro.size() > 0)
				{
					cout << "Maybe one of these commands is related:" << endl;
					for (size_t n = 0; n < apro.size(); n++)
					{
						cout << apro[n] << endl;
					}
					cout << "Use help(command_name) to retrieve help for specific command." << endl;
				}

				return;
			}
		}
		else
		{
			cout << "List of commands follows. Use help(command_name) to retrieve help for specific command." << endl << endl;

			cout << CommandList::list(true) << endl;
		}

	}

	void EchoCommandsCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		bool echo = pop<bool>(args);
		bool timing = pop<bool>(args);
		system->showCommands(echo, timing);
	}

	void DelayingCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		bool enable = pop<bool>(args);
		if(system->getDistributor())
			system->getDistributor()->delaying(enable);
	}

	void PrintTaskScriptsCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		bool enable = pop<bool>(args);
		if (system->getDistributor())
			system->getDistributor()->showScripts(enable);
	}

	void MaxMemoryCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		double maxMem = pop<double>(args);
		if (system->getDistributor())
			system->getDistributor()->allowedMemory(itl2::round(maxMem * 1024.0 * 1024.0));
	}

	void GetMaxMemoryCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		double* maxMem = pop<double*>(args);
		if (system->getDistributor())
			*maxMem = system->getDistributor()->allowedMemory() / (1024.0 * 1024.0);
		else
			*maxMem = 0;
	}

	void MaxJobsCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		size_t n = pop<size_t>(args);
		if (system->getDistributor())
			system->getDistributor()->maxJobs(n);
	}

	void ChunkSizeCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		Vec3c chunkSize = pop<Vec3c>(args);
		if (system->getDistributor())
			system->getDistributor()->chunkSize(chunkSize);
	}

	void DistributeCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string provider = pop<string>(args);
		system->distributedProcessing(provider);
	}


	void CommandReferenceCommand::run(vector<ParamVariant>& args) const
	{
		string path = pop<string>(args);

		// Get names of commands (without args)
		vector<string> names;
		const vector<unique_ptr<Command> >& commands = CommandList::all();
		for (size_t n = 0; n < commands.size(); n++)
		{
			if (!commands[n]->isInternal())
				names.push_back(commands[n]->name());
		}
		removeDuplicates(names);
		sort(names.begin(), names.end());

		for(const auto& name : names)
		{
			string help = CommandList::get<HelpCommand>().run(name, HelpFormat::Rst);
			if (help.length() > 0)
			{
				fs::path p = path;
				writeText((p / (name + ".txt")).string(), help);
			}
			else
			{
				cout << "WARNING: No help for command " << name << endl;
			}
		}


	}
}
