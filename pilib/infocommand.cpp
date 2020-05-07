

#include "infocommand.h"
#include "utilities.h"
#include "commandlist.h"

#include <iostream>

using namespace std;
using namespace itl2;

namespace pilib
{
	void addInfoCommands()
	{
		CommandList::add<InfoCommand>();
	}

	void InfoCommand::run(vector<ParamVariant>& args) const
	{
		cout << "Process Image = pi2" << endl;
		cout << "Copyright(C) 2011- Arttu Miettinen" << endl;
		cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
		cout << "This is free software, and you are welcome to redistribute it" << endl;
		cout << "under certain conditions; run `license()' command for more details." << endl;
		cout << endl;
		cout << "Contact: arttu.i.miettinen@jyu.fi" << endl;
		cout << endl;
		cout << "Based on work done at" << endl;
		cout << "Complex materials research group, Department of Physics, University of Jyvaskyla, Finland" << endl;
		cout << "X-ray tomography research group, TOMCAT beamline, Swiss Light Source (SLS), Paul Scherrer Institute (PSI), Switzerland" << endl;
		cout << "Institute for Biomedical Engineering, Swiss Federal Institute of Technology Zurich (ETH Zurich), Switzerland" << endl;
		cout << "Centre d'Imagerie BioMedicale (CIBM), Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland" << endl;
		cout << endl;

		if (sizeof(void*) == 4)
			cout << "32-bit version" << endl;
		else if (sizeof(void*) == 8)
			cout << "64-bit version" << endl;
		else
			cout << "Warning: Using neither 32-bit nor 64-bit version." << endl;

#if defined(NO_OPENCL)
		cout << "OpenCL support is disabled." << endl;
#else
		cout << "OpenCL support is enabled." << endl;
#endif

		cout << "Number of threads: " << omp_get_max_threads() << endl;
		cout << "Available RAM: " << bytesToString((double)memorySize()) << endl;

		// This defines VERSION variable
		#include "commit_info.txt"

		cout << "Version: " << VERSION << endl;
	}
}