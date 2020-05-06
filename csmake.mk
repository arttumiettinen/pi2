
# This file is used to interoperate between GNU makefiles for C++ compilation
# and xbuild/msbuild for C# compilation.

# When including this file, variables PROJECT_FILE, TARGET_FILE and CONFIG
# must be declared.


# Convert CONFIG to CS_CONFIG
ifeq ($(CONFIG), release)
    CS_CONFIG = Release
else ifeq ($(CONFIG), release-nocl)
    CS_CONFIG = Release no OpenCL
endif


#change this to the depth of the project folders
#if needed, add a prefix for a common project folder
CSHARP_SOURCE_FILES := $(wildcard *.cs *.resx *.ico Properties/*.cs Properties/*.resx)

OUTDIR:=../x64/$(CS_CONFIG)
OUTFILE:=$(OUTDIR)/$(TARGET_FILE)

# This is ugly hack to convert spaces " " to "\ "
empty:=
space:= $(empty) $(empty)
slashspace:=\ $(empty)
OUTDIR:=$(subst $(space),$(slashspace),$(OUTDIR))
OUTFILE:=$(subst $(space),$(slashspace),$(OUTFILE))

all: $(OUTFILE)

$(OUTFILE): $(CSHARP_SOURCE_FILES)
ifneq ($(shell which xbuild), )
	xbuild /p:Configuration="$(CS_CONFIG)" $(PROJECT_FILE)
	mkdir -p ../x64/$(CONFIG)
	cp $(OUTDIR)/* ../x64/$(CONFIG)
endif

clean:
ifneq ($(shell which xbuild), )
	xbuild /target:clean $(PROJECT_FILE)
endif

