#
# Makefile for GENIE test programs
#
# Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

#gtestDecay		 \
#gtestHadronization	 \

TGT = gAnalysis

all: $(TGT)


gAnalysis: FORCE
	$(CXX) $(CXXFLAGS) -c genie_analysis.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) genie_analysis.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gAnalysis

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gAnalysis 	

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gAnalysis 	

FORCE:

