###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
# All right reserved.
#
# Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

#MAKE_SETTING ?= ../make_setting

#include $(MAKE_SETTING)

TARGET = libFB.a

CPMLIB_DIR = ../CPMlib
PMLIB_DIR = ../PMlib-1.4

AR ?= ar
ARFLAGS ?= rv
RANLIB ?= ranlib

CXXSRCS = \
  FBUtility.C \
  BndOuter.C \
  VoxInfo.C \
  History.C \
  dfi.C \
  FileIO.C

F90SRCS = \
  FB_util.f90

SRCS  = $(CXXSRCS) $(F90SRCS)

.SUFFIXES: .o .C .f90

CXXOBJS = $(CXXSRCS:.C=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(CXXOBJS) $(F90OBJS)

$(TARGET):$(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.C.o:
	$(CXX) $(CXXFLAGS) -I$(CPMLIB_DIR)/include -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

#depend:
#	makedepend -I. -Y $(SRCS)

# DO NOT DELETE

FBUtility.o: FBUtility.h FB_Define.h mydebug.h
BndOuter.o: BndOuter.h FB_Define.h mydebug.h
VoxInfo.o: VoxInfo.h DomainInfo.h FB_Define.h mydebug.h FBUtility.h
VoxInfo.o: Component.h Medium.h BndOuter.h
History.o: History.h ItrCtl.h FB_Define.h mydebug.h Component.h FBUtility.h
dfi.o: dfi.h FB_Define.h mydebug.h
FileIO.o: FileIO.h DomainInfo.h FB_Define.h mydebug.h FB_Ffunc.h
