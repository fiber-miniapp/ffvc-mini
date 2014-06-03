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

FFVC_VERSION = 0.9.2
FFVC_MINI_VERSION = 1.0.0

#MAKE_SETTING ?= ../make_setting

#include $(MAKE_SETTING)

CPMLIB_DIR = ../CPMlib

TARGET = ../../bin/ffvc_mini

CXXSRCS = \
    main.C \
    ffv.C \
    ffv_Initialize.C \
    ffv_Alloc.C \
    ffv_SetBC.C \
    ffv_Loop.C \
    ffv_LS.C \
    ffv_Post.C \
    NS_FS_E_Binary.C

SRCS = $(CXXSRCS)

.SUFFIXES: .o .C
CXXOBJS = $(CXXSRCS:.C=.o)
OBJS  = $(CXXOBJS)

CXXFLAGS += -I../FB -I../F_CORE -I$(CPMLIB_DIR)/include
CXXFLAGS += -DFFVC_VERSION=\"$(FFVC_VERSION)\" -DFFVC_MINI_VERSION=\"$(FFVC_MINI_VERSION)\"
LDFLAGS += -L../FB -L../F_CORE -L$(CPMLIB_DIR)/lib
LIBS += -lFB -lFFVcore -lcpmlib

all: $(TARGET)

ifneq (, $(findstring -DPROF_MAPROF, $(CXXFLAGS)))

MAPROF_DIR = ma_prof/src
MAPROF_LIB = $(MAPROF_DIR)/libmaprof.a

CXXFLAGS += -I$(MAPROF_DIR)
LDFLAGS += -L$(MAPROF_DIR)
LIBS += -lmaprof

export
$(MAPROF_LIB):
	$(MAKE) -C $(MAPROF_DIR) cxx_mpi MAPROF_CXX="CXX CXXFLAGS" MAPROF_F="F90 F90FLAGS"

$(TARGET): $(MAPROF_LIB) $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

else

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

endif


.C.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJS)
ifneq (, $(findstring -DPROF_MAPROF, $(CXXFLAGS)))
	$(MAKE) -C $(MAPROF_DIR) clean
endif

test: $(TARGET)
	cd ../../test && ./go.sh

#depend:
#	makedepend -I. -I../IP -I../FB -I../F_CORE -Y $(CXXSRCS)

# DO NOT DELETE

main.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
main.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
main.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
main.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h ../FB/FB_Ffunc.h
main.o: ../FB/dfi.h CommandLine.h CLParser.h timing.h
ffv.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
ffv.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
ffv.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h ../FB/FB_Ffunc.h
ffv.o: ../FB/dfi.h CommandLine.h CLParser.h
ffv_Initialize.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Initialize.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h
ffv_Initialize.o: ../FB/BndOuter.h ../FB/Medium.h ../FB/VoxInfo.h
ffv_Initialize.o: ../FB/Component.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h
ffv_Initialize.o: ../FB/History.h ../FB/FileIO.h ../FB/FB_Ffunc.h ../FB/dfi.h
ffv_Initialize.o: CommandLine.h CLParser.h
ffv_Alloc.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Alloc.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
ffv_Alloc.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
ffv_Alloc.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h
ffv_Alloc.o: ../FB/FB_Ffunc.h ../FB/dfi.h CommandLine.h CLParser.h
ffv_SetBC.o: ffv_SetBC.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_SetBC.o: ../FB/BndOuter.h ../FB/Component.h ffv_Define.h
ffv_SetBC.o: ../F_CORE/ffv_Ffunc.h
ffv_Loop.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Loop.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
ffv_Loop.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
ffv_Loop.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h
ffv_Loop.o: ../FB/FB_Ffunc.h ../FB/dfi.h CommandLine.h CLParser.h timing.h
ffv_LS.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_LS.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
ffv_LS.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
ffv_LS.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h
ffv_LS.o: ../FB/FB_Ffunc.h ../FB/dfi.h CommandLine.h CLParser.h timing.h
ffv_Post.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Post.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h ../FB/BndOuter.h
ffv_Post.o: ../FB/Medium.h ../FB/VoxInfo.h ../FB/Component.h ffv_SetBC.h
ffv_Post.o: ../F_CORE/ffv_Ffunc.h ../FB/History.h ../FB/FileIO.h
ffv_Post.o: ../FB/FB_Ffunc.h ../FB/dfi.h CommandLine.h CLParser.h timing.h
NS_FS_E_Binary.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
NS_FS_E_Binary.o: ffv_Define.h ../FB/FBUtility.h ../FB/ItrCtl.h
NS_FS_E_Binary.o: ../FB/BndOuter.h ../FB/Medium.h ../FB/VoxInfo.h
NS_FS_E_Binary.o: ../FB/Component.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h
NS_FS_E_Binary.o: ../FB/History.h ../FB/FileIO.h ../FB/FB_Ffunc.h ../FB/dfi.h
NS_FS_E_Binary.o: CommandLine.h CLParser.h timing.h
