#MAKE_SETTING ?= ../../make_setting

#include $(MAKE_SETTING)

TARGET = ../lib/libcpmlib.a

AR ?= ar
ARFLAGS ?= rv
RANLIB ?= ranlib
 
SRCS = cpm_DomainInfo.cpp \
       cpm_ParaManager_frtIF.cpp \
       cpm_ParaManager_Alloc.cpp \
       cpm_ParaManager_MPI.cpp \
       cpm_VoxelInfo.cpp \
       cpm_ParaManager.cpp

FSRCS =

.SUFFIXES: .o .cpp

OBJS=$(SRCS:.cpp=.o)


$(TARGET):$(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../include -o $@ -c $< 

clean:
	$(RM) *.o $(TARGET)
