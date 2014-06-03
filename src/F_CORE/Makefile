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

CPMLIB_DIR = ../CPMlib

TARGET = libFFVcore.a

AR ?= ar
ARFLAGS ?= rv
RANLIB ?= ranlib

F90SRCS = \
  ffv_utility.f90 \
  ffv_vbc_outer.f90 \
  ffv_velocity_binary.f90 \
  ffv_poisson.f90 

SRCS = $(F90SRCS)

.SUFFIXES: .o .f90
OBJS = $(F90SRCS:.f90=.o)

$(TARGET):$(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.f90.o:
	$(F90) $(F90FLAGS) -I. -I../FB -I$(CPMLIB_DIR)/include -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

#depend:
#	makedepend -I. -I../FB -I../FFV -Y $(SRCS)

# DO NOT DELETE

ffv_poisson.o: ../FB/ffv_f_params.h
ffv_utility.o: ../FB/ffv_f_params.h
ffv_vbc_outer.o: ../FB/ffv_f_params.h
ffv_velocity_binary.o: ../FB/ffv_f_params.h load_var_stencil5.h muscl.h d_o_o_p.h
