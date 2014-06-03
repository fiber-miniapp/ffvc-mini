# Copyright (C) 2014 RIKEN AICS
# This library is released under the terms of the MIT license.
# http://fiber-miniapp.mit-license.org/

MAX_SECTIONS ?= 20
OUTPUT ?= maprof_output.yaml
#OUTPUT ?= stdout

PFLAGS :=

ifdef MAPROF_F
PFLAGS += -DMAPROF_FC='"$($(word 1,$(MAPROF_F)))"' \
          -DMAPROF_FC_VERSION='"$(shell ./check_version.sh $($(word 1,$(MAPROF_F))))"' \
		  -DMAPROF_FFLAGS='"$(subst \\,\,$(subst ",\",$($(word 2,$(MAPROF_F)))))"'
endif

ifdef MAPROF_C
PFLAGS += -DMAPROF_CC='"$($(word 1,$(MAPROF_C)))"' \
          -DMAPROF_CC_VERSION='"$(shell ./check_version.sh $($(word 1,$(MAPROF_C))))"' \
		  -DMAPROF_CFLAGS='"$(subst \\,\,$(subst ",\",$($(word 2,$(MAPROF_C)))))"'
endif

ifdef MAPROF_CXX
PFLAGS += -DMAPROF_CXX='"$($(word 1,$(MAPROF_CXX)))"' \
          -DMAPROF_CXX_VERSION='"$(shell ./check_version.sh $($(word 1,$(MAPROF_CXX))))"' \
		  -DMAPROF_CXXFLAGS='"$(subst \\,\,$(subst ",\",$($(word 2,$(MAPROF_CXX)))))"'
ifndef MAPROF_C
CC = $($(word 1,$(MAPROF_CXX)))
CFLAGS = $($(word 2,$(MAPROF_CXX)))
endif
endif

CFLAGS += -DMAPROF_MAX_SECTIONS=$(MAX_SECTIONS) -DMAPROF_OUTPUT='"$(OUTPUT)"'

CFLAGS += $(PFLAGS)

C_OBJS = maprof.o maprof_time.o maprof_yaml.o maprof_proc.o
F_OBJS = mod_maprof.o

all: c

c_mpi cxx_mpi f_mpi: CFLAGS += -DUSE_MPI
f_mpi: FFLAGS += -DUSE_MPI

c c_mpi cxx cxx_mpi: libmaprof.a
f f_mpi: libmaprof_f.a

AR ?= ar
ARFLAGS ?= rv
RANLIB ?= ranlib

libmaprof.a: $(C_OBJS)
	$(AR) $(ARFLAGS) $@ $(C_OBJS)
	$(RANLIB) $@

libmaprof_f.a: $(C_OBJS) $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $(C_OBJS) $(F_OBJS)
	$(RANLIB) $@

a.out: main.o libmaprof.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

proc_test: proc_test.o maprof_proc.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

.SUFFIXES:
.SUFFIXES: .F90 .c .o .mod

.c.o:
	$(CC) $(CFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod *.i90 *.a proc_test


maprof_yaml.o: maprof_yaml.h
maprof_proc.o: maprof_proc.h
maprof.o: maprof_yaml.h maprof_proc.h
maprof.o: maprof.h
maprof_time.o: maprof.h
main.o: maprof.h
proc_test.o: maprof_proc.h

