MAKE_SETTING ?= make_setting
include $(MAKE_SETTING)
export 
CPMLIB_DIR = CPMlib

.PHONY: CPMlib FB FFV F_CORE

all: FFV

FFV: CPMlib FB F_CORE

#BASE_DIR = $(PWD)
#export MAKE_SETTING := $(BASE_DIR)/$(MAKE_SETTING)

CPMlib:
	$(MAKE) -C $(CPMLIB_DIR)/src

FB:
	$(MAKE) -C FB

F_CORE:
	$(MAKE) -C F_CORE

FFV:
	$(MAKE) -C FFV

clean:
	$(MAKE) -C $(CPMLIB_DIR)/src clean
	$(MAKE) -C FB clean
	$(MAKE) -C F_CORE clean
	$(MAKE) -C FFV clean

test: FFV
	cd ../test && ./go.sh
