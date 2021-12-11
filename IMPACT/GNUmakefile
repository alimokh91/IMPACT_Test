.PHONY: clean tests all
.DEFAULT_GOAL = all
.PRECIOUS: %.F90

TOP_DIR := $(shell pwd)
SRC_DIR=$(TOP_DIR)/src
DST_DIR=$(TOP_DIR)/prog
TEST_DIR=$(TOP_DIR)/tests

VPATH = . $(SRC_DIR) $(TEST_DIR)

include $(PFUNIT)/include/base.mk

ifeq ($(USEMPI),YES)
   FC=mpifort
endif

EXE_TEST = tests$(EXE_EXT)

all: $(EXE_TEST)

SUT: IMPACT_TESTING=1
SUT:
	make -C $(SRC_DIR)
	make -C $(TEST_DIR) tests

tests: EXE_BUILD=1
tests:
	make -C $(TEST_DIR) tests
	make -C $(SRC_DIR)

ifeq ($(USEMPI),YES)
	mpirun -np 1 ./$(EXE_TEST)
else
	./$(EXE_TEST)
endif

$(EXE_TEST): testSuites.inc SUT
#	$(FC) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(DST_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(FFLAGS) $(FPPFLAGS)

clean: local-E0-clean

local-E0-clean:
	$(MAKE) -C $(SRC_DIR) clean
	$(MAKE) -C $(TEST_DIR) clean
	$(RM) -f $(EXE_TEST) *$(OBJ_EXT)

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export DST_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
export IMPACT_TESTING
export EXE_TEST
export EXE_BUILD
