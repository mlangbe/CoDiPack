#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
# Homepage: http://www.scicomp.uni-kl.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
#
# This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
#
# CoDiPack is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# CoDiPack is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU General Public License for more details.
# You should have received a copy of the GNU
# General Public License along with CoDiPack.
# If not, see <http://www.gnu.org/licenses/>.
#
# Authors:
#  - SciComp, TU Kaiserslautern:
#     Max Sagebaum
#     Tim Albring
#     Johannes Blühdorn
#

# names of the basic directories
TEST_DIR = tests
BUILD_DIR = build
DRIVER_DIR = drivers
RESULT_BASE_DIR = testResults
RESULT_DIR = results

#list all source files in TEST_DIR
TEST_FILES   = $(wildcard $(TEST_DIR)/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/**/Test**.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/Test**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/**/Test**.d)

FLAGS = -Wall -Wextra -pedantic -std=c++11 -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true -DCODI_EnableCombineJacobianArguments

# The default is to run all drives
DRIVERS?=ALL

ifeq ($(OPT), yes)
  CXX_FLAGS := -O3 $(FLAGS)
else
  CXX_FLAGS := -O0 -g $(FLAGS)
endif

ifeq ($(CXX), )
	CXX := g++
else
	CXX := $(CXX)
endif

ifeq ($(UPDATE), yes)
  UPDATE_FLAG=-u
endif

CODI_DIR := ..
CODI_INCLUDE = -I$(CODI_DIR)/include

# Tests will contain the complete list of all tests
TESTS =

# set default rule
all:

# disable the deletion of secondary targets
.SECONDARY:

### define general sets for tests
# Tests for all tapes
BASIC_TESTS = $(wildcard $(TEST_DIR)/basic/Test**.cpp) $(wildcard $(TEST_DIR)/expressions/Test**.cpp) $(wildcard $(TEST_DIR)/exceptions/Test**.cpp) $(wildcard $(TEST_DIR)/jacobiOptimizations/Test**.cpp) $(wildcard $(TEST_DIR)/helpers/Test**.cpp)
BASIC_TESTS_EH = $(filter-out $(TEST_DIR)/helpers/TestStatementPushHelper.cpp, $(BASIC_TESTS))
# Tests that run only for reverse mode tapes
REVERSE_TESTS = $(wildcard $(TEST_DIR)/external_functions/Test**.cpp) $(wildcard $(TEST_DIR)/io/Test**.cpp) $(wildcard $(TEST_DIR)/helpers/reverse/Test**.cpp)
# Tests that run for non vector mode tapes
REVERSE_VALUE_TESTS = $(wildcard $(TEST_DIR)/preaccumulation/Test**.cpp)

# The build rules for all drivers.
define DRIVER_RULE
$(BUILD_DIR)/$(DRIVER_MODE)/$(DRIVER_NAME).o : $(DRIVER_SRC)
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $(DRIVER_INC) $< -c -o $@
	@$(CXX) $(CXX_FLAGS) $(DRIVER_INC) $< -MM -MP -MT $@ -MF $@.d

$(BUILD_DIR)/$(DRIVER_MODE)/%_$(DRIVER_NAME).o : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $(DRIVER_INC) $< -c -o $@
	@$(CXX) $(CXX_FLAGS) $(DRIVER_INC) $< -MM -MP -MT $@ -MF $@.d

$(BUILD_DIR)/$(DRIVER_MODE)/%_$(DRIVER_NAME).exe : $(BUILD_DIR)/$(DRIVER_MODE)/%_$(DRIVER_NAME).o $(BUILD_DIR)/$(DRIVER_MODE)/$(DRIVER_NAME).o
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $(DRIVER_LIB) $^ -o $@

DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/$(DRIVER_MODE)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/$(DRIVER_MODE)/%.test: $(RESULT_DIR)/$(DRIVER_MODE)/%_$(DRIVER_NAME).out
endef

# The general template for each driver. It checks if the driver should be run
define DRIVER_INST
  ifeq (${DRIVER_NAME}, $(filter ${DRIVER_NAME}, ${DRIVERS}))
    $(eval $(value DRIVER_RULE))
  else
    ifeq (ALL, $(filter ALL, ${DRIVERS}))
      $(eval $(value DRIVER_RULE))
    endif
  endif
endef

include Makefile.primal
TESTS += $(TESTS_PRIMAL)

include Makefile.deriv1stOrder
TESTS += $(TESTS_1ST_ORDER)

include Makefile.deriv2ndOrder
TESTS += $(TESTS_2ND_ORDER)

# rules for generating the test files
$(RESULT_DIR)/%.out : $(BUILD_DIR)/%.exe
	@mkdir -p $(@D)
	$< > $@

$(BUILD_DIR)/%.exe : %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@
	@$(CXX) $(CXX_FLAGS) $< -MM -MP -MT $@ -MF $@.d

# rule for printing the results (dependencies are generated by the drivers)
$(RESULT_DIR)/%.test: $(BUILD_DIR)/compare.exe
	@sh ./compare.sh $(UPDATE_FLAG) -n $* -b $(RESULT_BASE_DIR)/$*.out $(filter-out $(BUILD_DIR)/compare.exe,$^)

all: $(TESTS)
	@mkdir -p $(BUILD_DIR)

all0th: $(TESTS_PRIMAL)
	@mkdir -p $(BUILD_DIR)

all1st: $(TESTS_1ST_ORDER)
	@mkdir -p $(BUILD_DIR)

all2nd: $(TESTS_2ND_ORDER)
	@mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)
	rm -fr $(RESULT_DIR)

-include $(DEP_FILES)
