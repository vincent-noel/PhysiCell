# Builds a small driver that writes full saves with the options set in a
# config file. Run the checks with: python3 test_full_save_options.py

PROGRAM_NAME := full_save_options

CC := g++
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

CFLAGS := -O2 -fopenmp -m64 -std=c++11
DEPFLAGS := -MMD -MP

DIR := ../..
# same object list as the main Makefile
BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o pugixml.o
PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o \
PhysiCell_signal_behavior.o PhysiCell_rules.o
PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o
OBJECTS := $(addprefix obj/,$(BioFVM_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS))

vpath %.cpp $(DIR)/BioFVM $(DIR)/core $(DIR)/modules

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): main.cpp $(OBJECTS)
	$(CC) $(CFLAGS) $(DEPFLAGS) -o $@ $(OBJECTS) main.cpp

obj/%.o: %.cpp | obj
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm -rf obj $(PROGRAM_NAME) *.d

-include $(OBJECTS:.o=.d)
