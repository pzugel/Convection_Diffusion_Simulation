#=============================================================================
# Define directories

UG_PATH = ${HOME}/ug4
PLUGIN_PATH = $(UG_PATH)/plugins/ConvectionDiffusion
LIB_PATH = $(PLUGIN_PATH)/labview_wrapper/LabView/lib
SOURCE_PATH = $(PLUGIN_PATH)/labview_wrapper/src

#=============================================================================

all:	
	@echo "Writing shared library ... "
	gcc -DUG_ALGEBRA -DUG_BRIDGE -DUG_CPU_1 -DUG_DEBUG -DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3 -DUG_DISC -DUG_FOR_LUA -DUG_GRID -DUG_PLUGINS -DUG_POSIX -DUG_PROFILER_SHINY_CHECK_CONSISTENCY -DUG_TARGET=\"ugshell\" -I$(UG_PATH)/externals/BoostForUG4 -I$(UG_PATH)/ugcore/ugbase -L$(UG_PATH)/lib -L$(UG_PATH)/bin/plugins -Wl,-rpath=$(UG_PATH)/lib -Wl,-rpath=$(UG_PATH)/bin/plugins -shared -o $(LIB_PATH)/libConvectionDiffusionWrapper.so -fPIC -w $(SOURCE_PATH)/convection_diffusion_wrapper.cpp -lug4 -lConvectionDiffusion
