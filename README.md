# ConvectionDiffusionWrapper

To compile the shared library open Terminal and run (currently not possible):

    cd $HOME/ug4/plugins/ConvectionDiffusion/
    git clone https://github.com/pzugel/Convection_Diffusion_Simulation
    git pull
    cd wrapper/build
    cmake --build .
    
Alternatively, if cmake fails (USE THIS FOR NOW):

    gcc -DUG_ALGEBRA -DUG_BRIDGE -DUG_CPU_1 -DUG_DEBUG -DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3 -DUG_DISC -DUG_FOR_LUA -DUG_GRID -DUG_PLUGINS -DUG_POSIX -DUG_PROFILER_SHINY_CHECK_CONSISTENCY -DUG_TARGET=\"ugshell\" -I/home/paul/ug4/externals/BoostForUG4 -I/home/paul/ug4/ugcore/ugbase -L/home/paul/ug4/lib -L/home/paul/ug4/bin/plugins -Wl,-rpath=/home/paul/ug4/lib -Wl,-rpath=/home/paul/ug4/bin/plugins -shared -o libConvectionDiffusionWrapper.so -fPIC -w convection_diffusion_wrapper.cpp -lug4 -lConvectionDiffusion


The compiled library will be located in _$HOME/ug4/bin/plugins_. Copy the _libConvectionDiffusionWrapper.so_ into _LabView/lib_ or use the pre-compiled library that is already present. Please note that this library will only work if your local version uf UG4 was compiled with the above compiler directives.

TODO: 
* Wrapper only chooses between FV1 and FE
* Implement SubsetHandler  
* Print out UG_LOG
* Generate Plot via VTU Output
* Neumann Boundary
* Choose Subsets Inner/Boundary
