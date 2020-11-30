# ConvectionDiffusionWrapper

To compile the shared library open Terminal and run (MISSING COMPILER FLAGS!):
```
cd /home/paul/ug4/plugins/ConvectionDiffusion/wrapper/build
cmake --build .
```
Alternatively (USE THIS FOR NOW):
```
gcc -DUG_ALGEBRA -DUG_BRIDGE -DUG_CPU_1 -DUG_DEBUG -DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3 -DUG_DISC -DUG_FOR_LUA -DUG_GRID -DUG_PLUGINS -DUG_POSIX -DUG_PROFILER_SHINY_CHECK_CONSISTENCY -DUG_TARGET=\"ugshell\" -I/home/paul/ug4/externals/BoostForUG4 -I/home/paul/ug4/ugcore/ugbase -L/home/paul/ug4/lib -L/home/paul/ug4/bin/plugins -Wl,-rpath=/home/paul/ug4/lib -Wl,-rpath=/home/paul/ug4/bin/plugins -shared -o libConvectionDiffusionWrapper.so -fPIC -w convection_diffusion_wrapper.cpp -lug4 -lConvectionDiffusion
```

The compiled library is located in /home/paul/ug4/bin/plugins.

The main.cpp runs the wrapper with some predefined values in a 2D Domain. To compile, run:
```
g++ -DUG_ALGEBRA -DUG_BRIDGE -DUG_CPU_1 -DUG_DEBUG -DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3 -DUG_DISC -DUG_FOR_LUA -DUG_GRID -DUG_PLUGINS -DUG_POSIX -DUG_PROFILER_SHINY_CHECK_CONSISTENCY -DUG_TARGET=\"ugshell\" -I/home/paul/ug4/externals/BoostForUG4 -I/home/paul/ug4/ugcore/ugbase -L/home/paul/ug4/lib -L/home/paul/ug4/bin/plugins -Wl,-rpath=/home/paul/ug4/lib -Wl,-rpath=/home/paul/ug4/bin/plugins -w main.cpp -lug4 -lConvectionDiffusion
```

TODO: 
* Wrapper only chooses between FV1 and FE
* Implement SubsetHandler  
* Print out UG_LOG
* Generate Plot via VTU Output
* Neumann Boundary
* Choose Subsets Inner/Boundary
