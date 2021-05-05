# ConvectionDiffusionWrapper

To compile the shared library open Terminal and run:

    cd $HOME/ug4/plugins/ConvectionDiffusion/
	mkdir labview_wrapper
	cd labview_wrapper
    git clone https://github.com/pzugel/Convection_Diffusion_Simulation
    git pull
    cd labview_wrapper
    make

TODO: 
* Wrapper only chooses between FV1 and FE
* Implement SubsetHandler  
* Print out UG_LOG
* Generate Plot via VTU Output
* Neumann Boundary
* Choose Subsets Inner/Boundary
