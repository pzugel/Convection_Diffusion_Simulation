# Convection_Diffusion_Simulation

This project contains a **LabVIEW** application to simulation convection/diffusion flows via **UG4**. 

## Setting up UG

To get stated first install UG via the [ughub](https://github.com/UG4/ughub) pakage manager. Follow the instructions, make sure to set up UG in your home directory as proposed, e.g.

```
$HOME/ug4 (Linux) 
%HOMEPATH%\ug4 (Windows)
```
We refer to this path as `$UGPATH`.

Don't forget to install the ConvectionDiffusion Plugin: 
```
cd $UGPATH
ughub install ConvectionDiffusion
```

When compiling UG please use the following options:
```
cd $UGPATH/build
cmake -DDIM=ALL -DCPU=1 -DTARGET=ugshell -DPOSIX=ON -DConvectionDiffusion=ON .
```

## Setting up LabVIEW

When you finished compiling UG clone this repository and build the required libraries for LabVIEW:

```
cd $UGPATH/plugins/ConvectionDiffusion/
mkdir labview_wrapper
cd labview_wrapper
git clone https://github.com/pzugel/Convection_Diffusion_Simulation .
make
```

This should create a library object (as **.dll** or **.so**) in your `labview_wrapper/LabView/lib` folder. Once the library is build, feel free to move the `labview_wrapper` folder to any location on your system and/or rename it. 

To start the application you need to install LabVIEW, which you can download from the *National Instruments* [web page](https://www.ni.com/de-de/support/downloads/software-products/download.labview.html). Make sure select at least Version 2019.

Finally open up the LabVIEW project file:

> Convection_diffusion.lvproj

and run the *main.vi*.



## TODO: 
* Wrapper only chooses between FV1 and FE
* Implement SubsetHandler
* Print out UG_LOG
* Generate Plot via VTU Output
* Neumann Boundary
* Choose Subsets Inner/Boundary
