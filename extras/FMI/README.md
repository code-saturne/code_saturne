What is FMI2 ?
==============

The __Functional Mock-up Interface standard__ (FMI2) provides an interface
standard for coupling of simulation tools in a co-simulation environment. A
master algorithm controls the data exchange and synchronization between all
simulation solvers (slaves).

What is a FMU ?
===============

A co-simulation slave is called a __FMU__ (Functional Mock-up Unit) and is
distributed in one ZIP file which contains several files:
* An XML file defining every variables in the FMU and other static information
(path for example);
* All required model equations or the access to co-simulation tools are provided
with a small set of C or C++ functions.

How to use code_saturne as a FMU
================================

Building a code_saturne FMU
---------------------------

A FMU Builder tool developed by Monentia and available at
https://bitbucket.org/simulage/c-fmu-builder/wiki/Home
can be used to generate an FMU based on code_saturne.

The main steps are the following:
* Open the FMU Builder;
* Define every needed variables and static informations and save
them in a .fpd file;
* Generate the project, creating a directory containing the XML file,
auto-generated C++ files that should not be modified, and a C++ file named
after the project (project.cpp) which much be filled;
* Once every methods are completed (mainly init, doStep, terminate, getters and
setters of every used variable types), the FMU can be built;

Note that when updating the FMU (for example adding variables), the
generated `code_saturne.h` file in the project should be removed before
re-building the FMU, as it will not be overwritten by default.

Using code_saturne as a FMU
---------------------------

When the code_saturne FMU is built, a .fmu file is created. In order to test it,
a few test environments are available. Most tests have been made using a free
python library to simulate FMUs called
[FMPy](https://github.com/CATIA-Systems/FMPy).

The .fmu file can be opened in FMPy and variables can be plotted during the
simulation.

Example
=======

An example is provided in the current directory. This example is based on a
code_saturne simulation on a stratified case (such as the stratified_junction
tutorial) where the cold inlet velocity is set and outlet temperatures are
retrievied using the FMU.

A .fpd file as well as the main project C++ file allow one to build a full
code_saturne FMU. In order to use it, one can follow the steps defined in the 
[Building a code_saturne FMU](#building-a-code_saturne-fmu) section. Paths in
the .fpd file have to be modified.
