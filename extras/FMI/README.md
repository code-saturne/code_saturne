What is FMI2 ?
==============

The [Functional Mock-up Interface standard](fmi-standard.org) (FMI) provides an interface
standard for coupling of simulation tools in a co-simulation environment. A
master algorithm controls the data exchange and synchronization between all
simulation solvers (slaves).

The current implementation for code_saturne implements version 2.0 of the FMI
specification.

What is a FMU ?
===============

A co-simulation slave is called a __FMU__ (Functional Mock-up Unit) and is
distributed in one ZIP file which contains several files:
* An XML file defining every variables in the FMU and other static information
(path for example);
* All required model equations or the access to co-simulation tools are provided
with a small set of C or C++ functions.

Client-server code_saturne FMU model
====================================

The FMU behaves as a client, controlling a standard code_saturne run as a server.

When the FMU starts, it executes the following steps:
* Call `code_saturne run --initialize` in the specified
  case directory, to prepare the run.
* Initialize a server socket.
* Generate a `control_file` in the case's execution directory, specifying
  the associated FMU host and port on which the communication will occur,
  along with a randomly-generated key (only accessible to the user running the
  computation, as others should not have read access to that file).
* Start the server in the specified execution directory, using the
  `run_solver` script (generated at the first step).
  -  When starting, the solver will read that `control_file` and attempt
     a connection with the FMU client.
* Accept the handshake from this server using the provided key.
* Send the lists of variables to be sent and received by the FMU;
  - All of these variables (defined in the FMU model) must have a matching
    notebook entry in the case setup.
  - The setup may optionally include additional standard notebook entries
    not used or controlled by the FMU.

After these steps, the computation starts, exchanging the FMI model variables with the
FMU at each time step.

The number of time steps or physical time that should be executed, as defined in
the case setup, is ignored when connected to the FMU. The FMU sends a stop message
to the server when the FMI environment instructs it to halt.

Building a code_saturne FMU
===========================

To build an FMU, the `modelDescription.xml` file can be copied from
the `example_xml` directory to the current directory, and adapted
to the desired model.

Note that in the model description, the following string variables must always
be defined:
- _code_saturne_ specifies the path of the main `code_saturne` script
  (from a valid code_saturne installation).
- _casename_ specifies the path of the case that will be controlled by the
  FMU client.
- _run_id_ specifies the matching run id (so that the actual execution
  directory will be: `<casename>/RESU/<run_id>`).

The following steps are then required:
* run `./fmu_generator.py`, which will generate a `code_saturne_fmu_variables.h`
  file from the XML model description.
* run `Make`, which wil generate the FMU.

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

Optional helper tools
=====================

Various third-party tools, including open-source ones, can be used to assist
in the generation of an FMU.

For example, the FMU Builder tool developed by Monentia and available at
https://bitbucket.org/simulage/c-fmu-builder/wiki/Home
can be used to generate an FMU.

When using such a tool, only the resulting `modelDescription.xml`. file is of
interest to us, as the skeleton FMU would generally require further adaptation,
but  so it is recommended
to run this step in a separate directory.

The generated model description must contain at least the _code_saturne_,
_casename_, and _run_id_ ScalarVariable definitions, as in the example file.

Logging and debugging
=====================

Several environment variables allow logging the FMU's behavior.

On the FMU side:
- `CS_FMU_COMM_TRACE=<file_name>` allows defining a file in which all
   communications with the connected code_saturne instance are logged.
- `CS_FMU_TRACE=1` allows logging all-low level FMU calls and actions
  to the standard output, except those who are already logged using the
  FMI environment.

On the code_saturne side:
- `CS_CONTROL_COMM_TRACE=<file_name>` allows defining a file in which all
   communications with the connected code_saturne FMU instance are logged
   (which allows checking that calls match those logged by the FMU using
  `CS_FMU_COMM_TRACE`).
- `CS_CONTROL_RECV_LOG=<file_name>` allows defining a file in which commands
   received by the FMU are logged. This allows generating a `control_file`
   which can be adapted or used to replay a computation in standalone mode.

