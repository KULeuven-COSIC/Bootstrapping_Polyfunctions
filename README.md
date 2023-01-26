# Bootstrapping_Polyfunctions

Improved digit extraction procedure for BGV and BFV bootstrapping based on polyfunctions theory. The implementation relies on the HElib library that is included in this repository. Note that some of the header and source files of HElib have been adapted, so those are not identical to the ones found in https://github.com/homenc/HElib. Reinstallation of HElib is therefore required.

Any bugs can be reported to robin.geelen@esat.kuleuven.be.

## Installation

* Compile HElib via HElib/build_install (it is assumed that GMP and NTL are already installed in /usr/local).
* Compile and run HElib/Polynomials/main.cpp which performs bootstrapping for toy parameters. One can also run the commented parameter sets in that file. For running with different parameters, the appropriate polynomials need to be generated first (see below).

## Generate digit extraction polynomials

* Obtain a Magma license from http://magma.maths.usyd.edu.au/magma.
* Start Magma and change to the root directory of this repository via the “ChangeDirectory” command.
* Change the parameters in Scripts/Find_digit_poly.m and load it via the “load” command. This will generate the required polynomials and write them to a file.