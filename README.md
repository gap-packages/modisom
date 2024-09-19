[![Build Status](https://github.com/gap-packages/modisom/workflows/CI/badge.svg?branch=master)](https://github.com/gap-packages/modisom/actions?query=workflow%3ACI+branch%3Amaster)
[![Code Coverage](https://codecov.io/github/gap-packages/modisom/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/modisom)


# ModIsom

ModIsom is a GAP4 package. It contains various methods for computing with 
nilpotent associative algebras. In particular, it contains a method to 
determine the automorphism group and to test isomorphis of such algebras 
over finite fields and of modular group algebras of finite p-groups, and 
it contains a nilpotent quotient algorithm for finitely presented associative 
algebras and a method to determine Kurosh algebras.
Some of the functions to compute with nilpotent algebras are tailored towards
application for the Modular Isomorphism Problem.
Moreover, the package allows to compute various group theoretical properties of groups
related to the Modular Isomorphism Problem.


## Installation

To install the package ModIsom move the file `modisom-X.Y.Z.tar.gz`
(or any other archive containing it) into the `pkg` directory.
Usually this will be the `pkg` subdirectory in your GAP4 installation.
However, it is also possible to have a `pkg` subdirectory in a different
place, see the section `Installing GAP Packages` of the GAP4 reference
manual for more information. Then simply unpack `modisom-XX.tar.gz` and
your installation is complete. To load the package, call

    gap> LoadPackage( "modisom" );

             
## Documentation

The manual of ModIsom is contained in `modisom/doc` and `modisom/htm`
directories.
# modisom3
