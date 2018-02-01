# Example files for IRC4Orca

In this directory, you'll find some usage examples for IRC4Orca:

1. Proton transfer in AcOH
2. Aza-Diels-Alder reaction

These examples are known to work with the current version of IRC4Orca (using python 3.6.x) and Orca-4.0.1.

## Proton Transfer in AcOH
The files provided for this example allow the optimization of the Transition State (TS) at the TPSS/Def2-TZVP level of theory, and the subsequent input files for calculating the IRC in the foward and bakwards direction. Also provided are the expected results from both IRC4Orca runs.
* __acoh-migration-ts-optts.inp__ - input file for Orca (TS search and optimization)
* __acoh-migration-ts-optts.hess__ - hess file from Orca needed for IRC
* __acoh-migration-ts-optts.gbw__ - gbw file from Orca 4.0.1 used as initial guess in the IRC calculations (optional, but highly recomended)
* __acoh-migration-ts-irc-f.inp__ - Input file for IRC4Orca: IRC in the forward direction
* __acoh-migration-ts-irc-r.inp__ - Input file for IRC4Orca: IRC in the reverse direction
* __acoh-migration-ts-irc-f_example.log__ - Example output from IRC4Orca: IRC in the forward direction
* __acoh-migration-ts-irc-f_example.trj__ - Example trajectory from IRC4Orca: IRC in the forward direction
* __acoh-migration-ts-irc-r_example.log__ - Example output from IRC4Orca: IRC in the backward direction
* __acoh-migration-ts-irc-r_example.trj__ - Example trajectory from IRC4Orca: IRC in the backward direction
* __acoh-migration-ts-irc-merged.trj__ - Merged trajectories from IRC4Orca using irc-concatenate (_see_ utils directory).

## Aza-Diels-Alder Reaction
This example was adapted from _RSC Adv._, __2015__, 5, 50729-50740 and refers to the rate-determining step of the 2+4 cycloaddition of a protonated imine with cyclopentadiene. This example is substantially more demanding on the computational resourses that the previous one. Still, the optimization of the TS and the IRC calculations should be done in about one afternoon using a modest Core i7 processor.

This example illustrates the usage of IRC4Orca when dealing with a TS located on a relatively flat region of the Potential Energy Surface (PES), specially when moving in the backwards direction. To overcome this, a larger #ircmaxd was used in both calculations. Additionally, we increase #ircdamp to 0.1 in order to mix 10% of the previous displacement in the initial guess displacement of the current step.

The following files are supplied (adopting the same nomenclature as the previous example):
* __ada-ts-optts.inp__ - input file for Orca (TS search and optimization)
* __ada-ts-optts.hess__ - hess file from Orca needed for IRC
* __ada-ts-optts.gbw__ - gbw file from Orca 4.0.1 used as initial guess in the IRC calculations (optional, but highly recomended)
* __ada-ts-irc-f.inp__ - Input file for IRC4Orca: IRC in the forward direction
* __ada-ts-irc-r.inp__ - Input file for IRC4Orca: IRC in the reverse direction
* __ada-ts-irc-f_example.log__ - Example output from IRC4Orca: IRC in the forward direction
* __ada-ts-irc-f_example.trj__ - Example trajectory from IRC4Orca: IRC in the forward direction
* __ada-ts-irc-r_example.log__ - Example output from IRC4Orca: IRC in the backward direction
* __ada-ts-irc-r_example.trj__ - Example trajectory from IRC4Orca: IRC in the backward direction
* __ada-ts-irc-merged.trj__ - Merged trajectories from IRC4Orca using irc-concatenate (_see_ utils directory).


