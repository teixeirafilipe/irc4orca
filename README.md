# IRC4Orca
## An Implementation of Morokuma's IRC method for the Orca Electronic Structure Software package

IRC4Orca is an implementation of Morokuma's IRC method (adapted from *J. Chem.  Phys.*, __1977__, 66, 2153-2156) for the Orca Electronic Structure Software package. This is a wrapper script written in Python that calls Orca to provide the necessary Single-Point Energy and Gradient calculations. Because the core of the IRC routines lies separate from the Orca interface, this program can be easily adapted to provide IRC capabilities to other ESS packages, as needed.

The need for an IRC module for Orca emerged at the end of my PhD work, in 2013, and was first applied in the calculation of the IRC for vanadium catalysed epoxidations (although later replaced by _Ab Initio_ Molecular Dynamics calculations in _Inorg. Chem._, __2017__, 56, 2124-2134). The original python script used python 2.7 and was designed to work with Orca 2.8. 

Due to some interest in this script capabilities, I decided to maintain the code in late 2017. This required adapting the interface to Orca as well as updating the code to be able to run in the more recent dialect of Python 3. In its current state (as of January 2018) the code is capable of performing Mororuma's IRC calculations in Cartesian coordinates using Orca 4.0.1.

## Requirements and Usage
IRC4Orca does not require much from the system, as it only performs the manipulations on a system's geometry. On the other hand, _Ab Initio_ calculations might be quite demanding. Because each step along the IRC trajectory requires at least 3 energy calculations and one gradient calculation, the process of finding the IRC can be quite lengthy, specially for large systems and basis sets.

The current version of IRC4Orca (as of January 2018) was designed to work with __Python__ 3.6.x and should work with any of the 3.x and 4.0.x versions of Orca. It also requires __Numpy__ (tested on numpy-1.14.0).

IRC4Orca requires one input file as its argument, the contents of which are explained in the next section. Its usage (assuming irc4orca.py is in your $PATH) is:

```
$ irc4orca.py file.inp
```

At the end of the calculations, you should be presented with two additional files: file.log and file.trj. The former is a human-readable text file containing the results from the calculations, while the latter is a trajectory file, in `xyz` format which can be used to visualise the trajectory using Molden.  

## Input File Structure
The input file is a normal Orca input file, although I recommend it not to bear any indication as to the type of calculation (Opt, enGrad, SP or MD keywords). The parameters to be read by IRC4Orca are given as comments, with one instruction per line:

### Mandatory instructions
* __#orcacmd *cmdname*__ -  location of an executable file that  takes an Orca input as argument and directs Orca's output to a filen with the same base name and the.out extension
* __#irchess *filename.hess*__ - location of the file hess file for the system
* __#ircmode *n*__ - Vibrational mode to follow (default:0, should be set to 6 in most cases)

### Optional Instructions
* __#ircrestart *[0/1]*__ - This is a restart from a point in the IRC that is not the TS (default:0/False) 
* __#ircguess *filename.gbw*__ - location of a GBW for the first point    
* __#ircalpha *x.xx*__ - Alpha parameter for scaling the initial displacement (default: 0.1).             
* __#ircdelta *x.xx*__ - Scaling parameter for finding the closest minimum at a given point (default: 0.05).
* __#ircdamp *x.xx*__ - Percentage of the previous displacement on the calculation of the current one    (default: 0.05).                         
* __#ircautodamp *[0/1]*__ - Update ircdamp. The update procedure reduces ircdump by 0.1/log(n) for n>1 (default: 0, False)                      
* __#irctol *x.xx*__ - Acceptable value of the RMS gradient for considering the IRC to have converged (default:1.0e-4).                        
* __#ircdir *[-1/+1]*__ - Direction of the initial displacement (default: +1).
* __#ircmaxd *x.xx*__ - Norm of the displacement vector between IRC poinnts (default: 0.01 angs/amu^(1/2)).          
* __#ircpts *n*__ - Maximum number of points in the IRC (default: 25)
* __#ircalg *[1/2]*__ - IRC algorithm: 1 (default) is the traditional Morokuma algorithm. 2 is an updated version that requires 3 additional energy evaluations per step.

## Known Issues 
* The Orca interface was only tested for DFT and HF calculations.  Because of the different keywords used to navigate Orca's output, the script cannot understand the output generated from semi-empirical methods (AM1, PM3, etc). It remains untested for IRC calculations using post-HF methods, although it should work fine with MP2.

* The defaults are a little bit conservative. IRC calculations on large systems (or in cases where the atomic motions associated with the TS are spread over a large number of atoms) might need a larger #ircalpha or #ircmaxd in order to work.

 * On very flat regions of the Potential Energy Surface, the IRC tracking can "de-rail", using a larger #ircdamp parameter usually ameliorates this issue. Because of this the IRC calculation might end with an increase in energy while still considerably far from the end-point.  A geometry optimisation starting from the last IRC geometry is thus highly recommended.

## Citation
The author highly recomends the citation of Morokuma's original paper (*J. Chem. Phys.*, __1977__, 66, 2153-2156) as well as any relevant papers describing the level of theory and implemntation decisions involved in the underlying Orca calculations (please refer to the Orca Manual for that purpose). As for the specific implementation of IRC4Orca, it may be cited using the following BibTeX entry:

```bibtex
@Electronic{TeixeiraIRC4Orca,
  author    = {Filipe Teixeira},
  title     = {IRC4Orca - An Implementation of Morokuma's IRC method for the Orca Electronic Structure Software package.},
  year      = {2018},
  date      = {2018-01-30},
  url       = {https://github.com/teixeirafilipe/irc4orca},
  urldate   = {2018-01-30}
}
```
