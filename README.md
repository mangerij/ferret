### We have moved from BitBucket! ###

As of July 23rd, 2019 the BitBucket repository is no longer supported. All official changes will be made here on GitHub and Ferret business will continue as usual.
The continuous integration, verification, and testing server (CIVET: https://civet.inl.gov/) ensures this version is current with the MOOSE master and that Ferret tests are continuously validated.

### New Features (January 2021) ###

Ferret can treat both ferroelectric phase field and micromagnetic modeling. Both sectors of phenomenological modeling can be implemented separately or coupled together. The new micromagnetic module is based on a Landau-Lifshitz-Bloch approach which can handle finite temperature simulations. Its details will be included in a forthcoming methods paper.

### Access and Policy on Use ###

Ferret still exists under the open-source GNU License but to access the repository you will need an RSA deploy key registered with our application. We freely provide these keys in order to track repository usage and measure impact. We can allow write-access on request.
Please contact one of the developers for assistance in generating one.

After obtaining a key, please cite Mangeri et al Nanoscale, 2017, 9, 1616-1624 in the *main text* of your paper if you use Ferret. A methods paper is in the works which will supercede this requirement in the future. Finally, one more comment about usage: any issues with legacy inputs not running on this version of Ferret should be resolved with corresponding authors of appropriate papers. The same goes for unofficial Ferret patches that assume old versions of MOOSE. As this code belongs under GNU licensing, we cannot possibly support all of the custom spin-offs and additions to the code. 

Before publishing results with Ferret, we strongly suggest that you update MOOSE and then update Ferret. Make a test in the /tests directory representative of your problem and submit a pull-request. This ensures your results will be reproducible on all current and future versions of the entire software stack (Ferret->MOOSE->libMesh->PETSc). You can link to these tests in your papers if you would like to support open science. A final note: assistance from the developers in testing, debugging, upgrading, and using Ferret-related code may constitute a reasonable assumption of co-authorship of any of the Ferret users/developers *and possibly their contributors*. While not strictly required, we still suggest that you clear this up with developers before you submit your manuscript. If you wish to have your paper appear/highlight/quoted on the Ferret website, please send us it after it is published and we will do our best to give you a much deserved spotlight.

### What is Ferret and how can I compile it? ###

Ferret is an application within the finite element [MOOSE](http://mooseframework.org) framework useful for simulating the ferroic nanostructure. More information about features is available at our website (under-construction). To build and use Ferret you will need to build MOOSE. It might be a good idea to learn a few things about MOOSE, too.

Generally, once you have a RSA-deploy key registered with us, you first need to authenticate to Ferret by typing
```
ssh -T git@github.com
```
Once you have access, you will be able to clone the repository with 
```
git clone git@github.com:mangerij/ferret.git
```
Those with admin access to the code can clone with a slightly different command
```
git clone https://github.com/mangerij/ferret.git
```
Make sure Ferret is cloned in the adjacent to the MOOSE folder. For example if MOOSE is in

```
~/projects/moose
```

then Ferret should be

```
~/projects/ferret
```

Ferret can be built with

```
cd <ferret>
./configure
make
```


### Additional information: ###

The ./configure step can include additional libraries such as BOOST which will allow compile of MOOSE objects with mathematical special functions.
Or it can include an advanced fast-multipole boundary element method (FMM-BEM) developed by Prof. Xikai Jiang and co-workers. In order to use the BEM, you need to compile ScalFMM with the following commands
```
cd <ferret>/contrib
./build_scalfmm
cd <ferret>
./configure --with-scalfmm=contrib/scalfmm
make -j2
```

For compilation of this O(N) solver on supercomputing resources, it is recommended to avoid the ./configure step and manually set the FERRET_HAVE_SCALFMM flag. The generation of the M2L compressors for the FMM-BEM should be done in serial and then can be read on any number of processors. The M2L compressor files do not depend on the mesh, but only the flags of the FMM-BEM method. 

More reading at X. Jiang et al J. Chem. Phys. 145, 064307 (2016).


### Who do I talk to? ###

* The following people are contributors to the Ferret project and can be reached by email for help:

* John Mangeri (john.mangeri@list.lu)
* Serge Nakhmanson (serge.nakhmanson@uconn.edu)
* Olle Heinonen (heinonen@anl.gov)
* Lukasz Kuna (lukasz.kuna@uconn.edu)

### Manual or Tutorials?

* A semi-updated tutorial and manual exists but is not in the repository at the present time due to its status as work-in-progress. If you would like a copy of it, please contact one of the contributors.

### MOOSE can be found at http://www.mooseframework.org ###
 * The GitHub discussion board for MOOSE (https://github.com/idaholab/moose/discussions) is a great resource for help with MOOSE or new implementations.

### Contributing ###
 * If you have an idea for a ferroelectric materials problem and would like to add it to the repository, then please submit a pull request *to the devel branch* with the relevant information. Generally, to submit a pull request, you should *fork* Ferret, type git checkout devel, git fetch origin, git rebase origin/devel, then make the changes/additions and a commit message, type git push. Then submit the pull request on the page (the next steps will be handled automatically).
