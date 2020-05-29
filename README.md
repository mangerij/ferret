### We have moved from BitBucket! ###

As of July 23rd, 2019 the BitBucket repository is no longer supported. All official changes will be made here on GitHub and Ferret business will continue as usual.
The continuous integration and testing server (CIVET: https://civet.inl.gov/) ensures this version is current with the MOOSE master.

Ferret still exists under the open-source GNU License but to access the repository you will need an RSA deploy key registered with our application. We freely provide these keys in order to track repository usage and measure impact. The deploy key itself allows write-access on request.
Please contact one of the developers for assistance in generating one.

After obtaining a key, please cite Mangeri et al Nanoscale, 2017, 9, 1616-1624 in the *main text* of your paper if you use Ferret. A methods paper is in the works. Finally, one more comment about usage: any issues with legacy inputs not running on this version of Ferret should be resolved with corresponding authors of appropriate papers. As this code belongs under GNU licensing, we cannot possibly support all of the custom spin-offs and additions to the code. We wish we could but sometimes cooperation and a mature interest in the open-source approach to science is less than optimal.

### What is Ferret and how can I compile it? ###

Ferret is an application within the finite element [MOOSE](http://mooseframework.org) framework useful for simulating the ferroic nanostructure. More information is available at our website https://ferretnano.weebly.com/ which is still under-construction. To build and use Ferret you will need to build MOOSE. It might be a good idea to learn a few things about MOOSE, too.

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
Once MOOSE is built, Ferret can be built simply with
```
cd <ferret>
./configure
make METHOD=opt MOOSE_DIR=<moose>
```
or
```
make METHOD=opt MOOSE_DIR=<moose> -j 4
```
to make the compile go faster.

Here `<ferret>` is the location of the Ferret clone, and `<moose>` is the MOOSE clone.
Naturally, use an appropriate `METHOD` -- one of those that was used to build MOOSE
(e.g., `dbg` or `opt`).

If libMesh is built in a non-standard location (i.e., NOT built inside the MOOSE tree 
using the `<moose>/scripts/update_and_rebuild_libmesh.sh`) then `LIBMESH_DIR` also has 
to be set:
```
make METHOD=opt MOOSE_DIR=<moose> LIBMESH_DIR=<libmesh>
```
### Additional information: ###

The ./configure step can include additional libraries such as BOOST which will allow compile of MOOSE objects with mathematical special functions.
Or it can include an advanced boundary element method (BEM) developed by Prof. Xikai Jiang. In order to use the BEM, you need to compile ScalFMM with the following commands
```
cd <ferret>/contrib
./build_scalfmm
./configure --with-scalfmm=contrib/scalfmm
make -j2
```
and then run the ./configure step as above.

### Who do I talk to? ###

* The following people are contributors to the Ferret project and can be reached by email for help:

* John Mangeri (mangeri@fzu.cz)
* Serge Nakhmanson (serge.nakhmanson@uconn.edu)
* Olle Heinonen (heinonen@anl.gov)
* Lukasz Kuna (lukasz.kuna@uconn.edu)

### Manual or Tutorials?

* A semi-updated tutorial and manual exists but is not in the repository at the present time due to its status as work-in-progress. If you would like a copy of it, please contact one of the contributors.

### MOOSE can be found at http://www.mooseframework.org ###
 * The moose-users list at https://groups.google.com/forum/#!forum/moose-users is a great link for help with MOOSE or new implementations.

### Contributing ###
 * If you have an idea for some ferroelectric materials problem and would like to add it to the repository, then please submit a pull request *to the devel branch* with the relevant information. Generally, to submit a pull request, you should *fork* Ferret, type git checkout devel, git fetch origin, git rebase origin/devel, then make the changes/additions and a commit message, type git push. Then submit the pull request on the page (the next steps will be handled automatically).
 * We have a google group now! https://groups.google.com/forum/#!forum/ferret-users/new
