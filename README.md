### We have moved from BitBucket! ###

As of July 23rd, 2019 the BitBucket repository is no longer supported or current with the MOOSE master. All official changes will be made here at GitHub and Ferret business will continue as usual on CIVET.

Ferret still exists under the open-source GNU License but to access the repository, please contact one of the developers.

### What is Ferret and how can I compile it? ###

Ferret is an application within the [MOOSE](http://mooseframework.org) framework.
To build and use Ferret you will need to build MOOSE.
It might be a good idea to learn a few things about MOOSE, too.

Generally, you will clone Ferret this way:
```
git clone https://github.com/mangerij/ferret.git
```
Once MOOSE is built, Ferret can usually be built simply as
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
NOTE: The ./configure step can alternatively include a BOOST directory which will allow compile of MOOSE objects with mathematical special functions.

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
 * At the moment, please cite Nanoscale, 2017, 9, 1616-1624 if you use Ferret. A methods paper is currently in the works.
 * We have a google group now! https://groups.google.com/forum/#!forum/ferret-users/new
