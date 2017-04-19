Ferret is an application within the [MOOSE](http://mooseframework.org) framework.
To build and use Ferret you will need to build MOOSE.
It might be a good idea to learn a few things about MOOSE, too.

Generally, you will clone Ferret this way:
```
git clone https://user@bitbucket.org/mesoscience/ferret.git
```
If you have write permissions, cloning this way will enable you to push to the repo.
Generally this means that you are a member of the `mesoscience` BitBucket team.
If you only have read permissions, you might want to clone this way:
```
git clone bitbucket.org/mesoscience/ferret.git
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

* John Mangeri (john.mangeri@uconn.edu)
* Serge Nakhmanson (serge.nakhmanson@uconn.edu)
* Olle Heinonen (heinonen@anl.gov)
* Lukasz Kuna (lukasz.kuna@uconn.edu

* MOOSE can be found at http://www.mooseframework.org
* The moose-users list at https://groups.google.com/forum/#!forum/moose-users is good for help with MOOSE!
