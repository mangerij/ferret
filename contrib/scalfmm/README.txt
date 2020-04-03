
ScalFmm, Inria, Please read the licence.

---------------------------------------------------
---------------------------------------------------

To compile:
==========
# Go to 
cd scalfmm/Build
# Use cmake first by
cmake ..
# Or if you want to use MPI
cmake .. -DSCALFMM_USE_MPI=ON


# Configure, for example with:
ccmake ..
# turn on/off the options you want
# For example, it is advised to build the tests and have a look in scalfmm/Tests/
# Finally you can build by
make
# And access executables in scalfmm/Build/Tests/{Release,Debug}/.....

To compile with StarPU you can pass, as an example:
cmake .. -DSCALFMM_USE_STARPU=ON -DSCALFMM_USE_CUDA=OFF -DSCALFMM_USE_OPENCL=OFF  \
               -DHWLOC_DIR=/home/berenger/Download/hwloc-1.10.0/install/      \
               -DSTARPU_DIR=/home/berenger/Download/starpu-work/StarPU/installwithfxt


---------------------------------------------------
---------------------------------------------------

Build the doc:
=============
In scalfmm/Doc you can find several pdf and .tex file about
the implementation, kernels and data structure.

# The developer documentation is generated using DOxygen: 
cd scalfmm/Build
cmake .. -DSCALFMM_BUILD_DOC=ON
# Tape 
make doc
# This will create a Html dir
browser scalfmm/Build/Doc/html/index.html


---------------------------------------------------
---------------------------------------------------

Getting help and having news from us:
====================================

You can subscribe to the scalfmm users mailing list ( scalfmm-public-users@lists.gforge.inria.fr,  http://lists.gforge.inria.fr/cgi-bin/mailman/listinfo/scalfmm-public-users ). Very low trafic (~ 2 mails per year) just to know when a new version or an improvement is available.

Contact the developers at : scalfmm-public-support@lists.gforge.inria.fr


---------------------------------------------------
---------------------------------------------------

What inside :
=============
- Src : The Core of Scalfmm is under the Src directory. Users should not need to modify the source.
One can want to implement its own kernel or even its own parallelization whithout modifying the sources.
- Data : example of particles distributions
- Examples : examples of very common usage of Scalfmm
- Doc : should contains the generated Doc
- UTests : contains some unit tests (it can be a good example to understand some features)
- Tests : examples to know how to use scalfmm/put particles in the tree/iterate on the tree...
- Utils : some scripts to work with the data files.


---------------------------------------------------
---------------------------------------------------
Concerning the 1.3 to 1.4 API update:
=====================================
For example the previous M2L interface was:
void M2L(CellClass* const FRestrict inLocal, const CellClass* inInteractions[], const int /*inSize*/, const int inLevel) {

But the new one is:
void M2L(CellClass* const FRestrict inLocal, const CellClass* inInteractions[],
 const int neighborPositions[], const int inSize, const int inLevel) override {

In the previous version we were iterating by
for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
 // if interaction exits
 if(inInteractions[idxNeigh]){
     inInteractions[idxNeigh]


In the new one we do it by
for(int idxExistingNeigh = 0 ; idxExistingNeigh < inSize ; ++idxExistingNeigh){
 const int idxNeigh = neighborPositions[idxExistingNeigh];
 inInteractions[idxExistingNeigh]
