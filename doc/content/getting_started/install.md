# Install FERRET

FERRET leverages the integration of multiple high-quality open-source scientific computing projects.

- FERRET is built upon [MOOSE](https://mooseframework.inl.gov/), an object oriented parallel FEM framework.
- MOOSE is built upon [libMesh](http://libmesh.github.io/), a C++ FEM library, to utilize its FEM basics.
- Both MOOSE and libMesh rely on many other tools, such as [PETSc](https://www.mcs.anl.gov/petsc/).

Therefore, it can be complicated to get FERRET to work -especially on a high-performance computing (HPC) cluster.

Fortunately, the installation process has been thoroughly tested for many different operating systems and HPC configurations. The installation process can be summarized in the following steps,
each of which should only take a handful of commands.

### Step 1: Install dependencies style=line-height:150%;

First, follow [these instructions](getting_started/conda.md) to install environment packages for MOOSE and FERRET.

### Step 2: Clone FERRET style=line-height:150%;

Next, clone both MOOSE and FERRET into your local projects directory:

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/idaholab/moose.git
git clone https://github.com/mangerij/ferret.git
cd ferret
git checkout master
./configure
```

!alert! note title=MOOSE is a sibling directory, not a submodule
FERRET does +not+ vendor MOOSE. The `Makefile` locates it with `MOOSE_DIR ?= $(shell dirname \`pwd\`)/moose`, so MOOSE must sit +next to+ FERRET (for example `~/projects/moose` and `~/projects/ferret`). If you keep MOOSE somewhere else, export `MOOSE_DIR` to point at it before compiling:

```bash
export MOOSE_DIR=/path/to/moose
```
!alert-end!

FERRET's only git submodule is ScalFMM, which is optional and covered in Step 4. You do not need `git submodule update` for a standard build.

> +\[Optional\]+ In Step 1, if you didn't choose to include moose-libmesh in your Conda environment (typically on HPC systems where Conda is not suggested), you need to compile PETSc and libMesh from within the MOOSE clone:
>
> ```bash
> cd ~/projects/moose
> ./scripts/update_and_rebuild_petsc.sh
> ./scripts/update_and_rebuild_libmesh.sh
> ```

!alert note title=Important!
If you are trying to clone Ferret (or MOOSE) in the Windows Sublayer for Linux (WSL), you may find that the handshake fails and the repository isn't initialized properly or outright fails. This is due to known issues with the HNS container networking protocols. See [here](https://stackoverflow.com/questions/38618885/error-rpc-failed-curl-transfer-closed-with-outstanding-read-data-remaining ) for some work-arounds.


### Step 3: Compile FERRET style=line-height:150%;

Next, from inside the FERRET directory, compile it using

```bash
make -j N
```

where `N` is the number of processors you want to use to compile FERRET in parallel. This produces the `ferret-opt` executable in the FERRET root directory.

> +\[Optional\]+ To make sure FERRET is working properly, run the regression tests:
>
> ```bash
>  ./run_tests -j N
> ```

### Step 4 (optional): Compile with ScalFMM to use the fast-multipole boundary element method


The `./configure` step can include additional libraries such as BOOST which will allow compile of MOOSE objects with mathematical special functions. Or it can include an advanced fast-multipole boundary element method (FMM-BEM) developed by Prof. Xikai Jiang and co-workers. In order to use the BEM, you need to compile ScalFMM with the following commands

```bash
cd <ferret>/scripts
sh ./build_scalfmm_local
(may encounter an error)
sh ./build_scalfmm_local
cd <ferret>
./configure --with-scalfmm=ScalFMM
make -j2
```

For compilation of this $O(\mathrm{N})$ solver on supercomputing resources, it is recommended to avoid the `./configure` step and manually set the `FERRET_HAVE_SCALFMM` flag. The generation of the M2L compressors for the FMM-BEM should be done in serial and then can be read on any number of processors. The `M2L` compressor files do not depend on the `Mesh`, but only the flags of the FMM-BEM method.

This is an experimental feature at the moment and we plan to extend this capability further. More reading at [!cite](Jiang2016).
