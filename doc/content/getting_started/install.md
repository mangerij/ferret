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

Next, clone FERRET to your local projects directory:

```bash
mkdir ~/projects
cd ~/projects
git clone https://github.com/mangerij/ferret.git
./configure
cd ferret
git checkout master
git submodule update --init --recursive
```

These commands should download a copy of FERRET and a copy of MOOSE (as a submodule) to your local projects directory.

> +\[Optional\]+ In Step 1, if you didn't choose to include moose-libmesh in your Conda environment (typically on HPC systems where Conda is not suggested), you need to compile PETSc and libMesh using
>
> ```bash
> ./moose/scripts/update_and_rebuild_petsc.sh
> ./moose/scripts/update_and_rebuild_libmesh.sh
> ```

### Step 3: Compile FERRET style=line-height:150%;

Next, you can compile FERRET using

```bash
make -j N
```

where `N` is the number of processors you want to use to compile FERRET in parallel.

> +\[Optional\]+ To make sure FERRET is working properly, run the regression tests:
>
> ```bash
> ./run_tests -j N
> ```


### Step 4 (optional): Compile with ScalFMM to use the fast-multipole boundary element method


The ./configure step can include additional libraries such as BOOST which will allow compile of MOOSE objects with mathematical special functions. Or it can include an advanced fast-multipole boundary element method (FMM-BEM) developed by Prof. Xikai Jiang and co-workers. In order to use the BEM, you need to compile ScalFMM with the following commands

```bash
cd <ferret>/contrib
./build_scalfmm
cd <ferret>
./configure --with-scalfmm=contrib/scalfmm
make -j2
```

For compilation of this $O(\mathrm{N})$ solver on supercomputing resources, it is recommended to avoid the ./configure step and manually set the FERRET_HAVE_SCALFMM flag. The generation of the M2L compressors for the FMM-BEM should be done in serial and then can be read on any number of processors. The M2L compressor files do not depend on the mesh, but only the flags of the FMM-BEM method.

This is an experimental feature. More reading at [!cite](Jiang2016).
