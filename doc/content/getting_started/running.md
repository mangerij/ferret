# Running your first model

## Understanding the code organization

Before running any model, it is important to understand what a "model" means in MOOSE. In the MOOSE ecosystem, the source code and the configuration of your model are separated.

The source codes of your model are written in C++, and they are compiled into an executable file. In FERRET the executable is called `ferret-opt`.

The configurations of your model are written in Hit format, and they are called "input files" with a `.i` extension. While your models may work for a wide variety of problems, an input file typically contains configurations for a *specific* problem. An input file may contain information about mesh, variables, weak form, constitutive models, boundary conditions, initial conditions, postprocessors, and so on.

!alert tip title=No hard-coding!
We highly recommend you to follow this idea of separating the model from its parameters. As this allows for more general code, less code duplication, maximum modularity/reusability, and minimum maintenance cost. And even better, it is suggested you use the automatic differentation tools within MOOSE to eliminate the need for complicated algebra required to generate the residual and jacobian terms. However, this particular advice was rather lost as FERRET was developed :)

## Command line interface

Before running a model, first make sure you have installed and compiled FERRET following [the installation guide](install.md).

To run a model, go to the directory of its input files, e.g.

```bash
cd ~/projects/ferret/tutorials/multidomain
```

and run the input file using the compiled executable:

```bash
../../ferret-opt -i multidomain.i
```

The first argument is the path to the FERRET executable. It is followed by a `-i` option meaning that we will supply our input file immediately after. `elasticity.i` is the path to the input file that we want to run.

To view the complete list of options, simply run the executable without any additional arguments:

```bash
../../ferret-opt
```

To run the model using multiple processors in parallel, use `mpiexec`, e.g.

```bash
mpiexec -n 2 ../../ferret-opt -i multidomain.i
```

where 2 is the number of processors you want to use.

After the calculation is finished, you can visualize the result following [the visualization guide](getting_started/paraview.md)

## Where to store your own models?

Since the implementation and the configuration of a model are separated, we need to consider them separately:

- The source code of the model, i.e. `.C` files, should be stored in `~/projects/ferret/src`.
- The header file of the model, i.e. `.h` files, should be stored in `~/projects/ferret/include`.
- The input files and all of their accompanying files shall be stored in a separate repository outside FERRET, unless you want to eventually merge into the official FERRET repository.

For your convenience, there is a symlink in the FERRET repository to a directory sitting in parallel. The symlink is

```bash
~/projects/ferret/examples -> ~/projects/ferret_examples
```

meaning that if you create a directory outside FERRET called `ferret_examples`, you will be able to access of its content from `~/projects/ferret/examples` as if they sit right there, but all content within that directory are not under version control (you can still set up one on your own, and it is actually highly recommended).

!alert tip title=Output naming
Be careful with your output naming convention! If multiple simulations are run in the same folder to compare parameters but the names are not modified then they will be overwritten by the executable! This can cause some headache if calculations are heavy so always make sure you change the names of outputs if you want to compare results.
