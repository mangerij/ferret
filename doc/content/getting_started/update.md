# Update FERRET

MOOSE changes frequently. Since FERRET is integrated into the Continuous-Integration and Verification Test Environment (CIVET, https://civet.inl.gov/), the MOOSE team may make periodic changes to FERRET to keep it up-to-date.

Therefore, to make sure your changes are not going to be made deprecated or conflict with newer commits of MOOSE, you will need to update your Conda environment as well as your local repository frequently to keep up with the remote.

To update the Conda environment, use

```bash
conda activate moose
conda update --all
```

To update FERRET, use

```bash
git pull
git submodule update --init --recursive
```

Then, recompile FERRET using

```bash
make -j N
```

> Sometimes after updating FERRET and all submodules, you may get compile errors. Most likely that is due to stale objects floating around in the repository. To clean stale objects, +make sure you have everything committed+, and do
>
> ```bash
> git submodule foreach --recursive git clean -xfd
> git clean -xfd
> ```
>
> After which you should be able to recompile FERRET
>
> ```bash
> make -j N
> ```
