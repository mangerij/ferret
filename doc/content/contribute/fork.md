## Code Development

FERRET exists under the GNU General Public License. As such, it is *permissible* to make changes to the source code of FERRET, but +DO NOT+ make changes to `mangerij/ferret` as you don't have write access. To propose changes properly, you need to fork FERRET under your own GitHub account, make modifications in your fork, and submit a pull request to the relevant branch (in this case: `mangerij/ferret` ).

#### Fork FERRET

1. Log into your GitHub account, navigate to the [mangerij/ferret](https://github.com/mangerij/ferret) repository.
2. In the top-right corner of the page, click +Fork+

   !media media/fork_button.jpg
          style=width:300px;padding:auto;

3. After the forking process is done, you should be able to find your fork under your GitHub account, i.e. `https://github.com/YourGitHubUsername/ferret`.

#### Create a local clone of your fork

Now that you have a fork of `mangerij/ferret` under `YourGithubUsername/ferret`, you can setup the workflow for your fork. First, open a terminal, navigate to the FERRET repository, and verify your current remotes:

```bash
cd ~/projects/ferret
git remote -v
```

Your remote should look like the following:

```bash
origin	https://github.com/mangerij/ferret.git (fetch)
origin	https://github.com/mangerij/ferret.git (push)
```

Next, let's add your fork of FERRET to the remote:

```bash
git remote add YourRemoteName git@github.com:YourGitHubUsername/ferret.git
git remote -v
```

and the output should look like the following:

```bash
origin	https://github.com/mangerij/ferret.git (fetch)
origin	https://github.com/mangerij/ferret.git (push)
YourRemoteName	git@github.com:YourGitHubUsername/ferret.git (fetch)
YourRemoteName	git@github.com:YourGitHubUsername/ferret.git (push)
```

Next, you may checkout the `devel` branch of FERRET to your own development branch and let it track your own fork of FERRET, i.e.

```bash
git checkout origin/devel
git switch -c YourBranchName
git push -u YourRemoteName YourBranchName
```

Now, you may begin your code development. Just remember to add and commit often, i.e.

```bash
git add -A .
git commit -m 'Some commit message'
git push
```

Often times, during your development, the `devel` branch of `mangerij/ferret` may have new commits since your last checkout. To integrate new commits into your development branch, do:

```bash
git rebase -i origin devel
```

Additionally, please use the issues system to label your commit message. Github will automatically assign the commit log to the issue which allows for ease of tracking code changes.

If you would like to make changes to the documentation (and/or this website), follow [these instructions](getting_started/docs.md)
