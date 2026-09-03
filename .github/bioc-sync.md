# Keeping GitHub and Bioconductor in sync

The GitHub default branch (`master`, until it is renamed) and Bioconductor's
`devel` branch must point to the same commit.  Bioconductor is the deployment
remote; GitHub is the collaboration mirror.  Do not force-push either branch.

## One-time local setup

```sh
git remote add upstream git@git.bioconductor.org:packages/cosmosR.git
git fetch --prune origin
git fetch --prune upstream
```

`origin` should be the GitHub repository and `upstream` the Bioconductor
repository.  Confirm that the SSH key used by this computer is registered with
Bioconductor before attempting a push.

## Before every package update

Start from a clean working tree and fetch both remotes.  Merge the current
Bioconductor `devel` branch before making the package change, then bump the
patch (`z`) component of `Version` in `DESCRIPTION` for the final commit.

```sh
git fetch --prune origin
git fetch --prune upstream
git switch master
git merge --ff-only origin/master
git merge upstream/devel
```

After tests pass, publish the exact same commit to both remotes:

```sh
git push upstream HEAD:devel
git push origin HEAD:master
```

The `Verify Bioconductor synchronization` workflow checks strict commit-ID
equality after pushes and daily.  A failure means that a Bioconductor release
transition or a GitHub-only commit has not yet been mirrored; resolve it with a
normal merge and an ordinary (non-force) push.

## Recommended follow-up

Bioconductor recommends renaming the GitHub default branch to `devel`.  Once
that is done in GitHub's branch settings, replace `master` with `devel` in the
commands above and update branch-protection rules and any branch-specific CI
conditions.  Matching branch names removes the `HEAD:devel` mapping and makes
accidental divergence less likely.
