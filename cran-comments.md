# mobr 2.0.0

This is a new package submission.

## Test environments
* windows X (local), R 4.0.2
* ubuntu 14.04.6 (personal server), R 3.6.2
* ubuntu 16.04.6 (on travis-ci), R 4.0.0
* win-builder (devel & release)

## R CMD check results

### Windows X
* No warnings or errors

### ubuntu (personal server)
* No warnings or errors

### ubuntu server (travis-ci)
* No warnings or errors

### on win-builder (devel & release)
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
* checking PDF version of manual without hyperrefs or index ... ERROR

The warning and error on win-builder are both related to a latex error in the
building of the package manual. On my local windows X environment, I successfully
built the pdf from source using R CMD Rd2pdf which suggests that this warning
may be erroneous.