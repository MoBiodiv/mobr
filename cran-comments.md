# mobr 2.0.0

This is a new package submission.

## Test environments
* windows X (local), R 4.0.2
* ubuntu 14.04.6 (personal server), R 3.6.2
* ubuntu 16.04.6 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

### Windows X
* No warnings or errors

### ubuntu (personal server)
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
* checking PDF version of manual without hyperrefs or index ... ERROR

### ubuntu server (travis-ci)
* No warnings or errors

### on win-builder (devel)
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
* checking PDF version of manual without hyperrefs or index ... ERROR

The warning and error encountered on the personal ubuntu server and win-builder
both related to the latex errors in the building of the package manual. These
errors did not occur in the local windows X or ubuntu travis ci tests which
suggests that this has less to do with the documentation of the package and more
to do with third party dependencies that are not installed on these other
testing environments.


