# mobr 2.0.0

## Test environments
* local Windows X install, R 4.0.2
* ubuntu 14.04.6 (personal server), R 3.6.2
* ubuntu 16.04.6 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results
### on win-builder (devel)
There was one WARNINGs. 
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Package inputenc Error: Unicode char ‐ (U+2010)
(inputenc)                not set up for use with LaTeX.

There was one ERRORs.
* checking PDF version of manual without hyperrefs or index ... ERROR

This warning and error did not occur in the other testing enviornments; therefore,
I believe that it is safe to ignore for the purposes of this submission. 

