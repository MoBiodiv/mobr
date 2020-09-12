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
* No warnings or errors

## CRAN package review response

Thanks,

* Please omit the redundant "in R" from the title.

done

* If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

done


* Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or
'F' as vector names.)

done

* You have examples for unexported functions.
    groups_panel2() in:
       ind_rare_perm.Rd
    plotStacked() in:
       sphere_dist.Rd
  Please either omit these examples or export the functions.

The examples were removed.

* \dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.

`\dontrun{}` has been changed to `\donttest{}` in all cases in which
either the example takes more than 5 seconds to run or the example
requires an external R package that is not a direct package dependancy.

* You write information messages to the console that cannot be easily
suppressed.
It is more R like to generate objects that can be used to extract the
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning()  or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console.
(except for print, summary, interactive functions)

I have replaced `print()` with `message()` where necessary.

* Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)       # code line i
on.exit(par(oldpar))                    # code line i + 1
...
par(mfrow=c(2,2))                       # somewhere after
...
e.g.:
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.

done

* Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

done

Thank you for the constructive and helpful package review. 