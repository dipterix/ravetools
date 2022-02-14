## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Notes, fixes

```
If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")
```

Thanks, two citations have been added to the `inst/CITATION` file, and the major one has been added to `DESCRIPTION`.


```
Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags:
     pwelch.Rd: \value

Please fix and resubmit.
```

Thanks, detailed value description has been added to the `pwelch` function.


In addition, 
* All examples and tests are using 2 cores to comply with the CRAN policy.
* A critical `C++` bug has been fixed on `Win-32` system.
