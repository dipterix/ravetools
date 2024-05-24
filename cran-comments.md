## R CMD check results

0 errors | 0 warnings | 0 note

## Address CRAN submission check issues


Dear CRAN maintainer, I received Prof. Ripley's email on May 20, 2024:

```
Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_ravetools.html>.

Please correct before 2024-06-03 to safely retain your package on CRAN.

The CRAN Team
```

I was able to reproduce this issue.

The C++ code that caused this issue was completely removed from this new version.

The new version was tested under the same environment and no more `valgrind` error was reported.

Also checked with R-devel compiled with clang-ASAN, passed without error.

All windows (oldrel, release, and devel) are checked.
