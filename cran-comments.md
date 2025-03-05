## R CMD check results

0 errors | 0 warnings | 0 note

In addition, the issues listed on CRAN has been fixed.

Problem:

```
Version: 0.2.1
Check: whether package can be installed
Result: WARN
  Found the following significant warnings:
    vcglib/wrap/callback.h:79:12: warning: address of stack memory associated with local variable 'formatted' returned [-Wreturn-stack-address]
  See ‘/data/gannet/ripley/R/packages/tests-clang/ravetools.Rcheck/00install.out’ for details.
  * used C++ compiler: ‘clang version 20.1.0-rc2’
Flavor: r-devel-linux-x86_64-fedora-clang
```

Solution: removed the function that caused this issue
