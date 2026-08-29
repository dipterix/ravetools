## R CMD check results

0 errors | 0 warnings | 0 note

Currently CRAN's servers report "OK" on all flavors, but self-checking found some critical bugs.

This patch aims to improve the robustness by fixing:

* A hidden `valgrind` issue reported by `rhub`;
* A syntax issue that prevents the code from compiling on `C++11`;
* A bug that corrupted mesh in `vcg_smooth_implicit` may lead to segfault.


