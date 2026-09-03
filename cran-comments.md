## R CMD check results

0 errors | 0 warnings | 0 note

This patch addresses the `clang`-`UBSAN` report in `build_vignettes.log`
("left shift of negative value" in the bundled `VCG` library). The cause was a
bit-flag allocator leak in `VCG`'s `CountNonManifoldEdgeFF`, which never
released the three flags it acquired; after a few calls in one session the
allocator's counter overran the sign bit. The same leak also silently inflated
the non-manifold edge count reported by `vcg_mesh_volume`, so this is a
correctness fix as well.

The bundled `VCG` sources are patched locally (the defect is still present
upstream) and the allocator's shift is now performed on an unsigned type.


