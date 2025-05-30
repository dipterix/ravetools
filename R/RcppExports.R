# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

baselineArray <- function(x, bl, dims, bldims, tidx, per, rest, method) {
    .Call(`_ravetools_baselineArray`, x, bl, dims, bldims, tidx, per, rest, method)
}

bucketFillVolume <- function(volume, x, y, z, fill) {
    .Call(`_ravetools_bucketFillVolume`, volume, x, y, z, fill)
}

collapser_cplx <- function(x, keep, method = 1L, average = 0L) {
    .Call(`_ravetools_collapser_cplx`, x, keep, method, average)
}

collapser_real <- function(x, keep, method = 1L, average = 0L) {
    .Call(`_ravetools_collapser_real`, x, keep, method, average)
}

columnQuantile <- function(x, prob, naRm) {
    .Call(`_ravetools_columnQuantile`, x, prob, naRm)
}

columnMedian <- function(x, naRm) {
    .Call(`_ravetools_columnMedian`, x, naRm)
}

fastColMeans <- function(x, col, xDim) {
    .Call(`_ravetools_fastColMeans`, x, col, xDim)
}

quickQuantile <- function(x, prob, naRm) {
    .Call(`_ravetools_quickQuantile`, x, prob, naRm)
}

quickMedian <- function(x, naRm) {
    .Call(`_ravetools_quickMedian`, x, naRm)
}

fastcov <- function(x1, x2, col1 = NULL, col2 = NULL, df = -1.0) {
    .Call(`_ravetools_fastcov`, x1, x2, col1, col2, df)
}

fftw_r2c <- function(data, HermConj = 1L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_r2c`, data, HermConj, fftwplanopt, ret)
}

fftw_c2c <- function(data, inverse = 0L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_c2c`, data, inverse, fftwplanopt, ret)
}

fftw_c2r <- function(data, HermConj = 1L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_c2r`, data, HermConj, fftwplanopt, ret)
}

mvfftw_r2c <- function(data, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_mvfftw_r2c`, data, fftwplanopt, ret)
}

mvfft_c2r <- function(data, fftwplanopt = 0L, retrows = 0L, ret = NULL) {
    .Call(`_ravetools_mvfft_c2r`, data, fftwplanopt, retrows, ret)
}

fftw_r2c_2d <- function(data, HermConj = 1L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_r2c_2d`, data, HermConj, fftwplanopt, ret)
}

fftw_c2c_2d <- function(data, inverse = 0L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_c2c_2d`, data, inverse, fftwplanopt, ret)
}

fftw_r2c_3d <- function(data, HermConj = 1L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_r2c_3d`, data, HermConj, fftwplanopt, ret)
}

fftw_c2c_3d <- function(data, inverse = 0L, fftwplanopt = 0L, ret = NULL) {
    .Call(`_ravetools_fftw_c2c_3d`, data, inverse, fftwplanopt, ret)
}

conjugate <- function(data) {
    .Call(`_ravetools_conjugate`, data)
}

cpp_filter <- function(b, a, x, z) {
    .Call(`_ravetools_cpp_filter`, b, a, x, z)
}

Matrix4__new <- function() {
    .Call(`_ravetools_Matrix4__new`)
}

Matrix4__to_array <- function(self) {
    .Call(`_ravetools_Matrix4__to_array`, self)
}

Matrix4__from_array <- function(self, array, offset = 0L) {
    invisible(.Call(`_ravetools_Matrix4__from_array`, self, array, offset))
}

Matrix4__identity <- function(self) {
    invisible(.Call(`_ravetools_Matrix4__identity`, self))
}

Matrix4__copy <- function(self, m) {
    invisible(.Call(`_ravetools_Matrix4__copy`, self, m))
}

Matrix4__copy_position <- function(self, m) {
    invisible(.Call(`_ravetools_Matrix4__copy_position`, self, m))
}

Matrix4__extract_basis <- function(self, x, y, z) {
    invisible(.Call(`_ravetools_Matrix4__extract_basis`, self, x, y, z))
}

Matrix4__make_basis <- function(self, x, y, z) {
    invisible(.Call(`_ravetools_Matrix4__make_basis`, self, x, y, z))
}

Matrix4__extract_rotation <- function(self, m) {
    invisible(.Call(`_ravetools_Matrix4__extract_rotation`, self, m))
}

Matrix4__look_at <- function(self, eye, target, up) {
    invisible(.Call(`_ravetools_Matrix4__look_at`, self, eye, target, up))
}

Matrix4__multiply_matrices <- function(self, a, b) {
    invisible(.Call(`_ravetools_Matrix4__multiply_matrices`, self, a, b))
}

Matrix4__multiply_scalar <- function(self, s) {
    invisible(.Call(`_ravetools_Matrix4__multiply_scalar`, self, s))
}

Matrix4__determinant <- function(self) {
    .Call(`_ravetools_Matrix4__determinant`, self)
}

Matrix4__transpose <- function(self) {
    invisible(.Call(`_ravetools_Matrix4__transpose`, self))
}

Matrix4__set_position <- function(self, x, y, z) {
    invisible(.Call(`_ravetools_Matrix4__set_position`, self, x, y, z))
}

Matrix4__invert <- function(self) {
    invisible(.Call(`_ravetools_Matrix4__invert`, self))
}

Matrix4__scale <- function(self, v) {
    invisible(.Call(`_ravetools_Matrix4__scale`, self, v))
}

Matrix4__get_max_scale_on_axis <- function(self) {
    .Call(`_ravetools_Matrix4__get_max_scale_on_axis`, self)
}

Matrix4__make_translation <- function(self, x, y, z) {
    invisible(.Call(`_ravetools_Matrix4__make_translation`, self, x, y, z))
}

Matrix4__make_rotation_x <- function(self, theta) {
    invisible(.Call(`_ravetools_Matrix4__make_rotation_x`, self, theta))
}

Matrix4__make_rotation_y <- function(self, theta) {
    invisible(.Call(`_ravetools_Matrix4__make_rotation_y`, self, theta))
}

Matrix4__make_rotation_z <- function(self, theta) {
    invisible(.Call(`_ravetools_Matrix4__make_rotation_z`, self, theta))
}

Matrix4__make_rotation_axis <- function(self, axis, angle) {
    invisible(.Call(`_ravetools_Matrix4__make_rotation_axis`, self, axis, angle))
}

Matrix4__make_scale <- function(self, x, y, z) {
    invisible(.Call(`_ravetools_Matrix4__make_scale`, self, x, y, z))
}

Matrix4__make_shear <- function(self, xy, xz, yx, yz, zx, zy) {
    invisible(.Call(`_ravetools_Matrix4__make_shear`, self, xy, xz, yx, yz, zx, zy))
}

Matrix4__make_perspective <- function(self, left, right, top, bottom, near, far) {
    invisible(.Call(`_ravetools_Matrix4__make_perspective`, self, left, right, top, bottom, near, far))
}

Matrix4__make_orthographic <- function(self, left, right, top, bottom, near, far) {
    invisible(.Call(`_ravetools_Matrix4__make_orthographic`, self, left, right, top, bottom, near, far))
}

Quaternion__new <- function() {
    .Call(`_ravetools_Quaternion__new`)
}

Quaternion__set <- function(self, x, y, z, w) {
    invisible(.Call(`_ravetools_Quaternion__set`, self, x, y, z, w))
}

Quaternion__copy <- function(self, quaternion) {
    invisible(.Call(`_ravetools_Quaternion__copy`, self, quaternion))
}

Quaternion__to_array <- function(self) {
    .Call(`_ravetools_Quaternion__to_array`, self)
}

Quaternion__getX <- function(self) {
    .Call(`_ravetools_Quaternion__getX`, self)
}

Quaternion__setX <- function(self, v) {
    invisible(.Call(`_ravetools_Quaternion__setX`, self, v))
}

Quaternion__getY <- function(self) {
    .Call(`_ravetools_Quaternion__getY`, self)
}

Quaternion__setY <- function(self, v) {
    invisible(.Call(`_ravetools_Quaternion__setY`, self, v))
}

Quaternion__getZ <- function(self) {
    .Call(`_ravetools_Quaternion__getZ`, self)
}

Quaternion__setZ <- function(self, v) {
    invisible(.Call(`_ravetools_Quaternion__setZ`, self, v))
}

Quaternion__getW <- function(self) {
    .Call(`_ravetools_Quaternion__getW`, self)
}

Quaternion__setW <- function(self, v) {
    invisible(.Call(`_ravetools_Quaternion__setW`, self, v))
}

Quaternion__set_from_axis_angle <- function(self, axis, angle) {
    invisible(.Call(`_ravetools_Quaternion__set_from_axis_angle`, self, axis, angle))
}

Quaternion__set_from_rotation_matrix <- function(self, m) {
    invisible(.Call(`_ravetools_Quaternion__set_from_rotation_matrix`, self, m))
}

Quaternion__set_from_unit_vectors <- function(self, v_from, v_to) {
    invisible(.Call(`_ravetools_Quaternion__set_from_unit_vectors`, self, v_from, v_to))
}

Quaternion__angle_to <- function(self, q) {
    .Call(`_ravetools_Quaternion__angle_to`, self, q)
}

Quaternion__rotate_towards <- function(self, q, step) {
    invisible(.Call(`_ravetools_Quaternion__rotate_towards`, self, q, step))
}

Quaternion__slerp <- function(self, qb, t) {
    invisible(.Call(`_ravetools_Quaternion__slerp`, self, qb, t))
}

Quaternion__identity <- function(self) {
    invisible(.Call(`_ravetools_Quaternion__identity`, self))
}

Quaternion__invert <- function(self) {
    invisible(.Call(`_ravetools_Quaternion__invert`, self))
}

Quaternion__conjugate <- function(self) {
    invisible(.Call(`_ravetools_Quaternion__conjugate`, self))
}

Quaternion__dot <- function(self, v) {
    .Call(`_ravetools_Quaternion__dot`, self, v)
}

Quaternion__length_squared <- function(self) {
    .Call(`_ravetools_Quaternion__length_squared`, self)
}

Quaternion__length <- function(self) {
    .Call(`_ravetools_Quaternion__length`, self)
}

Quaternion__normalize <- function(self) {
    invisible(.Call(`_ravetools_Quaternion__normalize`, self))
}

Quaternion__multiply <- function(self, q) {
    invisible(.Call(`_ravetools_Quaternion__multiply`, self, q))
}

Quaternion__premultiply <- function(self, q) {
    invisible(.Call(`_ravetools_Quaternion__premultiply`, self, q))
}

Quaternion__multiply_quaternions <- function(self, a, b) {
    invisible(.Call(`_ravetools_Quaternion__multiply_quaternions`, self, a, b))
}

Vector3__new <- function() {
    .Call(`_ravetools_Vector3__new`)
}

Vector3__from_array <- function(self, array, offset = 0L, n_elems = -1L) {
    invisible(.Call(`_ravetools_Vector3__from_array`, self, array, offset, n_elems))
}

Vector3__resize <- function(self, n_elems) {
    invisible(.Call(`_ravetools_Vector3__resize`, self, n_elems))
}

Vector3__get_size <- function(self) {
    .Call(`_ravetools_Vector3__get_size`, self)
}

Vector3__to_array <- function(self, n_skip = 0L, max_n_elems = -1L) {
    .Call(`_ravetools_Vector3__to_array`, self, n_skip, max_n_elems)
}

Vector3__set_scalar <- function(self, value) {
    invisible(.Call(`_ravetools_Vector3__set_scalar`, self, value))
}

Vector3__set_x <- function(self, value) {
    invisible(.Call(`_ravetools_Vector3__set_x`, self, value))
}

Vector3__set_y <- function(self, value) {
    invisible(.Call(`_ravetools_Vector3__set_y`, self, value))
}

Vector3__set_z <- function(self, value) {
    invisible(.Call(`_ravetools_Vector3__set_z`, self, value))
}

Vector3__get_x <- function(self, i) {
    .Call(`_ravetools_Vector3__get_x`, self, i)
}

Vector3__get_y <- function(self, i) {
    .Call(`_ravetools_Vector3__get_y`, self, i)
}

Vector3__get_z <- function(self, i) {
    .Call(`_ravetools_Vector3__get_z`, self, i)
}

Vector3__get_item <- function(self, i) {
    .Call(`_ravetools_Vector3__get_item`, self, i)
}

Vector3__copy <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__copy`, self, v))
}

Vector3__add <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__add`, self, v))
}

Vector3__add_scalar <- function(self, s) {
    invisible(.Call(`_ravetools_Vector3__add_scalar`, self, s))
}

Vector3__add_vectors <- function(self, a, b) {
    invisible(.Call(`_ravetools_Vector3__add_vectors`, self, a, b))
}

Vector3__add_scaled <- function(self, v, s) {
    invisible(.Call(`_ravetools_Vector3__add_scaled`, self, v, s))
}

Vector3__sub <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__sub`, self, v))
}

Vector3__sub_scalar <- function(self, s) {
    invisible(.Call(`_ravetools_Vector3__sub_scalar`, self, s))
}

Vector3__sub_vectors <- function(self, a, b) {
    invisible(.Call(`_ravetools_Vector3__sub_vectors`, self, a, b))
}

Vector3__multiply <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__multiply`, self, v))
}

Vector3__multiply_scalar <- function(self, s) {
    invisible(.Call(`_ravetools_Vector3__multiply_scalar`, self, s))
}

Vector3__multiply_vectors <- function(self, a, b) {
    invisible(.Call(`_ravetools_Vector3__multiply_vectors`, self, a, b))
}

Vector3__apply_matrix3 <- function(self, m) {
    invisible(.Call(`_ravetools_Vector3__apply_matrix3`, self, m))
}

Vector3__apply_matrix4 <- function(self, m) {
    invisible(.Call(`_ravetools_Vector3__apply_matrix4`, self, m))
}

Vector3__apply_quaternion <- function(self, q) {
    invisible(.Call(`_ravetools_Vector3__apply_quaternion`, self, q))
}

Vector3__transform_direction <- function(self, m) {
    invisible(.Call(`_ravetools_Vector3__transform_direction`, self, m))
}

Vector3__divide <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__divide`, self, v))
}

Vector3__divide_scalar <- function(self, s) {
    invisible(.Call(`_ravetools_Vector3__divide_scalar`, self, s))
}

Vector3__min <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__min`, self, v))
}

Vector3__max <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__max`, self, v))
}

Vector3__clamp <- function(self, min, max) {
    invisible(.Call(`_ravetools_Vector3__clamp`, self, min, max))
}

Vector3__floor <- function(self) {
    invisible(.Call(`_ravetools_Vector3__floor`, self))
}

Vector3__ceil <- function(self) {
    invisible(.Call(`_ravetools_Vector3__ceil`, self))
}

Vector3__round <- function(self) {
    invisible(.Call(`_ravetools_Vector3__round`, self))
}

Vector3__round_to_zero <- function(self) {
    invisible(.Call(`_ravetools_Vector3__round_to_zero`, self))
}

Vector3__negate <- function(self) {
    invisible(.Call(`_ravetools_Vector3__negate`, self))
}

Vector3__dot <- function(self, v) {
    .Call(`_ravetools_Vector3__dot`, self, v)
}

Vector3__length_squared <- function(self) {
    .Call(`_ravetools_Vector3__length_squared`, self)
}

Vector3__length <- function(self) {
    .Call(`_ravetools_Vector3__length`, self)
}

Vector3__length_manhattan <- function(self) {
    .Call(`_ravetools_Vector3__length_manhattan`, self)
}

Vector3__normalize <- function(self) {
    invisible(.Call(`_ravetools_Vector3__normalize`, self))
}

Vector3__set_length <- function(self, length) {
    invisible(.Call(`_ravetools_Vector3__set_length`, self, length))
}

Vector3__lerp <- function(self, v, alpha) {
    invisible(.Call(`_ravetools_Vector3__lerp`, self, v, alpha))
}

Vector3__lerp_vectors <- function(self, v1, v2, alpha) {
    invisible(.Call(`_ravetools_Vector3__lerp_vectors`, self, v1, v2, alpha))
}

Vector3__cross <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__cross`, self, v))
}

Vector3__cross_vectors <- function(self, a, b) {
    invisible(.Call(`_ravetools_Vector3__cross_vectors`, self, a, b))
}

Vector3__project_on_vector <- function(self, v) {
    invisible(.Call(`_ravetools_Vector3__project_on_vector`, self, v))
}

Vector3__project_on_plane <- function(self, planeNormal) {
    invisible(.Call(`_ravetools_Vector3__project_on_plane`, self, planeNormal))
}

Vector3__reflect <- function(self, normal) {
    invisible(.Call(`_ravetools_Vector3__reflect`, self, normal))
}

Vector3__angle_to <- function(self, v) {
    .Call(`_ravetools_Vector3__angle_to`, self, v)
}

Vector3__distance_to <- function(self, v) {
    .Call(`_ravetools_Vector3__distance_to`, self, v)
}

Vector3__distance_to_squared <- function(self, v) {
    .Call(`_ravetools_Vector3__distance_to_squared`, self, v)
}

Vector3__distance_to_manhattan <- function(self, v) {
    .Call(`_ravetools_Vector3__distance_to_manhattan`, self, v)
}

Vector3__set_from_spherical_coords <- function(self, radius, phi, theta) {
    invisible(.Call(`_ravetools_Vector3__set_from_spherical_coords`, self, radius, phi, theta))
}

Vector3__set_from_matrix_position <- function(self, m) {
    invisible(.Call(`_ravetools_Vector3__set_from_matrix_position`, self, m))
}

Vector3__set_from_matrix_scale <- function(self, m) {
    invisible(.Call(`_ravetools_Vector3__set_from_matrix_scale`, self, m))
}

rawToUInt8 <- function(x) {
    .Call(`_ravetools_rawToUInt8`, x)
}

rawToInt8 <- function(x) {
    .Call(`_ravetools_rawToInt8`, x)
}

rawToUInt16 <- function(x) {
    .Call(`_ravetools_rawToUInt16`, x)
}

rawToInt16 <- function(x) {
    .Call(`_ravetools_rawToInt16`, x)
}

rawToUInt32 <- function(x) {
    .Call(`_ravetools_rawToUInt32`, x)
}

rawToInt32 <- function(x) {
    .Call(`_ravetools_rawToInt32`, x)
}

rawToInt64 <- function(x) {
    .Call(`_ravetools_rawToInt64`, x)
}

rawToFloat <- function(x) {
    .Call(`_ravetools_rawToFloat`, x)
}

rawToString <- function(x) {
    .Call(`_ravetools_rawToString`, x)
}

resample3D <- function(arrayDim, fromArray, newVoxToWorldTransposed, oldVoxToWorldTransposed, na) {
    .Call(`_ravetools_resample3D`, arrayDim, fromArray, newVoxToWorldTransposed, oldVoxToWorldTransposed, na)
}

shiftArray <- function(x, alongIdx, unitIdx, shiftAmount) {
    .Call(`_ravetools_shiftArray`, x, alongIdx, unitIdx, shiftAmount)
}

getDefaultNumThreads <- function() {
    .Call(`_ravetools_getDefaultNumThreads`)
}

vcgIsoSurface <- function(array_, thresh) {
    .Call(`_ravetools_vcgIsoSurface`, array_, thresh)
}

vcgSmoothImplicit <- function(vb_, it_, lambda_, useMassMatrix, fixBorder, useCotWeight, degree, lapWeight_, SmoothQ) {
    .Call(`_ravetools_vcgSmoothImplicit`, vb_, it_, lambda_, useMassMatrix, fixBorder, useCotWeight, degree, lapWeight_, SmoothQ)
}

vcgSmooth <- function(vb_, it_, iter, method, lambda, mu, delta_) {
    .Call(`_ravetools_vcgSmooth`, vb_, it_, iter, method, lambda, mu, delta_)
}

vcgUniformResample <- function(vb_, it_, voxelSize, offsetThr, discretizeFlag, multiSampleFlag, absDistFlag, mergeCloseVert, silent) {
    .Call(`_ravetools_vcgUniformResample`, vb_, it_, voxelSize, offsetThr, discretizeFlag, multiSampleFlag, absDistFlag, mergeCloseVert, silent)
}

vcgUpdateNormals <- function(vb_, it_, select, pointcloud, silent) {
    .Call(`_ravetools_vcgUpdateNormals`, vb_, it_, select, pointcloud, silent)
}

vcgEdgeSubdivision <- function(vb_, it_) {
    .Call(`_ravetools_vcgEdgeSubdivision`, vb_, it_)
}

vcgVolume <- function(mesh_) {
    .Call(`_ravetools_vcgVolume`, mesh_)
}

vcgSphere <- function(subdiv, normals) {
    .Call(`_ravetools_vcgSphere`, subdiv, normals)
}

vcgDijkstra <- function(vb_, it_, source, maxdist_) {
    .Call(`_ravetools_vcgDijkstra`, vb_, it_, source, maxdist_)
}

vcgRaycaster <- function(vb_, it_, rayOrigin, rayDirection, maxDistance, bothSides, threads = 1L) {
    .Call(`_ravetools_vcgRaycaster`, vb_, it_, rayOrigin, rayDirection, maxDistance, bothSides, threads)
}

vcgKDTreeSearch <- function(target_, query_, k, nPointsPerCell = 16L, maxDepth = 64L) {
    .Call(`_ravetools_vcgKDTreeSearch`, target_, query_, k, nPointsPerCell, maxDepth)
}

vcgSubset <- function(vb_, it_, selector_) {
    .Call(`_ravetools_vcgSubset`, vb_, it_, selector_)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call(`_ravetools_RcppExport_registerCCallable`)
})
