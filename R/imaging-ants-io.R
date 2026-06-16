# ANTs / ITK-compatible transform I/O ----------------------------------------
#
# Helpers to read and write the transforms produced by `register_volume3d` in
# the file formats used by 'ANTs' / 'ITK', so they interoperate with
# `antsApplyTransforms`, 'ITK-SNAP', 'FreeSurfer', etc.:
#
#   * the linear part  -> an 'ITK' affine `.mat` (a MATLAB level-4 binary file,
#     LPS convention), hand-written here so the package gains no dependency.
#   * the deformation  -> an 'ANTs' 5-D warp 'NIfTI' (`.nii.gz`, LPS vectors),
#     via the (suggested) 'freesurferformats' package.
#
# Conventions: 'ravetools' works in RAS while 'ITK'/'ANTs' work in LPS. The two
# are related by the symmetric conjugation  M_lps = F %*% M_ras %*% F  with
# F = diag(c(-1, -1, 1, 1)) (F is its own inverse), and for displacement vectors
# by negating the x and y components.


# Convert a 4x4 transform between \verb{RAS} and \verb{LPS}
ras_lps_conjugate <- function(m) {
  m <- as.matrix(m)
  if (nrow(m) == 3L && ncol(m) == 4L) {
    m <- rbind(m, c(0, 0, 0, 1))
  }
  if (!all(dim(m) == c(4L, 4L))) {
    stop("`ras_lps_conjugate`: input must be a 4x4 (or 3x4) matrix.")
  }
  f <- diag(c(-1, -1, 1, 1))
  f %*% m %*% f
}


# ---- MATLAB level-4 (v4) primitives, specialized to 'ITK' transforms --------
# 'ITK' / 'ANTs' write transform `.mat` files in the *level-4* MAT format (not
# level-5): there is no 128-byte file header; each variable is a 20-byte header
# (MOPT type, mrows, ncols, imagf, name length incl. null) followed by the
# null-terminated name and the column-major real data. MOPT = M*1000 + O*100 +
# P*10 + T, where M is the byte order (0 little, 1 big), P the precision (0
# double, 1 single), T the matrix type (0 full). 'ANTs' stores
# `AffineTransform_double_3_3` (12 x 1) and `fixed` (3 x 1), both `double`.

# Is `t` a plausible level-4 MOPT value (M in 0:4, O == 0, P in 0:5, T in 0:2)?
mat4_valid_mopt <- function(t) {
  if (is.na(t) ||
      t < 0L) {
    return(FALSE)
  }
  m <- t %/% 1000L
  rem <- t %% 1000L
  o <- rem %/% 100L
  p <- (rem %/% 10L) %% 10L
  tt <- t %% 10L

  m %in% 0:4 &&
    o == 0L &&
    p %in% 0:5 &&
    tt %in% 0:2
}

# Detect file byte order from the first variable's MOPT word.
mat4_endian <- function(raw) {
  t_le <- readBin(raw[1:4], "integer", n = 1L, size = 4L, endian = "little")
  if (mat4_valid_mopt(t_le)) return("little")
  t_be <- readBin(raw[1:4], "integer", n = 1L, size = 4L, endian = "big")
  if (mat4_valid_mopt(t_be)) return("big")
  "little"
}

# Serialize one real double variable in level-4 format (little-endian).
mat4_write_var <- function(name, values) {
  name_raw <- c(charToRaw(name), as.raw(0L))               # null-terminated
  header <- writeBin(
    as.integer(c(0L, length(values), 1L, 0L, length(name_raw))),
    raw(), size = 4L, endian = "little")                   # MOPT=0, m, n=1, imagf=0, namlen
  c(header, name_raw, writeBin(as.double(values), raw(), size = 8L, endian = "little"))
}

# Parse an entire level-4 buffer into a named list of numeric variables.
mat4_parse <- function(raw) {
  n <- length(raw)
  if (n < 20L)
    stop("not a valid MATLAB level-4 file (too short).")
  endian <- mat4_endian(raw)
  cur <- 1L
  vars <- list()
  while (cur + 19L <= n) {
    hdr <- readBin(raw[cur:(cur + 19L)],
                   "integer",
                   n = 5L,
                   size = 4L,
                   endian = endian)
    type <- hdr[1L]
    mrows <- hdr[2L]
    ncols <- hdr[3L]
    imagf <- hdr[4L]
    namlen <- hdr[5L]
    if (!mat4_valid_mopt(type) ||
        namlen < 1L || mrows < 0L || ncols < 0L) {
      break
    }
    # precision code
    p <- (type %/% 10L) %% 10L
    # drop trailing null
    name <- rawToChar(raw[(cur + 20L):(cur + 20L + namlen - 2L)])
    dstart <- cur + 20L + namlen
    nval <- mrows * ncols

    # double vs single
    sz <- if (p == 0L) {
      8L
    } else {
      4L
    }
    bytes <- nval * sz
    vals <- readBin(raw[dstart:(dstart + bytes - 1L)],
                    "double",
                    n = nval,
                    size = sz,
                    endian = endian)
    vars[[name]] <- vals
    # skip imaginary part if any
    cur <- dstart + bytes * (1L + (imagf == 1L))
  }
  vars
}


#' @title Read and write an \verb{'ITK'}/\verb{'ANTs'} \verb{affine} transform
#' @description
#' \code{write_ants_transform} stores a \eqn{4\times 4} \verb{RAS}-to-\verb{RAS}
#' \verb{affine} (such as the \code{transform} returned by
#' \code{\link{register_volume3d}}) as an \pkg{'ITK'} \code{.mat} file (a
#' \code{'MATLAB'} level-5 binary, \verb{LPS} convention) compatible with
#' \code{antsApplyTransforms}. \code{read_ants_transform} reads such a file back
#' into a \eqn{4\times 4} \verb{RAS} matrix. The reader folds a non-zero
#' \verb{ITK} center of rotation into the translation, so transforms written by
#' \pkg{'ANTs'} are read correctly.
#' @param transform a \eqn{4\times 4} (or \eqn{3\times 4}) \verb{RAS} \verb{affine}
#' @param file path to the \code{.mat} file
#' @returns \code{write_ants_transform} returns \code{file} invisibly;
#' \code{read_ants_transform} returns a \eqn{4\times 4} \verb{RAS} matrix.
#' @seealso \code{\link{register_volume3d}}, \code{\link{save_registration}}
#' @examples
#' tf <- tempfile(fileext = ".mat")
#' m <- diag(4); m[1:3, 4] <- c(3, -2, 1)
#' write_ants_transform(m, tf)
#' read_ants_transform(tf)
#' @export
write_ants_transform <- function(transform, file) {
  m_ras <- as.matrix(transform)
  if (nrow(m_ras) == 3L && ncol(m_ras) == 4L) m_ras <- rbind(m_ras, c(0, 0, 0, 1))
  if (!all(dim(m_ras) == c(4L, 4L))) {
    stop("`write_ants_transform`: `transform` must be a 4x4 (or 3x4) matrix.")
  }
  m_lps <- ras_lps_conjugate(m_ras)
  a <- m_lps[1:3, 1:3]
  tr <- m_lps[1:3, 4]
  # 'ITK' stores the 3x3 matrix row-major, then the translation; center = 0.
  params <- c(as.vector(t(a)), tr)
  out <- c(mat4_write_var("AffineTransform_double_3_3", params),
           mat4_write_var("fixed", c(0, 0, 0)))
  writeBin(out, file)
  invisible(file)
}

#' @rdname write_ants_transform
#' @export
read_ants_transform <- function(file) {
  raw <- readBin(file, "raw", n = file.size(file))
  vars <- mat4_parse(raw)
  anm <- grep("^AffineTransform_(double|float)_[0-9]+_[0-9]+$", names(vars), value = TRUE)
  if (!length(anm)) {
    stop("`read_ants_transform`: no 'AffineTransform_*' variable found in '", file, "'.")
  }
  p <- vars[[anm[1]]]
  if (length(p) < 12L) {
    stop("`read_ants_transform`: '", anm[1], "' has fewer than 12 parameters.")
  }
  a <- matrix(p[1:9], 3L, 3L, byrow = TRUE)
  tr <- p[10:12]
  ctr <- if (!is.null(vars[["fixed"]])) vars[["fixed"]][1:3] else c(0, 0, 0)
  # 'ITK': offset = translation + (I - A) %*% center.
  w <- tr + (diag(3) - a) %*% ctr
  m_lps <- rbind(cbind(a, as.vector(w)), c(0, 0, 0, 1))
  ras_lps_conjugate(m_lps)
}


# ---- hidden-affine codec for the warp 'NIfTI' header ------------------------
# The 4x4 RAS affine is encoded (lossily, opt-in on read) into two standard
# header text fields so a single warp `.nii.gz` can recover it as a last resort.
#   descrip  (80 chars): the RAS 3x3, column-major, 9 x "%+08.5f" comma-joined.
#   aux_file (24 chars): the RAS translation, 3 x "%+07.2f" comma-joined (23),
#                        then a 1-char direction flag "F" (forward) / "I"
#                        (inverse). 'ANTs' ignores both fields.

encode_affine_header <- function(affine, direction = c("forward", "inverse")) {
  direction <- match.arg(direction)
  flag <- if (direction == "forward") "F" else "I"
  m <- as.matrix(affine)
  rot <- sprintf("%+08.5f", as.vector(m[1:3, 1:3]))    # column-major, 9 values
  tr <- sprintf("%+07.2f", m[1:3, 4])
  if (any(nchar(rot) != 8L) || any(nchar(tr) != 7L)) {
    warning("ravetools: affine values exceed the fixed-width header encoding; ",
            "not embedding the affine in the NIfTI header (the .mat keeps full precision).")
    return(list(descrip = "", aux_file = ""))
  }
  list(descrip = paste(rot, collapse = ","),
       aux_file = paste0(paste(tr, collapse = ","), flag))
}

decode_affine_header <- function(descrip, aux_file) {
  descrip <- descrip %||% ""
  aux_file <- aux_file %||% ""
  if (!nzchar(descrip) || !nzchar(aux_file)) return(NULL)
  rot <- tryCatch(scan(text = descrip, sep = ",", quiet = TRUE), error = function(e) numeric(0))
  if (length(rot) != 9L || anyNA(rot)) return(NULL)
  last <- substr(aux_file, nchar(aux_file), nchar(aux_file))
  if (last %in% c("F", "I")) {
    flag <- last
    trtext <- substr(aux_file, 1L, nchar(aux_file) - 1L)
  } else {
    flag <- "F"
    trtext <- aux_file
  }
  tr <- tryCatch(scan(text = trtext, sep = ",", quiet = TRUE), error = function(e) numeric(0))
  if (length(tr) != 3L || anyNA(tr)) return(NULL)
  m <- diag(4)
  m[1:3, 1:3] <- matrix(rot, 3L, 3L)                   # column-major
  m[1:3, 4] <- tr
  list(transform = m, direction = if (identical(flag, "I")) "inverse" else "forward")
}

# Construct a complete NIfTI-1 header list for a 5-D warp field, using only
# exported 'freesurferformats' helpers (no ':::').
build_warp_header <- function(data_dim, vox2ras, descrip = "", aux_file = "") {
  v <- as.matrix(vox2ras)
  sp <- sqrt(colSums(v[1:3, 1:3]^2))                   # voxel spacing
  list(
    endian = "little", sizeof_hdr = 348L,
    dim = freesurferformats::nifti.datadim.to.dimfield(data_dim),
    uses_freesurfer_hack = FALSE,
    intent_p1 = 0, intent_p2 = 0, intent_p3 = 0,
    intent_code = 1007L,                               # NIFTI_INTENT_VECTOR
    datatype = 16L, bitpix = 32L,                      # float32
    slice_start = 0L,
    pix_dim = c(1, sp[1], sp[2], sp[3], 0, 0, 0, 0),
    vox_offset = 352, scl_slope = 0, scl_inter = 0,
    slice_end = 0L, slice_code = 0L,
    xyzt_units = 2L,                                   # NIFTI_UNITS_MM
    cal_max = 0, cal_min = 0, slice_duration = 0, toffset = 0,
    glmax = 0L, glmin = 0L,
    descrip = descrip, aux_file = aux_file,
    qform_code = 0L, sform_code = 1L,                  # use sform (SCANNER_ANAT)
    quatern_b = 0, quatern_c = 0, quatern_d = 0,
    qoffset_x = 0, qoffset_y = 0, qoffset_z = 0,
    srow_x = v[1, ], srow_y = v[2, ], srow_z = v[3, ],
    intent_name = "", magic = "n+1"
  )
}


#' @title Read and write an \verb{'ANTs'} deformation (warp) field
#' @description
#' \code{write_ants_warp} stores a dense displacement field (such as the
#' \code{forward_field} / \code{inverse_field} from
#' \code{\link{register_volume3d}}) as an \pkg{'ANTs'} 5-D warp \code{'NIfTI'}
#' (\code{intent_code} \code{1007}, \verb{LPS} vectors) readable by
#' \code{antsApplyTransforms}; \code{read_ants_warp} reads it back. Requires the
#' suggested \pkg{freesurferformats} package.
#'
#' The field is stored on disk with displacement vectors in \verb{LPS} (the x
#' and y components are negated relative to the \verb{RAS} field) and an
#' \code{sform} equal to \code{vox2ras}, so the file describes the fixed grid's
#' physical space exactly as \pkg{'ANTs'} expects. Optionally the \verb{RAS}
#' \verb{affine} can be hidden in the header text fields
#' (\code{descrip}/\code{aux_file}); \pkg{'ANTs'} ignores those, and recovery on
#' read is opt-in.
#' @param field a \code{(nx, ny, nz, 3)} array of \verb{RAS} displacements; the
#' \code{vox2ras} attribute is used when \code{vox2ras} is not supplied
#' @param file path to the warp \code{.nii} / \code{.nii.gz} file
#' @param vox2ras the fixed-grid \eqn{4\times 4} voxel-to-\verb{RAS} matrix
#' @param affine optional \eqn{4\times 4} \verb{RAS} \verb{affine} to embed (with
#' reduced precision) in the header; \code{NULL} (default) embeds nothing
#' @param direction \code{"forward"} or \code{"inverse"}; recorded as the
#' header direction flag when \code{affine} is embedded
#' @param recover_affine logical; if \code{TRUE}, attempt to decode a hidden
#' \verb{affine} from the header and attach it as the \code{"transform"} attribute
#' (with a \code{"direction"} attribute). Default \code{FALSE}
#' @returns \code{write_ants_warp} returns \code{file} invisibly;
#' \code{read_ants_warp} returns the \code{(nx, ny, nz, 3)} \verb{RAS} field with
#' a \code{"vox2ras"} attribute (and, if recovered, \code{"transform"} /
#' \code{"direction"} attributes).
#' @seealso \code{\link{register_volume3d}}, \code{\link{save_registration}}
#' @export
write_ants_warp <- function(field, file, vox2ras = attr(field, "vox2ras"),
                            affine = NULL, direction = c("forward", "inverse")) {
  if (!requireNamespace("freesurferformats", quietly = TRUE)) {
    stop("`write_ants_warp` requires the 'freesurferformats' package; please install it.")
  }
  direction <- match.arg(direction)
  field <- as.array(field)
  storage.mode(field) <- "double"
  d <- dim(field)
  if (length(d) != 4L || d[4] != 3L) {
    stop("`write_ants_warp`: `field` must be a 4D array with dim (nx, ny, nz, 3).")
  }
  if (is.null(vox2ras)) {
    stop("`write_ants_warp`: `vox2ras` is required (pass it, or set a 'vox2ras' attribute on `field`).")
  }
  v <- as.matrix(vox2ras)
  # RAS -> LPS: negate x and y displacement components.
  field[, , , 1] <- -field[, , , 1]
  field[, , , 2] <- -field[, , , 2]
  # 'ANTs' warps are 5-D: (nx, ny, nz, 1, 3).
  arr <- array(0.0, c(d[1], d[2], d[3], 1L, 3L))
  arr[, , , 1, ] <- field
  txt <- if (is.null(affine)) list(descrip = "", aux_file = "") else
    encode_affine_header(affine, direction)
  h <- build_warp_header(dim(arr), v, txt$descrip, txt$aux_file)
  freesurferformats::write.nifti1(file, arr, niiheader = h)
  invisible(file)
}

#' @rdname write_ants_warp
#' @export
read_ants_warp <- function(file, recover_affine = FALSE) {
  if (!requireNamespace("freesurferformats", quietly = TRUE)) {
    stop("`read_ants_warp` requires the 'freesurferformats' package; please install it.")
  }
  h <- freesurferformats::read.nifti1.header(file)
  arr <- freesurferformats::read.nifti1.data(file, drop_empty_dims = FALSE, header = h)
  storage.mode(arr) <- "double"
  d <- dim(arr)
  if (length(d) == 5L && d[4] == 1L && d[5] == 3L) {
    dim(arr) <- d[c(1L, 2L, 3L, 5L)]
  } else if (!(length(d) == 4L && d[4] == 3L)) {
    stop("`read_ants_warp`: unexpected warp dimensions (", paste(d, collapse = " x "),
         "); expected (nx, ny, nz, 1, 3).")
  }
  # LPS -> RAS: negate x and y displacement components.
  arr[, , , 1] <- -arr[, , , 1]
  arr[, , , 2] <- -arr[, , , 2]
  attr(arr, "vox2ras") <- rbind(h$srow_x, h$srow_y, h$srow_z, c(0, 0, 0, 1))
  if (isTRUE(recover_affine)) {
    dec <- decode_affine_header(h$descrip, h$aux_file)
    if (!is.null(dec)) {
      attr(arr, "transform") <- dec$transform
      attr(arr, "direction") <- dec$direction
    }
  }
  arr
}


# ---- DCF manifest helpers ---------------------------------------------------

mat_to_str <- function(m) {
  m <- as.matrix(m)
  paste(apply(m, 1L, function(r) paste(formatC(r, format = "g", digits = 12L), collapse = " ")),
        collapse = "; ")
}

str_to_mat <- function(s) {
  rows <- strsplit(s, ";", fixed = TRUE)[[1]]
  do.call(rbind, lapply(rows, function(r) as.numeric(strsplit(trimws(r), "\\s+")[[1]])))
}


#' @title Save or load a registration result in \verb{'ANTs'}-compatible files
#' @description
#' \code{save_registration} writes a \code{\link{register_volume3d}} result to
#' disk as \pkg{'ANTs'}-style files: an \verb{'ITK'} \verb{affine} \code{.mat}, the
#' forward/inverse warp \code{'NIfTI'}(s) when present, and a small
#' \verb{'DCF'} manifest (the same key-value format as \code{'DESCRIPTION'})
#' recording what each file is plus the registration parameters.
#' \code{load_registration} reads any of those back: a \code{.mat} returns a
#' \eqn{4\times 4} \verb{RAS} matrix, a warp \code{.nii}/\code{.nii.gz} returns
#' the field, and a \code{.dcf} manifest reassembles the full object.
#' @param x a \code{ravetools_register_volume3d} object (or, for the default
#' method, a \eqn{4\times 4} matrix)
#' @param path output directory, or the manifest (\code{.dcf}) path
#' @param prefix file-name prefix; defaults to \code{"registration"}
#' @param compress logical; write warp fields as \code{.nii.gz} (default) or
#' \code{.nii}
#' @param file a \code{.mat}, \code{.nii}/\code{.nii.gz}, or \code{.dcf} path
#' @param recover_affine_from_header logical; only used as a last resort when
#' loading. If \code{TRUE} and the \verb{affine} \code{.mat} is missing (or a warp
#' file is loaded directly), the \verb{affine} is recovered from the warp
#' header's hidden encoding. Off by default because that encoding is
#' reduced-precision
#' @param ... passed to methods
#' @returns \code{save_registration} returns the manifest path invisibly.
#' \code{load_registration} returns a \eqn{4\times 4} matrix, a field array
#' (with a \code{"vox2ras"} attribute), or a \code{ravetools_register_volume3d}
#' object, depending on the file type.
#' @seealso \code{\link{register_volume3d}}, \code{\link{write_ants_transform}},
#' \code{\link{write_ants_warp}}
#' @export
save_registration <- function(x, path, ...) {
  UseMethod("save_registration")
}

#' @rdname save_registration
#' @export
save_registration.ravetools_register_volume3d <- function(
    x, path, prefix = "registration", compress = TRUE, ...) {

  if (grepl("\\.dcf$", path, ignore.case = TRUE)) {
    manifest <- path
    dir <- dirname(path)
    prefix <- sub("\\.dcf$", "", basename(path), ignore.case = TRUE)
  } else {
    dir <- path
    manifest <- file.path(dir, paste0(prefix, ".dcf"))
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  ext <- if (isTRUE(compress)) ".nii.gz" else ".nii"
  geom <- x$geometry

  affine_file <- paste0(prefix, "0GenericAffine.mat")
  write_ants_transform(x$transform, file.path(dir, affine_file))

  rec <- list(
    Generator = "ravetools",
    ObjectClass = "ravetools_register_volume3d",
    TransformType = x$type %||% "",
    Metric = paste(x$metric, collapse = ","),
    AffineTransform = affine_file)

  if (!is.null(x$forward_field)) {
    fwd <- paste0(prefix, "1Warp", ext)
    write_ants_warp(
      x$forward_field, file.path(dir, fwd),
      vox2ras = attr(x$forward_field, "vox2ras") %||% geom$target_vox2ras,
      affine = x$transform, direction = "forward")
    rec$ForwardField <- fwd
  }
  if (!is.null(x$inverse_field)) {
    inv <- paste0(prefix, "1InverseWarp", ext)
    write_ants_warp(
      x$inverse_field, file.path(dir, inv),
      vox2ras = attr(x$inverse_field, "vox2ras") %||% geom$target_vox2ras,
      affine = solve(x$transform), direction = "inverse")
    rec$InverseField <- inv
  }
  if (!is.null(geom)) {
    if (!is.null(geom$source_vox2ras)) rec$SourceVox2RAS <- mat_to_str(geom$source_vox2ras)
    if (!is.null(geom$target_vox2ras)) rec$TargetVox2RAS <- mat_to_str(geom$target_vox2ras)
    if (!is.null(geom$target_dim))     rec$TargetDim <- paste(geom$target_dim, collapse = " ")
  }

  write.dcf(as.data.frame(rec, stringsAsFactors = FALSE, check.names = FALSE), manifest)
  invisible(manifest)
}

#' @rdname save_registration
#' @export
save_registration.default <- function(x, path, prefix = "registration", ...) {
  if (is.matrix(x) || (is.numeric(x) && length(x) %in% c(12L, 16L))) {
    file <- if (grepl("\\.mat$", path, ignore.case = TRUE)) path else
      file.path(path, paste0(prefix, "0GenericAffine.mat"))
    if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
    write_ants_transform(matrix(x, 4L, 4L), file)
    return(invisible(file))
  }
  stop("`save_registration`: don't know how to save an object of class '",
       paste(class(x), collapse = "/"), "'.")
}

#' @rdname save_registration
#' @export
load_registration <- function(file, recover_affine_from_header = FALSE) {
  lf <- tolower(file)
  if (grepl("\\.mat$", lf)) {
    return(read_ants_transform(file))
  }
  if (grepl("\\.nii(\\.gz)?$", lf)) {
    return(read_ants_warp(file, recover_affine = recover_affine_from_header))
  }
  if (grepl("\\.dcf$", lf)) {
    return(load_registration_manifest(file, recover_affine_from_header))
  }
  stop("`load_registration`: unsupported file '", file,
       "'. Expected a .mat, .nii/.nii.gz, or .dcf manifest.")
}

load_registration_manifest <- function(file, recover_affine_from_header) {
  d <- read.dcf(file)
  rec <- as.list(d[1L, ])
  dir <- dirname(file)

  res <- list(
    type = rec$TransformType,
    metric = if (!is.null(rec$Metric)) strsplit(rec$Metric, ",", fixed = TRUE)[[1]] else NULL)

  af <- if (!is.null(rec$AffineTransform)) file.path(dir, rec$AffineTransform) else NULL
  if (!is.null(af) && file.exists(af)) {
    res$transform <- read_ants_transform(af)
  } else if (isTRUE(recover_affine_from_header) && !is.null(rec$ForwardField)) {
    w <- read_ants_warp(file.path(dir, rec$ForwardField), recover_affine = TRUE)
    res$transform <- attr(w, "transform")
  }

  geom <- list()
  if (!is.null(rec$SourceVox2RAS)) geom$source_vox2ras <- str_to_mat(rec$SourceVox2RAS)
  if (!is.null(rec$TargetVox2RAS)) geom$target_vox2ras <- str_to_mat(rec$TargetVox2RAS)
  if (!is.null(rec$TargetDim)) geom$target_dim <- as.integer(strsplit(trimws(rec$TargetDim), "\\s+")[[1]])
  if (length(geom)) res$geometry <- geom

  if (!is.null(rec$ForwardField)) {
    res$forward_field <- read_ants_warp(file.path(dir, rec$ForwardField))
  }
  if (!is.null(rec$InverseField)) {
    res$inverse_field <- read_ants_warp(file.path(dir, rec$InverseField))
  }

  class(res) <- "ravetools_register_volume3d"
  res
}

#' @export
print.ravetools_register_volume3d <- function(x, ...) {
  cat("<ravetools registration>\n")
  cat("  type:   ", x$type %||% "?", "\n")
  cat("  metric: ", paste(x$metric, collapse = ", "), "\n")
  has <- c(
    if (!is.null(x$transform)) "affine",
    if (!is.null(x$forward_field)) "forward_field",
    if (!is.null(x$inverse_field)) "inverse_field")
  cat("  parts:  ", paste(has, collapse = ", "), "\n")
  if (!is.null(x$geometry$target_dim)) {
    cat("  grid:   ", paste(x$geometry$target_dim, collapse = " x "), "\n")
  }
  invisible(x)
}
