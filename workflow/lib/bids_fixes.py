# workflow/lib/bids_fixes.py

"""
bids_fixes.py — Library of dataset fix functions.
Each fix function operates on a Path and a fix specification dict.
They are auto-registered via the @register_fix decorator.
"""

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

# --- global registry ---
FIX_REGISTRY = {}


def register_fix(name=None):
    """Decorator to register a fix function automatically."""

    def decorator(func):
        fix_name = name or func.__name__
        FIX_REGISTRY[fix_name] = func
        return func

    return decorator


# --- fix implementations ---


@register_fix("remove")
def remove_file(path: Path, spec: dict) -> bool:
    """Remove the file entirely."""
    path.unlink(missing_ok=True)
    return True


@register_fix("update_json")
def update_json(path: Path, spec: dict) -> bool:
    """Update JSON file fields according to `updates` dict."""
    if path.suffix != ".json":
        return False
    updates = spec.get("updates", {})
    with open(path) as f:
        data = json.load(f)
    data.update(updates)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
    return True


def _axcodes2aff(axcodes, scale, translate, labels=None):
    """Create a homogeneous affine from axis codes.

    Uses the provided scale and translate to set diag and offset.

    Parameters
    ----------
    axcodes : sequence of length p
        Axis codes, e.g. ('R','A','S') or (None, 'L', 'S').
    scale: (3,) list of scaling values for X Y Z
    translate: (3,) list of translation values for X Y Z
    labels : sequence of (2,) label tuples, optional
        Same semantics as for axcodes2ornt / ornt2axcodes.  If None, defaults
        to (('L','R'), ('P','A'), ('I','S')).

    Returns
    -------
    aff : (p+1, p+1) ndarray
        Homogeneous affine implementing the permutation and flips implied by
        `axcodes`, with provided translation and scaling.

    Notes
    -----
    - If an axis code is None (a dropped axis), the corresponding column in
      the linear part is left all zeros.
    """
    from nibabel.orientations import axcodes2ornt

    ornt = axcodes2ornt(axcodes, labels)
    p = ornt.shape[0]
    aff = np.zeros((p + 1, p + 1), dtype=float)
    # Fill linear part: for each input axis (column), put a 1 or -1 in the
    # output-axis row indicated by ornt[:,0]
    for in_idx, (out_ax, flip) in enumerate(ornt):
        if np.isnan(out_ax):
            # dropped axis -> leave column zero
            continue
        out_idx = int(out_ax)
        aff[out_idx, in_idx] = float(flip) * scale[in_idx]
        aff[out_idx, p] = translate[in_idx]
    aff[p, p] = 1.0
    return aff


@register_fix("fix_orientation_quadruped")
def fix_orientation_quadruped(path: Path, spec: dict) -> bool:
    """Robust and minimal reorientation of quadruped (sphinx) data."""
    if not any(path.name.endswith(ext) for ext in [".nii", ".nii.gz"]):
        return False

    img = nib.load(path)
    scale = img.header.get_zooms()
    old_affine = img.affine

    # quadruped orientation requires reorder ([0,2,1]) then flip ([1,-1,1])
    quad_ornt = np.array([[0, 2, 1], [1, -1, 1]]).T  #  (e.g. RPI to RSP)

    # apply these transformations to the original orientation
    init_orient = nib.orientations.aff2axcodes(old_affine)

    out_orient_reordered = [init_orient[i] for i in quad_ornt[:, 0]]

    out_orient_flipped = []
    flip_lut = dict(zip("RASLPI", "LPIRAS", strict=False))
    for ax, flip in zip(out_orient_reordered, quad_ornt[:, 1], strict=False):
        if flip == 1:
            out_orient_flipped.append(ax)
        else:
            out_orient_flipped.append(flip_lut[ax])

    out_orient = "".join(out_orient_flipped)

    # get the voxel coordinates of origin (magnet isocentre)
    # using original affine (we want to ensure these same voxel
    # coordinates also get mapped to the magnet isocentre

    origin_ras = np.zeros((4, 1))
    origin_ras[-1, 0] = 1
    origin_old_vox = np.linalg.inv(old_affine) @ origin_ras

    # make an initial affine with zero translation offset
    # (will add the offset later based on magnet isocentre)
    affine = _axcodes2aff(out_orient, scale=scale, translate=np.zeros((3, 1)))

    offset = affine @ origin_old_vox

    # adjust offset to obtain phys origin (ie scanner isocenter) in identical vox location
    affine[:, -1] = -offset[:, 0]

    out_img = nib.Nifti1Image(img.dataobj, affine=affine, header=img.header)
    out_img.to_filename(path)
    return True


def _compute_nifti_hash(path: Path) -> str:
    """Compute hash of NIfTI file data for duplicate detection.

    Loads the NIfTI data and computes MD5 hash of the array data.
    This is more reliable than file-level hashing as it ignores
    minor header differences.
    """
    img = nib.load(path)
    data = np.asanyarray(img.dataobj)
    return hashlib.md5(data.tobytes()).hexdigest()


@register_fix("remove_duplicate_niftis")
def remove_duplicate_niftis(path: Path, spec: dict) -> bool:
    """Remove duplicate NIfTI files keeping the first one (alphanum sorted).

    This fix should be called with a directory path (e.g., pattern: "*" or ".").
    It will find all .nii/.nii.gz files in the directory tree, identify
    duplicates based on image data content, and remove all but the first
    occurrence (alphanumerically sorted). Corresponding .json sidecars are
    also removed.

    Note: This fix uses a class variable to ensure it only runs once even
    if called multiple times (e.g., when pattern matches multiple files).
    """
    # Use a class-level flag to ensure we only run once
    if not hasattr(remove_duplicate_niftis, "_already_run"):
        remove_duplicate_niftis._already_run = set()

    # Get the root directory (subject/session directory)
    if path.is_dir():
        root_dir = path
    else:
        # If called on a file, get the root subject/session directory
        # by going up to find the directory containing this file
        root_dir = path.parent
        # Go up until we find a directory that looks like a BIDS subject/session dir
        while root_dir.parent != root_dir and not (
            any(p.startswith("sub-") for p in root_dir.parts)
        ):
            root_dir = root_dir.parent

    # Convert to string for hashable key
    root_dir_str = str(root_dir)

    # Skip if we've already processed this directory
    if root_dir_str in remove_duplicate_niftis._already_run:
        return False

    remove_duplicate_niftis._already_run.add(root_dir_str)

    # Find all NIfTI files in the directory tree
    nifti_files = []
    for pattern in ["**/*.nii", "**/*.nii.gz"]:
        nifti_files.extend(root_dir.glob(pattern))

    if not nifti_files:
        return False

    # Sort files alphanumerically
    nifti_files = sorted(nifti_files, key=lambda p: str(p))

    # Group files by their content hash
    hash_to_files = {}
    for nifti_file in nifti_files:
        try:
            file_hash = _compute_nifti_hash(nifti_file)
            if file_hash not in hash_to_files:
                hash_to_files[file_hash] = []
            hash_to_files[file_hash].append(nifti_file)
        except Exception as e:
            # If we can't read a file, skip it
            print(f"Warning: Could not compute hash for {nifti_file}: {e}")
            continue

    # Remove duplicates (keep first, remove rest)
    files_removed = 0
    for _file_hash, files in hash_to_files.items():
        if len(files) > 1:
            # Keep the first file (already sorted), remove the rest
            for duplicate_file in files[1:]:
                # Remove the NIfTI file
                duplicate_file.unlink()
                files_removed += 1

                # Remove corresponding JSON sidecar if it exists
                # Handle both .nii and .nii.gz extensions
                if duplicate_file.suffix == ".gz" and str(duplicate_file).endswith(
                    ".nii.gz"
                ):
                    # For .nii.gz, replace with .json
                    json_file = Path(str(duplicate_file)[:-7] + ".json")
                elif duplicate_file.suffix == ".nii":
                    # For .nii, replace with .json
                    json_file = duplicate_file.with_suffix(".json")
                else:
                    # Unknown extension, skip
                    continue

                if json_file.exists():
                    json_file.unlink()
                    files_removed += 1

    return files_removed > 0


def describe_available_fixes():
    """Return a markdown list of all registered fixes and their docstrings."""
    lines = ["### Available Fixes:"]
    for name, func in FIX_REGISTRY.items():
        doc = (func.__doc__ or "").strip().split("\n")[0]
        lines.append(f"- **{name}** — {doc}")
    return "\n".join(lines)
