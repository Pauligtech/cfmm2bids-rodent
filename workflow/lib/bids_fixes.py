# workflow/lib/bids_fixes.py

"""
bids_fixes.py — Library of dataset fix functions.
Each fix function operates on a Path and a fix specification dict.
They are auto-registered via the @register_fix decorator.
"""

import json
from pathlib import Path

import nibabel as nib

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


@register_fix("fix_orientation")
def fix_orientation(path: Path, spec: dict) -> bool:
    """Reorient NIfTI file to canonical (RAS+) orientation."""
    if not any(path.name.endswith(ext) for ext in [".nii", ".nii.gz"]):
        return False
    img = nib.load(str(path))
    img = nib.as_closest_canonical(img)
    nib.save(img, str(path))
    return True


def describe_available_fixes():
    """Return a markdown list of all registered fixes and their docstrings."""
    lines = ["### Available Fixes:"]
    for name, func in FIX_REGISTRY.items():
        doc = (func.__doc__ or "").strip().split("\n")[0]
        lines.append(f"- **{name}** — {doc}")
    return "\n".join(lines)
