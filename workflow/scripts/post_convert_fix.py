# scripts/curate_bids.py
from pathlib import Path
import shutil, os, json
from lib import bids_fixes

print(bids_fixes.describe_available_fixes())

src = Path(snakemake.input.bids_subj_dir)
dst = Path(snakemake.output.bids_subj_dir)
fixes = snakemake.params.fixes

print(f"Fixing dataset for {src.name} → {dst}")

try:
    shutil.copytree(src, dst, copy_function=os.link)
except Exception:
    print("Hardlink copy failed; using normal copy.")
    shutil.copytree(src, dst)

num_changes = 0
for fix in fixes:
    name = fix.get("name", "unnamed_fix")
    pattern = fix["pattern"]
    action = fix["action"]
    handler = bids_fixes.FIX_REGISTRY.get(action)

    if handler is None:
        raise ValueError(f"Unknown fix type: {action}")

    print(f"\n=== Applying fix: {name} ({action}) ===")
    matches = list(dst.rglob(pattern))
    if not matches:
        print(f"  ⚠️ No matches for pattern: {pattern}")
        continue

    for path in matches:
        if handler(path, fix):
            num_changes += 1

print(f"✅ Done with {src.name}: {num_changes} files modified.")

# Optional lightweight provenance per subject/session
Path(snakemake.output.prov_json).write_text(json.dumps({
    "subject": snakemake.wildcards.subject,
    "session": snakemake.wildcards.session,
    "fixes_used": [f["action"] for f in fixes],
    "files_modified": num_changes,
}, indent=2))



