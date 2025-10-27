from pathlib import Path
import shutil, os, json
from lib import bids_fixes, utils
import stat

log_file = snakemake.log[0] if snakemake.log else None
logger = utils.setup_logger(log_file)

logger.info(bids_fixes.describe_available_fixes())

src = Path(snakemake.input.bids_subj_dir)
dst = Path(snakemake.output.bids_subj_dir)
fixes = snakemake.params.fixes

logger.info(f"Fixing dataset for {src.name} → {dst}")


# --- Helper functions ---
def make_tree_writable(path: Path):
    """Recursively make files and dirs under `path` writable by owner."""
    for p in [path] + list(path.rglob("*")):
        try:
            mode = p.stat().st_mode
            p.chmod(mode | stat.S_IWUSR)
        except Exception as e:
            logger.info(f"⚠️ Could not make writable: {p} ({e})")

# --- Always do a normal copy ---
logger.info("📂 Copying source tree...")
shutil.copytree(src, dst)  # full, independent copy

# --- Make writable before applying fixes ---
logger.info("🔧 Making copied tree writable...")
make_tree_writable(dst)


num_changes = 0
for fix in fixes:
    name = fix.get("name", "unnamed_fix")
    pattern = fix["pattern"]
    action = fix["action"]
    handler = bids_fixes.FIX_REGISTRY.get(action)

    if handler is None:
        raise ValueError(f"Unknown fix type: {action}")

    logger.info(f"\n=== Applying fix: {name} ({action}) ===")
    matches = list(dst.rglob(pattern))
    if not matches:
        logger.info(f"  ⚠️ No matches for pattern: {pattern}")
        continue

    for path in matches:
        if handler(path, fix):
            num_changes += 1

logger.info(f"✅ Done with {src.name}: {num_changes} files modified.")



# Optional lightweight provenance per subject/session
Path(snakemake.output.prov_json).write_text(json.dumps({
    "subject": snakemake.wildcards.subject,
    "session": snakemake.wildcards.session,
    "fixes_used": [f["action"] for f in fixes],
    "files_modified": num_changes,
}, indent=2))



