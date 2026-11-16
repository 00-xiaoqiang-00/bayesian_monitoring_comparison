"""Resolve output destinations outside the distributed source or submission."""
from pathlib import Path

def protected_root(code_root: Path) -> Path:
    code_root = Path(code_root).expanduser().resolve()
    submission = code_root.parent
    return submission if (submission / '04_Data_and_Results').is_dir() else code_root

def require_external_output(destination: Path, code_root: Path) -> Path:
    destination = Path(destination).expanduser().resolve()
    protected = protected_root(code_root)
    if destination == protected or protected in destination.parents:
        raise ValueError('Choose an output destination outside the submission/code package.')
    return destination
