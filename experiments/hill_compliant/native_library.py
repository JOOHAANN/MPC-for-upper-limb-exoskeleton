"""Build immutable native caches, safe for concurrent Windows processes."""

import ctypes
import hashlib
import os
from pathlib import Path
import shutil
import subprocess


def load_library(source, dependencies, root):
    flags = ["-O3", "-std=c99", "-ffp-contract=off", "-shared"]
    flags += ["-static-libgcc"] if os.name == "nt" else ["-fPIC"]
    digest = hashlib.sha256(
        (os.name + " ".join(flags)).encode()
        + b"".join(Path(p).read_bytes() for p in dependencies)
    ).hexdigest()[:20]
    directory = Path(root) / "results/native"
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / (
        Path(source).stem + "_" + digest + (".dll" if os.name == "nt" else ".so")
    )
    if not target.exists():
        compiler = (
            os.environ.get("HILL_GCC")
            or shutil.which("gcc")
            or "C:/Library/CPP/mingw64-posix/bin/gcc.exe"
        )
        temporary = target.with_name(target.stem + f"_{os.getpid()}" + target.suffix)
        subprocess.run(
            [compiler, *flags, str(source), "-o", str(temporary), "-lm"], check=True
        )
        try:
            if not target.exists():
                temporary.replace(target)
        except OSError:
            if not target.exists():
                raise
        finally:
            if temporary.exists():
                temporary.unlink()
    return ctypes.CDLL(str(target))
