import os
import subprocess
from pathlib import Path

llvm_root = Path(os.environ["LOCALAPPDATA"]) / "LLVM"
vcpkg_root = Path(os.environ["USERPROFILE"]) / "vcpkg"
os.chdir(Path(__file__).resolve().parent)
if not (llvm_root / "bin/clang++.exe").exists():
    subprocess.run(["winget", "install", "--id", "LLVM.LLVM", "-e", "--location", str(llvm_root), "--accept-package-agreements", "--accept-source-agreements"], check=True)
if not (llvm_root / "bin/clang++.exe").exists(): raise SystemExit(f"Clang installation failed at {llvm_root}")
if not (vcpkg_root / "vcpkg.exe").exists():
    if not (vcpkg_root / "bootstrap-vcpkg.bat").exists(): subprocess.run(["git", "clone", "https://github.com/microsoft/vcpkg.git", str(vcpkg_root)], check=True)
    subprocess.run(["cmd", "/c", str(vcpkg_root / "bootstrap-vcpkg.bat"), "-disableMetrics"], check=True)
subprocess.run([str(vcpkg_root / "vcpkg.exe"), "install", "eigen3", "gsl", "--triplet", "x64-mingw-static", "--host-triplet", "x64-mingw-static"], check=True)
