import os
import re
import subprocess
from pathlib import Path

TEST_FILES = ["test_hfb_cylindrical_01.cpp"]
test_names = [Path(filename).stem for filename in TEST_FILES]
os.chdir(Path(__file__).resolve().parent)
compiler = Path(os.environ["LOCALAPPDATA"]) / "LLVM/bin/clang++.exe"
try:
    subprocess.run(["cmake", "-S", ".", "-B", "build", "-G", "MinGW Makefiles", f"-DCMAKE_CXX_COMPILER={compiler}", "-DCMAKE_CXX_COMPILER_TARGET=x86_64-w64-mingw32", "-DVCPKG_TARGET_TRIPLET=x64-mingw-static"], check=True)
    subprocess.run(["cmake", "--build", "build", "--target", *test_names, "copy_compile_commands"], check=True)
    subprocess.run(["ctest", "--test-dir", "build", "--verbose", "--no-tests=error", "-R", "^(" + "|".join(map(re.escape, test_names)) + ")$" if test_names else "^$"], check=True)
except subprocess.CalledProcessError as error:
    raise SystemExit(error.returncode)
