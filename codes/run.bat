@echo off
cd /d "%~dp0"

cmake -S . -B build -G "MinGW Makefiles" ^
    -DCMAKE_CXX_COMPILER="%LOCALAPPDATA%\LLVM\bin\clang++.exe" ^
    -DCMAKE_CXX_COMPILER_TARGET=x86_64-w64-mingw32 ^
    -DVCPKG_TARGET_TRIPLET=x64-mingw-static
if errorlevel 1 exit /b %errorlevel%

cmake --build build
if errorlevel 1 exit /b %errorlevel%

ctest --test-dir build --verbose
exit /b %errorlevel%
