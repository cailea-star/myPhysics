@echo off
rem Install LLVM/Clang, vcpkg, and project dependencies.
cd /d "%~dp0"

set "LLVM_ROOT=%LOCALAPPDATA%\LLVM"
set "VCPKG_ROOT=%USERPROFILE%\vcpkg"
set "TRIPLET=x64-mingw-static"

if not exist "%LLVM_ROOT%\bin\clang++.exe" (
    winget install --id LLVM.LLVM -e --location "%LLVM_ROOT%" --accept-package-agreements --accept-source-agreements
    if errorlevel 1 exit /b 1

    if not exist "%LLVM_ROOT%\bin\clang++.exe" (
        echo Error: Clang installation failed at %LLVM_ROOT%
        exit /b 1
    )
)

if not exist "%VCPKG_ROOT%\vcpkg.exe" (
    if not exist "%VCPKG_ROOT%\bootstrap-vcpkg.bat" (
        git clone https://github.com/microsoft/vcpkg.git "%VCPKG_ROOT%"
        if errorlevel 1 exit /b 1
    )

    call "%VCPKG_ROOT%\bootstrap-vcpkg.bat" -disableMetrics
    if errorlevel 1 exit /b 1
)

"%VCPKG_ROOT%\vcpkg.exe" install eigen3 gsl --triplet %TRIPLET% --host-triplet %TRIPLET%
if errorlevel 1 exit /b %errorlevel%
