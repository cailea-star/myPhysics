# Codes

C++20 computational-physics tools. Each numerical module pairs a header implementation with a concise theory document, while `test/` provides usage, reference output, and acceptance checks.

```text
codes/
├── utils/             # general numerical tools
├── utils_spherical/   # spherical-coordinate tools
├── utils_axial/       # axial-coordinate tools
├── test/              # usage and acceptance tests
├── CMakeLists.txt     # build configuration
├── install.bat        # compiler and dependency installation
├── run.bat            # build and CTest entry point
└── README.md
```

## Commands

Run from the repository root:

```powershell
codes\install.bat   # Install LLVM/Clang, vcpkg, Eigen, and GSL
codes\run.bat       # Configure, build, and run CTest
```

Use [codes-rules](../.codex/skills/codes-rules/SKILL.md) when adding, migrating, reviewing, or modifying numerical tools.
