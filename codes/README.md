# Codes

C++20 computational-physics tools. Each numerical module pairs a header implementation with a concise theory document, while `test/` provides usage, reference output, and acceptance checks.

```text
codes/
├── utils/             # general numerical tools
├── utils_spherical/   # spherical-coordinate tools
├── utils_cylindrical/ # cylindrical-coordinate tools
├── test/              # usage and acceptance tests
├── CMakeLists.txt     # build configuration
├── install.py         # compiler and dependency installation
├── run.py             # selected-test build and CTest entry point
└── README.md
```

## Commands

On Windows, run from the repository root:

```powershell
python codes/install.py   # Install LLVM/Clang, vcpkg, Eigen3, and GSL
python codes/run.py       # Configure, build, and run selected tests
```

Use [codes-rules](../.codex/skills/codes-rules/SKILL.md) when adding, migrating, reviewing, or modifying numerical tools.
