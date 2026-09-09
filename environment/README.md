# Exact analysis environment

These files record the exact Python/package/runtime environment used for the
analysis:

- Python: 3.12.3
- OS: Linux 6.8.0-138-generic, x86_64, glibc 2.39
- PyTorch runtime: 2.13.0+cu130
- CUDA build: 13.0
- GPU: NVIDIA A2

`pip_freeze.txt` is the complete environment snapshot and includes transitive
CUDA/NVIDIA packages. `requirements.txt` at repository root is the smaller
analysis dependency list for recreating the environment.
