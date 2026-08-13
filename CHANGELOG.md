# Changelog

All notable changes to LSQSolver.Complex are documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers are intended to follow semantic versioning.

## [Unreleased]

### Added

- Added initial complex-valued least-squares support based on `System.Numerics.Complex`.
- Added realification adapters for complex vectors and column-major complex matrices.
  - Complex vectors are converted to interleaved real storage: `Re(z0), Im(z0), Re(z1), Im(z1), ...`.
  - An `m × n` complex matrix is converted to the equivalent `2m × 2n` real block representation.
  - Added reverse conversion helpers for vectors and structurally valid realified matrices.
- Added `Solver` overloads for single and multiple complex right-hand sides.
  - Complex problems are realified and delegated to the real-valued LSQSolver kernel.
  - Temporary realified arrays are passed to the kernel in overwrite mode to avoid unnecessary additional copies.
- Added complex input validation consistent with the status-based error reporting used by LSQSolver.
- Added `ComplexLSQSolverResult`.
  - Exposes the complex solution, original complex dimensions, right-hand-side count, status, and residual norm.
  - Exposes the underlying real-valued `LSQSolverResult` through `KernelResult` for advanced inspection.
  - Keeps QR intermediate data in the realified representation rather than presenting it as native complex QR data.
- Added support for underdetermined and rank-deficient complex problems through the minimum-norm capabilities of the parent LSQSolver engine.
- Added support for multiple right-hand sides through the corresponding LSQSolver kernel API.

### Notes

- The current implementation is an adapter around the real-valued LSQSolver engine rather than a native complex QR implementation.
- Explicit realification increases coefficient-matrix storage and computational work compared with a dedicated complex-valued factorization.
- Possible future extensions include a complex matrix object or view and, if justified by performance or numerical requirements, a dedicated complex numerical engine.
