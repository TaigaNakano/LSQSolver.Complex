# LSQSolver.Complex

[<u>日本語版</u>](https://github.com/TaigaNakano/LSQSolver.Complex/blob/master/README_jp.md)

A lightweight complex-valued extension for [LSQSolver](https://github.com/TaigaNakano/LSQSolver).

`LSQSolver.Complex` solves dense complex least-squares problems by converting them to equivalent real-valued problems and delegating the numerical computation to the original `LSQSolver` engine.

The extension uses `System.Numerics.Complex` and supports single and multiple right-hand sides.

---

## Overview

The parent project, **LSQSolver**, is a lightweight dense least-squares solver for .NET based on column-pivoted QR decomposition (CPQR), numerical rank detection, and Cholesky-based minimum-norm reconstruction. It supports overdetermined, underdetermined, and rank-deficient systems without requiring SVD or an external numerical library.

`LSQSolver.Complex` adds complex-valued input and output while keeping the numerical kernel real-valued. The extension is intentionally thin: complex matrices and vectors are realified, solved by LSQSolver, and converted back to `System.Numerics.Complex`.

Parent project: <https://github.com/TaigaNakano/LSQSolver>

---

## Installation

```bash
dotnet add package LSQSolver.Complex
```

`LSQSolver.Complex` uses the original `LSQSolver` package as its numerical kernel.

---

## Usage

### Solve a complex least-squares problem

Matrices are supplied as one-dimensional arrays in column-major order.

```csharp
using ComplexNumber = System.Numerics.Complex;
using LSQSolver;
using LSQSolver.Complex;

ComplexNumber[] A =
[
    new(1.0, 1.0),
    new(2.0, 0.0),
    new(0.0, -1.0),
    new(1.0, 2.0)
];

ComplexNumber[] b =
[
    new(0.0, 4.0),
    new(-1.0, 3.0)
];

var result = Solver.Solve(A, rows: 2, cols: 2, b);

if (result.Status == LSQSolverStatus.Success)
{
    ComplexNumber[] x = result.Solution;
    Console.WriteLine($"Residual norm: {result.ResidualNorm}");
}
```

The array `A` represents

```text
[ 1+i    -i   ]
[ 2      1+2i ]
```

because matrix elements are stored in column-major order. For an `m × n` matrix, element `(i, j)` is stored at `j * m + i`.

### `Solve` parameters

Single right-hand side:

```csharp
var result = Solver.Solve(
    columnMajorMatrix,
    rows,
    cols,
    b,
    store_intermediates: false,
    rank_tolerance: 2.22044604925032e-16,
    check_finite: true);
```

| Parameter | Description |
| --- | --- |
| `columnMajorMatrix` | Coefficient matrix `A` as a `Complex[]` in column-major order. Its length must be `rows * cols`. |
| `rows` | Number of rows of the original complex matrix `A`. |
| `cols` | Number of columns of the original complex matrix `A`. |
| `b` | Complex right-hand side vector. Its length must be `rows`. |
| `store_intermediates` | If `true`, asks the underlying LSQSolver kernel to retain QR-related intermediate data. These intermediates refer to the realified problem, not to a native complex QR factorization. |
| `rank_tolerance` | Relative tolerance used by the underlying LSQSolver for numerical rank detection. The default is the unit relative rounding error used by LSQSolver. |
| `check_finite` | If `true`, checks the real and imaginary parts of the input matrix and right-hand side for `NaN` or `Infinity` before solving. |

Invalid solver input is reported through `ComplexLSQSolverResult.Status`, following the status-based error contract of the parent LSQSolver project.

### Multiple right-hand sides

Multiple right-hand sides are also stored in column-major order.

```csharp
var result = Solver.Solve(
    A,
    rows,
    cols,
    B,
    rhs_count: 4,
    store_intermediates: false,
    rank_tolerance: 2.22044604925032e-16,
    check_finite: true);
```

For `rhs_count = k`, `B` represents an `rows × k` complex matrix in column-major order and must have length `rows * k`.

The returned `Solution` represents a `cols × k` complex matrix in column-major order and has length `cols * k`.

For example, the entries of the `r`-th solution column begin at

```text
r * cols
```

in `Solution`.

The `rhs_count` parameter must be positive and consistent with the length of `B`.

### Result object

`Solver.Solve()` returns a `ComplexLSQSolverResult`. Check `Status` before using `Solution`.

| Property | Type | Description |
| --- | --- | --- |
| `Status` | `LSQSolverStatus` | Status of the solve operation. On a successful solve, this is `Success`. |
| `Solution` | `Complex[]` | Computed complex solution. For multiple right-hand sides, this is a `Cols × RHSCount` column-major matrix. |
| `Rows` | `int` | Number of rows of the original complex coefficient matrix. |
| `Cols` | `int` | Number of columns of the original complex coefficient matrix. |
| `RHSCount` | `int` | Number of right-hand sides. |
| `ResidualNorm` | `double` | Euclidean norm `||Ax-b||₂` for one right-hand side, or Frobenius norm `||AX-B||F` for multiple right-hand sides. |
| `KernelResult` | `LSQSolverResult?` | Result object produced by the underlying real-valued LSQSolver. It is `null` when validation fails before the kernel is called. |

For rank-deficient or underdetermined problems, a successful solve returns the minimum-2-norm solution inherited from the parent LSQSolver algorithm.

`KernelResult` provides access to lower-level information such as numerical rank, pivot information, `RArray`, and `Qtb`. These values describe the **realified** problem. They must not be interpreted as the rank, pivoting, `R` factor, or transformed right-hand side of a native complex QR factorization.

### Conversion helpers

`Adapter` also exposes conversion functions independently of the solver:

```csharp
double[] realVector = Adapter.AsRealifiedVector(complexVector);
ComplexNumber[] complexVector2 = Adapter.AsComplexVector(realVector);

double[] realMatrix = Adapter.AsColumnMajorRealifiedMatrix(complexMatrix, rows, cols);
ComplexNumber[] complexMatrix2 = Adapter.AsColumnMajorComplexMatrix(realMatrix, 2 * rows, 2 * cols);
```

These methods are representation-conversion utilities. Unlike `Solver.Solve()`, invalid arguments to the adapter methods are reported by normal .NET argument exceptions.

---

## Principle

For a complex least-squares problem

```text
A x ≈ b,
```

write

```text
A = Ar + i Ai,
x = xr + i xi,
b = br + i bi.
```

The problem is equivalent to the real-valued system

```text
[ Ar  -Ai ] [ xr ]   [ br ]
[ Ai   Ar ] [ xi ] ≈ [ bi ].
```

`LSQSolver.Complex` stores the real and imaginary components in an interleaved form. A complex matrix entry

```text
a(i,j) = u + iv
```

is represented by the real `2 × 2` block

```text
[  u  -v ]
[  v   u ].
```

Thus an `m × n` complex coefficient matrix is converted to a `2m × 2n` real matrix. A complex right-hand side with `m` entries is converted to a real vector of length `2m` as

```text
Re(b0), Im(b0), Re(b1), Im(b1), ...
```

For multiple right-hand sides, each complex RHS column is realified in the same way; the number of right-hand sides itself does not double.

The resulting real problem is passed to LSQSolver. The real solution is then reconstructed as a complex solution.

The Euclidean norm of the interleaved real vector is identical to the usual complex 2-norm, so the minimum-norm solution of the realified problem corresponds to the minimum-norm solution of the original complex problem.

---

## Limitations

This package is an adapter around the real-valued LSQSolver kernel, not a native complex QR implementation.

- Explicit realification expands an `m × n` complex matrix into a `2m × 2n` real matrix. Relative to storing the original complex matrix as two doubles per entry, the realified coefficient matrix requires approximately twice the raw coefficient storage.
- A generic real QR factorization of the `2m × 2n` system performs more work than a dedicated complex-valued QR implementation. In leading-order operation counts, the realified approach is expected to be roughly twice as expensive as a comparable native complex QR algorithm, although actual runtime depends strongly on implementation and hardware.
- `KernelResult.Rank`, `Pivot`, `RArray`, and `Qtb` describe the realified real-valued problem. They are intentionally not converted into complex QR intermediates.
- Numerical rank detection is performed by the real-valued kernel on the realified system. Near a numerical rank threshold, its reported rank should therefore be interpreted as the rank estimate of that real representation rather than a separately computed complex rank estimate.
- The current public API uses flattened column-major arrays and does not provide a dedicated complex matrix object.
- Because the input is explicitly realified into newly allocated `double[]` work arrays, this package does not expose the parent solver's `overwrite` option for the original `Complex[]` inputs. The original complex input arrays are not overwritten.

For these reasons, this extension is intended primarily as a simple way to use the existing LSQSolver engine with complex-valued data rather than as a replacement for a highly optimized native complex linear-algebra library.

---

## Possible Future Development

Future development will depend on practical needs and usage.

Possible extensions include:

- **Complex matrix object** — a matrix container or lightweight view may be added to make complex-valued problems easier to construct and inspect while preserving the column-major representation used by the solver.
- **Dedicated complex numerical engine** — if performance, memory use, or complex-specific numerical behavior becomes important, a native complex solver may be considered instead of explicit realification. Such an implementation could operate directly on complex data or interleaved real/imaginary storage while preserving complex column structure throughout the factorization.

The current realification-based implementation is intentionally kept small so that these options can be introduced later without unnecessarily complicating the initial package.

---

## License

MIT License
