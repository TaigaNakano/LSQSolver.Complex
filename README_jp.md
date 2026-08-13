# LSQSolver.Complex

[<u>English version</u>](README.md)

[LSQSolver](https://github.com/TaigaNakano/LSQSolver) の複素数対応拡張です。

`LSQSolver.Complex` は、複素数の密行列最小二乗問題を等価な実数問題へ変換し、数値計算を既存の `LSQSolver` に委譲して解きます。

複素数型には `System.Numerics.Complex` を使用し、単一右辺および複数右辺に対応します。

---

## 概要

母プロジェクトである **LSQSolver** は、列ピボット付きQR分解（CPQR）、数値ランク判定、およびCholesky分解による最小ノルム解の構成を用いた .NET 向けの軽量な密行列最小二乗ソルバーです。過決定、劣決定、ランク落ち問題を、SVDや外部数値計算ライブラリを必須とせずに扱います。

`LSQSolver.Complex` は、その実数値ソルバーを数値計算エンジンとして利用し、複素数の入出力を追加する薄い拡張層です。複素行列と複素ベクトルを実数化し、LSQSolverで解いた後に解を `System.Numerics.Complex` へ戻します。

母プロジェクト: <https://github.com/TaigaNakano/LSQSolver>

---

## インストール

```bash
dotnet add package LSQSolver.Complex
```

`LSQSolver.Complex` は、数値計算カーネルとして元の `LSQSolver` パッケージを利用します。

---

## 使い方

### 複素最小二乗問題を解く

行列は1次元のcolumn-major配列として与えます。

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

上記の `A` はcolumn-major配列なので、行列としては

```text
[ 1+i    -i   ]
[ 2      1+2i ]
```

を表します。`m × n` 行列の要素 `(i, j)` は、配列上では `j * m + i` に格納します。

### `Solve` の引数

単一右辺の場合は次の形です。

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

| 引数 | 説明 |
| --- | --- |
| `columnMajorMatrix` | 係数行列 `A` をcolumn-major順に格納した `Complex[]`。長さは `rows * cols` である必要があります。 |
| `rows` | 元の複素係数行列 `A` の行数。 |
| `cols` | 元の複素係数行列 `A` の列数。 |
| `b` | 複素右辺ベクトル。長さは `rows` である必要があります。 |
| `store_intermediates` | `true` の場合、内部で利用するLSQSolverにQR関連の中間情報を保存させます。これらは実数化問題の中間情報であり、複素QR分解そのものの中間結果ではありません。 |
| `rank_tolerance` | 内部のLSQSolverが数値ランク判定に利用する相対許容値。既定値はLSQSolverで利用している単位丸め誤差です。 |
| `check_finite` | `true` の場合、係数行列と右辺の実部・虚部に `NaN` や `Infinity` が含まれていないかを解法実行前に検査します。 |

Solverへの不正な入力は、母プロジェクトLSQSolverと同様に `ComplexLSQSolverResult.Status` を通じて通知されます。

### 複数右辺

複数右辺もcolumn-major配列で与えます。

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

`rhs_count = k` のとき、`B` は `rows × k` の複素行列をcolumn-major順に格納した配列で、長さは `rows * k` である必要があります。

返される `Solution` は `cols × k` の複素行列をcolumn-major順に格納した配列で、長さは `cols * k` です。

例えば、第 `r` 右辺に対応する解ベクトルは `Solution` の

```text
r * cols
```

から始まります。

`rhs_count` は正の値で、`B` の長さと整合している必要があります。

### 戻り値 `ComplexLSQSolverResult`

`Solver.Solve()` は `ComplexLSQSolverResult` を返します。`Solution` を利用する前に `Status` を確認してください。

| プロパティ | 型 | 説明 |
| --- | --- | --- |
| `Status` | `LSQSolverStatus` | 解法の実行結果。正常終了時は `Success` です。 |
| `Solution` | `Complex[]` | 計算された複素最小二乗解。複数右辺の場合は `Cols × RHSCount` のcolumn-major行列です。 |
| `Rows` | `int` | 元の複素係数行列の行数。 |
| `Cols` | `int` | 元の複素係数行列の列数。 |
| `RHSCount` | `int` | 右辺の本数。 |
| `ResidualNorm` | `double` | 単一右辺では `||Ax-b||₂`、複数右辺では残差行列のFrobeniusノルム `||AX-B||F`。 |
| `KernelResult` | `LSQSolverResult?` | 内部の実数版LSQSolverが返した結果。カーネルを呼び出す前の入力検証で失敗した場合は `null` です。 |

ランク落ち問題や劣決定問題では、解法が成功した場合、母プロジェクトLSQSolverのアルゴリズムに基づく最小2-ノルム解を返します。

`KernelResult` からは、数値ランク、ピボット情報、`RArray`、`Qtb` などの詳細を参照できます。ただし、これらはすべて**実数化後の問題**に対する情報です。ネイティブな複素QR分解におけるランク、ピボット、`R` 因子、変換後右辺として解釈することはできません。

### 変換用Adapter

`Adapter` はSolverとは独立して、複素数表現と実数化表現の変換にも利用できます。

```csharp
double[] realVector = Adapter.AsRealifiedVector(complexVector);
ComplexNumber[] complexVector2 = Adapter.AsComplexVector(realVector);

double[] realMatrix = Adapter.AsColumnMajorRealifiedMatrix(complexMatrix, rows, cols);
ComplexNumber[] complexMatrix2 = Adapter.AsColumnMajorComplexMatrix(realMatrix, 2 * rows, 2 * cols);
```

これらは表現変換用のutilityです。`Solver.Solve()` とは異なり、Adapterへ不正な引数を渡した場合は通常の.NETの引数例外を送出します。

---

## 原理

複素最小二乗問題

```text
A x ≈ b
```

に対し、

```text
A = Ar + i Ai,
x = xr + i xi,
b = br + i bi
```

とすると、等価な実数問題

```text
[ Ar  -Ai ] [ xr ]   [ br ]
[ Ai   Ar ] [ xi ] ≈ [ bi ]
```

へ変換できます。

`LSQSolver.Complex` では、実部と虚部を交互に並べる形式を採用しています。複素行列の要素

```text
a(i,j) = u + iv
```

は、実数行列上では

```text
[  u  -v ]
[  v   u ]
```

という `2 × 2` ブロックとして表されます。

したがって、`m × n` の複素係数行列は `2m × 2n` の実数行列へ変換されます。長さ `m` の複素右辺ベクトルは、

```text
Re(b0), Im(b0), Re(b1), Im(b1), ...
```

という長さ `2m` の実数ベクトルへ変換されます。

複数右辺でも各右辺列を同じ方法で実数化します。右辺の本数そのものは2倍にはなりません。

この実数問題をLSQSolverへ渡し、得られた実数解を再び複素数へ復元します。

実部・虚部を並べた実数ベクトルのEuclideanノルムは通常の複素2-ノルムと一致するため、実数化問題の最小ノルム解は元の複素問題の最小ノルム解に対応します。

---

## 限界

本パッケージは実数版LSQSolverを利用するアダプターであり、複素数専用のQR分解を実装したものではありません。

- `m × n` の複素行列を明示的に `2m × 2n` の実数行列へ展開します。元の複素行列を1要素2個の `double` として保持する場合と比べ、係数行列の生の格納量はおよそ2倍になります。
- `2m × 2n` の実数問題を汎用の実QR分解で解くため、専用の複素QR分解より演算量が増えます。主要項の演算量では、同等の複素数専用QRに対して概ね2倍程度になることが予想されますが、実際の実行時間は実装やハードウェアに依存します。
- `KernelResult.Rank`、`Pivot`、`RArray`、`Qtb` は実数化問題に対する情報です。そのため、複素QR分解のランク、ピボット、R因子などへ変換して公開することはしていません。
- 数値ランク判定も実数化後の問題に対して実数版カーネルが行います。ランク判定閾値付近では、`KernelResult.Rank` は独立に計算した複素ランクではなく、実数表現に対する数値ランク推定値として解釈する必要があります。
- 現在の公開APIはflattenされたcolumn-major配列を利用しており、複素数専用のMatrixObjectは提供していません。
- 入力は新しく確保した `double[]` の作業配列へ明示的に実数化するため、母ソルバーの `overwrite` オプションは複素数入力に対して公開していません。元の `Complex[]` 入力配列は上書きされません。

したがって本拡張は、高度に最適化された複素線形代数ライブラリを置き換えるものではなく、既存のLSQSolverを複素数データから簡単に利用するための拡張として位置付けています。

---

## 今後の対応可能性

今後の開発は、実際の利用状況や必要性に応じて検討します。

候補としては以下があります。

- **複素数用MatrixObject** — 複素問題をより扱いやすくするため、column-major表現を維持した複素行列コンテナ、または軽量なviewを追加する可能性があります。
- **複素数専用エンジン** — 性能、メモリ使用量、あるいは複素数特有の数値的性質が重要になった場合には、明示的な実数化を行わず、複素数または実部・虚部のinterleaved storageを直接扱う専用ソルバーを検討する可能性があります。

現在の実数化方式は、これらの選択肢を将来追加できるよう、できるだけ小さな実装として維持する方針です。

---

## ライセンス

MIT License
