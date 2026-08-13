using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using ComplexNumber = System.Numerics.Complex;

namespace LSQSolver.Complex
{
    /// <summary>
    /// Provides conversion helpers between complex-valued data and
    /// the realified representations used by LSQSolver.
    /// </summary>
    public static class Adapter
    {
        /// <summary>
        /// Converts a complex vector to an interleaved real vector.
        /// </summary>
        /// <param name="complexVector">The complex vector.</param>
        /// <returns>An array ordered as Re(z0), Im(z0), Re(z1), Im(z1), ...</returns>
        public static double[] AsRealifiedVector(ComplexNumber[] complexVector)
        {
            ArgumentNullException.ThrowIfNull(complexVector);

            int outputLength = checked(2 * complexVector.Length);
            double[] result = new double[outputLength];
            ref double result0 = ref MemoryMarshal.GetArrayDataReference(result);

            for (int i = 0; i < complexVector.Length; i++)
            {
                ComplexNumber value = complexVector[i];
                int offset = 2 * i;

                Unsafe.Add(ref result0, offset) = value.Real;
                Unsafe.Add(ref result0, offset + 1) = value.Imaginary;
            }

            return result;
        }

        /// <summary>
        /// Converts an interleaved real vector to a complex vector.
        /// </summary>
        /// <param name="realifiedVector">The interleaved real vector.</param>
        /// <returns>The reconstructed complex vector.</returns>
        public static ComplexNumber[] AsComplexVector(double[] realifiedVector)
        {
            ArgumentNullException.ThrowIfNull(realifiedVector);

            if ((realifiedVector.Length & 1) != 0)
                throw new ArgumentException("The length of the realified vector must be even.", nameof(realifiedVector));

            int outputLength = realifiedVector.Length / 2;
            ComplexNumber[] result = new ComplexNumber[outputLength];
            ref double input0 = ref MemoryMarshal.GetArrayDataReference(realifiedVector);

            for (int i = 0; i < outputLength; i++)
            {
                int offset = 2 * i;
                result[i] = new ComplexNumber(Unsafe.Add(ref input0, offset), Unsafe.Add(ref input0, offset + 1));
            }

            return result;
        }

        /// <summary>
        /// Converts a column-major complex matrix to its realified matrix.
        /// </summary>
        /// <param name="columnMajorComplexMatrix">The complex matrix in column-major order.</param>
        /// <param name="rows">The number of complex rows.</param>
        /// <param name="cols">The number of complex columns.</param>
        /// <returns>A column-major real matrix with twice as many rows and columns.</returns>
        public static double[] AsColumnMajorRealifiedMatrix(ComplexNumber[] columnMajorComplexMatrix, int rows, int cols)
        {
            ArgumentNullException.ThrowIfNull(columnMajorComplexMatrix);
            ArgumentOutOfRangeException.ThrowIfNegative(rows);
            ArgumentOutOfRangeException.ThrowIfNegative(cols);

            if (columnMajorComplexMatrix.LongLength != (long)rows * cols)
                throw new ArgumentException("The array length must be equal to rows * cols.", nameof(columnMajorComplexMatrix));

            int realifiedRows = checked(2 * rows);
            int realifiedCols = checked(2 * cols);
            int outputLength = checked(realifiedRows * realifiedCols);

            double[] result = new double[outputLength];
            ref double result0 = ref MemoryMarshal.GetArrayDataReference(result);

            for (int j = 0; j < cols; j++)
            {
                int inputColumnOffset = j * rows;
                int realColumnOffset = (2 * j) * realifiedRows;
                int imaginaryColumnOffset = realColumnOffset + realifiedRows;

                for (int i = 0; i < rows; i++)
                {
                    ComplexNumber value = columnMajorComplexMatrix[inputColumnOffset + i];
                    int outputRow = 2 * i;

                    Unsafe.Add(ref result0, realColumnOffset + outputRow) = value.Real;
                    Unsafe.Add(ref result0, realColumnOffset + outputRow + 1) = value.Imaginary;
                    Unsafe.Add(ref result0, imaginaryColumnOffset + outputRow) = -value.Imaginary;
                    Unsafe.Add(ref result0, imaginaryColumnOffset + outputRow + 1) = value.Real;
                }
            }

            return result;
        }

        /// <summary>
        /// Converts a realified column-major matrix back to a complex matrix.
        /// </summary>
        /// <param name="columnMajorRealifiedMatrix">The realified matrix in column-major order.</param>
        /// <param name="realifiedRows">The number of realified rows.</param>
        /// <param name="realifiedCols">The number of realified columns.</param>
        /// <returns>The reconstructed complex matrix in column-major order.</returns>
        public static ComplexNumber[] AsColumnMajorComplexMatrix(double[] columnMajorRealifiedMatrix, int realifiedRows, int realifiedCols)
        {
            ArgumentNullException.ThrowIfNull(columnMajorRealifiedMatrix);
            ArgumentOutOfRangeException.ThrowIfNegative(realifiedRows);
            ArgumentOutOfRangeException.ThrowIfNegative(realifiedCols);

            if (columnMajorRealifiedMatrix.LongLength != (long)realifiedRows * realifiedCols)
                throw new ArgumentException("The array length must be equal to realifiedRows * realifiedCols.", nameof(columnMajorRealifiedMatrix));

            if ((realifiedRows & 1) != 0)
                throw new ArgumentException("The number of realified rows must be even.", nameof(realifiedRows));

            if ((realifiedCols & 1) != 0)
                throw new ArgumentException("The number of realified columns must be even.", nameof(realifiedCols));

            int rows = realifiedRows / 2;
            int cols = realifiedCols / 2;
            int outputLength = checked(rows * cols);

            ComplexNumber[] result = new ComplexNumber[outputLength];
            ref double input0 = ref MemoryMarshal.GetArrayDataReference(columnMajorRealifiedMatrix);

            for (int j = 0; j < cols; j++)
            {
                int realColumnOffset = (2 * j) * realifiedRows;
                int imaginaryColumnOffset = realColumnOffset + realifiedRows;
                int outputColumnOffset = j * rows;

                for (int i = 0; i < rows; i++)
                {
                    int outputRow = 2 * i;

                    double real = Unsafe.Add(ref input0, realColumnOffset + outputRow);
                    double imaginary = Unsafe.Add(ref input0, realColumnOffset + outputRow + 1);
                    double negativeImaginary = Unsafe.Add(ref input0, imaginaryColumnOffset + outputRow);
                    double repeatedReal = Unsafe.Add(ref input0, imaginaryColumnOffset + outputRow + 1);

                    if (negativeImaginary != -imaginary || repeatedReal != real)
                        throw new ArgumentException("The matrix does not have a valid realified complex-matrix structure.", nameof(columnMajorRealifiedMatrix));

                    result[outputColumnOffset + i] = new ComplexNumber(real, imaginary);
                }
            }

            return result;
        }
    }

    /// <summary>
    /// Solves complex least-squares problems by realifying them and
    /// delegating the numerical computation to LSQSolver.
    /// </summary>
    public static class Solver
    {
        /// <summary>
        /// Unit relative rounding error.
        /// </summary>
        const double EPS = 2.22044604925032e-16;

        /// <summary>
        /// Solves the complex least-squares problem Ax ≈ b.
        /// </summary>
        /// <param name="columnMajorMatrix">The complex coefficient matrix in column-major order.</param>
        /// <param name="rows">The number of complex matrix rows.</param>
        /// <param name="cols">The number of complex matrix columns.</param>
        /// <param name="b">The complex right-hand side vector.</param>
        /// <param name="store_intermediates">If true, stores intermediate data for the realified problem.</param>
        /// <param name="rank_tolerance">Relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, validates the input for NaN and Infinity.</param>
        /// <returns>The complex least-squares result.</returns>
        public static ComplexLSQSolverResult Solve(ComplexNumber[] columnMajorMatrix, int rows, int cols, ComplexNumber[] b, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            return Solve(columnMajorMatrix, rows, cols, b, 1, store_intermediates, rank_tolerance, check_finite);
        }

        /// <summary>
        /// Solves the complex multiple-right-hand-side least-squares problem AX ≈ B.
        /// </summary>
        /// <param name="columnMajorMatrix">The complex coefficient matrix in column-major order.</param>
        /// <param name="rows">The number of complex matrix rows.</param>
        /// <param name="cols">The number of complex matrix columns.</param>
        /// <param name="B">The complex right-hand sides stored as a column-major matrix.</param>
        /// <param name="rhs_count">The number of right-hand sides.</param>
        /// <param name="store_intermediates">If true, stores intermediate data for the realified problem.</param>
        /// <param name="rank_tolerance">Relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, validates the input for NaN and Infinity.</param>
        /// <returns>The complex least-squares result.</returns>
        public static ComplexLSQSolverResult Solve(ComplexNumber[] columnMajorMatrix, int rows, int cols, ComplexNumber[] B, int rhs_count, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            LSQSolverStatus status = ValidateSolveInput(columnMajorMatrix, rows, cols, B, rhs_count, check_finite);
            if (status != LSQSolverStatus.Success)
                return ComplexLSQSolverResult.FromValidationFailure(status, rows, cols, rhs_count);

            double[] realifiedMatrix = Adapter.AsColumnMajorRealifiedMatrix(columnMajorMatrix, rows, cols);
            double[] realifiedB = Adapter.AsRealifiedVector(B);

            LSQSolverResult kernelResult = LSQSolver.Solve(
                realifiedMatrix, 2 * rows, 2 * cols, realifiedB, rhs_count,
                overwrite: true,
                store_intermediates: store_intermediates,
                rank_tolerance: rank_tolerance,
                check_finite: false);

            return ComplexLSQSolverResult.FromRealifiedResult(kernelResult);
        }

        /// <summary>
        /// Validates the complex coefficient matrix and right-hand sides before realification.
        /// </summary>
        /// <returns>The validation status for the supplied inputs.</returns>
        private static LSQSolverStatus ValidateSolveInput(ComplexNumber[]? columnMajorMatrix, int rows, int cols, ComplexNumber[]? B, int rhs_count, bool check_finite)
        {
            if (columnMajorMatrix is null) return LSQSolverStatus.NullMatrix;
            if (rows <= 0 || cols <= 0) return LSQSolverStatus.EmptyMatrix;
            if ((long)rows * cols != columnMajorMatrix.LongLength) return LSQSolverStatus.InvalidMatrixStorage;
            if (B is null) return LSQSolverStatus.NullVector;
            if (rhs_count <= 0 || (long)rows * rhs_count != B.LongLength) return LSQSolverStatus.DimensionMismatch;
            if ((long)cols * rhs_count > int.MaxValue / 2) return LSQSolverStatus.DimensionMismatch;

            if (check_finite)
            {
                for (int i = 0; i < B.Length; i++)
                    if (!double.IsFinite(B[i].Real) || !double.IsFinite(B[i].Imaginary))
                        return LSQSolverStatus.InvalidVector;

                for (int i = 0; i < columnMajorMatrix.Length; i++)
                    if (!double.IsFinite(columnMajorMatrix[i].Real) || !double.IsFinite(columnMajorMatrix[i].Imaginary))
                        return LSQSolverStatus.InvalidMatrix;
            }

            return LSQSolverStatus.Success;
        }
    }

    /// <summary>
    /// Represents the result of a complex least-squares solve.
    /// </summary>
    public sealed class ComplexLSQSolverResult
    {
        private readonly LSQSolverStatus status;
        private readonly int rows;
        private readonly int cols;
        private readonly int rhsCount;

        /// <summary>
        /// Gets the status of the solve operation.
        /// </summary>
        public LSQSolverStatus Status => KernelResult?.Status ?? status;

        /// <summary>
        /// Gets the complex solution in column-major order.
        /// </summary>
        public ComplexNumber[] Solution { get; }

        /// <summary>
        /// Gets the number of rows of the original complex coefficient matrix.
        /// </summary>
        public int Rows => KernelResult is null ? rows : KernelResult.Rows / 2;

        /// <summary>
        /// Gets the number of columns of the original complex coefficient matrix.
        /// </summary>
        public int Cols => KernelResult is null ? cols : KernelResult.Cols / 2;

        /// <summary>
        /// Gets the number of right-hand sides.
        /// </summary>
        public int RHSCount => KernelResult?.RHSCount ?? rhsCount;

        /// <summary>
        /// Gets the residual norm of the complex least-squares solution.
        /// </summary>
        public double ResidualNorm => KernelResult?.ResidualNorm ?? 0.0;

        /// <summary>
        /// Gets the result returned by the underlying real-valued LSQSolver.
        /// Intermediate data in this result refers to the realified problem.
        /// </summary>
        public LSQSolverResult? KernelResult { get; }

        /// <summary>
        /// Initializes a result from a completed realified solve.
        /// </summary>
        private ComplexLSQSolverResult(LSQSolverResult kernelResult, ComplexNumber[] solution)
        {
            KernelResult = kernelResult;
            Solution = solution;
        }

        /// <summary>
        /// Initializes a result for an input validation failure.
        /// </summary>
        private ComplexLSQSolverResult(LSQSolverStatus status, int rows, int cols, int rhsCount)
        {
            this.status = status;
            this.rows = rows;
            this.cols = cols;
            this.rhsCount = rhsCount;
            Solution = [];
        }

        /// <summary>
        /// Creates a complex result from the result of a realified least-squares solve.
        /// </summary>
        /// <param name="kernelResult">The result returned by the underlying real-valued LSQSolver.</param>
        /// <returns>The corresponding complex least-squares result.</returns>
        public static ComplexLSQSolverResult FromRealifiedResult(LSQSolverResult kernelResult)
        {
            ArgumentNullException.ThrowIfNull(kernelResult);

            if ((kernelResult.Rows & 1) != 0)
                throw new ArgumentException("The row count of the realified result must be even.", nameof(kernelResult));

            if ((kernelResult.Cols & 1) != 0)
                throw new ArgumentException("The column count of the realified result must be even.", nameof(kernelResult));

            ComplexNumber[] solution = kernelResult.Solution.Length == 0 ? [] : Adapter.AsComplexVector(kernelResult.Solution);
            int expectedSolutionLength = checked((kernelResult.Cols / 2) * kernelResult.RHSCount);

            if (solution.Length != 0 && solution.Length != expectedSolutionLength)
                throw new ArgumentException("The solution length is inconsistent with the dimensions of the realified result.", nameof(kernelResult));

            return new ComplexLSQSolverResult(kernelResult, solution);
        }

        /// <summary>
        /// Creates a complex result for an input validation failure.
        /// </summary>
        internal static ComplexLSQSolverResult FromValidationFailure(LSQSolverStatus status, int rows, int cols, int rhsCount)
        {
            return new ComplexLSQSolverResult(status, rows, cols, rhsCount);
        }
    }
}