﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;

//TODO: need special 2x2 Matrix class
//TODO: once we know that an exception will be thrown, try to pinpoint the error: wrong node order, clockwise node order, the  
//      element's shape is too distorted, midpoints are too close to corners in quadratic elements, etc...
namespace ISAAR.MSolve.FEM.Interpolation.Jacobians
{
    /// <summary>
    /// This class encapsulates the determinant and inverse of the Jacobian matrix for a 2D isoparametric mapping.
    /// Let f be a mapping: x \in R^2 -> f(x) \in R^2. The Jacobian matrix of the mapping is (in numerator layout): 
    /// J = [df_1/dx_1 df_1/dx_2; df_2/dx_1 df_2/dx_2]. 
    /// Note that some sources call the transpose of this matrix as J. In FEM we are usually interested in the determinant and
    /// inverse of the Jacobian matrix. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class IsoparametricJacobian2D 
    {
        private const double determinantTolerance = 1E-8; // This needs to be in a static settings class.

        /// <summary>
        /// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
        /// derivatives.
        /// </summary>
        /// <param name="nodes">The nodes used for the interpolation.</param>
        /// <param name="naturalDerivatives">The shape function derivatives at a specific integration point.</param>
        public IsoparametricJacobian2D(IReadOnlyList<Node2D> nodes, Matrix2D naturalDerivatives)
        {
            // The original matrix is not stored. Only the inverse and the determinant
            DirectMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
            (InverseMatrix, DirectDeterminant) = InvertAndDeterminant(DirectMatrix);
            if (DirectDeterminant < determinantTolerance)
            {
                throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
                    + $" ({DirectDeterminant} < {determinantTolerance}). Check the order of nodes or the element geometry.");
            }
        }

        /// <summary>
        /// The determinant of the direct Jacobian matrix <see cref="DirectMatrix"/>.
        /// </summary>
        public double DirectDeterminant { get; }

        /// <summary>
        /// The Jacobian matrix of the direct mapping. Numerator layout is used:
        /// J = [df_1/dx_1 df_1/dx_2; df_2/dx_1 df_2/dx_2].
        /// </summary>
        public Matrix2D DirectMatrix { get; }

        /// <summary>
        /// The inverse of the Jacobian matrix. Numerator layout used is used:
        /// inv(J) = [dx_1/df_1 dx_1/df_2 ; dx_2/df_1 dx_2/df_2]
        /// </summary>
        public Matrix2D InverseMatrix { get; }

        /// <summary>
        /// Transforms the gradient of a vector-valued function from the natural to the global cartesian coordinate system.
        /// </summary>
        /// <param name="naturalGradient">The gradient of a vector-valued function in the natural coordinate system. Each row 
        ///     corresponds to the gradient of a single component of the vector function. Each column corresponds to the 
        ///     derivatives of all components with respect to a single coordinate.</param>
        public Matrix2D TransformNaturalDerivativesToCartesian(Matrix2D naturalGradient) => naturalGradient * InverseMatrix;

        public double[] TransformNaturalDerivativesToCartesian(double[] naturalGradient)
        {
            // naturalGradient * inverseJ = 1-by-2 * 2-by-2 = 1-by-2
            var result = new double[2];
            result[0] = naturalGradient[0] * InverseMatrix[0, 0] + naturalGradient[1] * InverseMatrix[1, 0];
            result[1] = naturalGradient[0] * InverseMatrix[0, 1] + naturalGradient[1] * InverseMatrix[1, 1];
            return result;
        }

        public double[] TransformNaturalDerivativesToCartesian(double derivativeXi, double derivativeEta)
        {
            // naturalGradient * inverseJ = 1-by-2 * 2-by-2 = 1-by-2
            var result = new double[2];
            result[0] = derivativeXi * InverseMatrix[0, 0] + derivativeEta * InverseMatrix[1, 0];
            result[1] = derivativeXi * InverseMatrix[0, 1] + derivativeEta * InverseMatrix[1, 1];
            return result;
        }

        private static Matrix2D CalculateJacobianMatrix(IReadOnlyList<Node2D> nodes, Matrix2D naturalDerivatives)
        {
            //TODO: describe this as a matrix operation
            var J = new double[2, 2];
            for (int nodeIndex = 0; nodeIndex < nodes.Count; ++nodeIndex)
            {
                double x = nodes[nodeIndex].X;
                double y = nodes[nodeIndex].Y;
                double N_xi = naturalDerivatives[nodeIndex, 0];
                double N_eta = naturalDerivatives[nodeIndex, 1];

                J[0, 0] += N_xi * x;
                J[0, 1] += N_eta * x;
                J[1, 0] += N_xi * y;
                J[1, 1] += N_eta * y;
            }
            return new Matrix2D(J);
        }

        private static (Matrix2D inverse, double det) InvertAndDeterminant(Matrix2D directMatrix)
        {
            // Leibniz formula:
            double det = directMatrix[0, 0] * directMatrix[1, 1] - directMatrix[0, 1] * directMatrix[1, 0];
            if (Math.Abs(det) < determinantTolerance) throw new Exception(
                $"|Determinant| = {Math.Abs(det)} < tolerance = {determinantTolerance}. The matrix is singular");

            // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
            double[,] inverse = new double[,] 
            { 
                { directMatrix[1, 1] / det, -directMatrix[0, 1] / det }, 
                { -directMatrix[1, 0] / det, directMatrix[0, 0] / det }
            };

            return (new Matrix2D(inverse), det);
        }
    }
}
