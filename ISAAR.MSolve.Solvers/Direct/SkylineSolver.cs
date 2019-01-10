﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: factorization should be done during Solve() to time correctly. Alternatively, an observer should record the durations.
//TODO: investigate if it is possible to avoid casting the matrix provided by the analyzer/assembler into skyline. Perhaps the
//      the matrix could be obtained directly from the assembler and the analyzer could provide delegates with the operations it 
//      wants done on the matrix, instead of doing them itself.
//TODO: try to abstract the subdomain logic from ther analyzers. I doubt it is possible though.
//TODO: directly pass the single linear system instead of a list that must be checked. The same holds for all solvers and 
//      assemblers.
namespace ISAAR.MSolve.Solvers.Direct
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on sparse symmetric positive definite 
    /// matrices stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver : SingleSubdomainSolverBase<SkylineMatrix>
    {
        private readonly double factorizationPivotTolerance;

        private bool factorizeInPlace = true;
        private bool mustFactorize = true;
        private CholeskySkyline factorizedMatrix;

        private SkylineSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer):
            base(model, dofOrderer, new SkylineAssembler(), "SkylineSolver")
        {
            this.factorizationPivotTolerance = factorizationPivotTolerance;
        }

        public override void Initialize() { }

        public override void OnMatrixSetting()
        {
            mustFactorize = true;
            factorizedMatrix = null;
        }

        public override void PreventFromOverwrittingMatrix() => factorizeInPlace = false;

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public override void Solve()
        {
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            if (mustFactorize)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(factorizeInPlace, factorizationPivotTolerance); 
                mustFactorize = false;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        //TODO: Copied from Stavroulakis code. Find out what the purpose of this is. I suspect he wanted to compare with some 
        //      old solution that used single precision.
        //private void DestroyAccuracy(ILinearSystem_v2 linearSystem)
        //{
        //    if (AccuracyDigits < 1) return;

        //    for (int i = 0; i < linearSystem.RhsVector.Length; i++)
        //    {
        //        //ScientificDouble s = ScientificDouble.GetScientificDouble(subdomain.RHS[i]);
        //        //s.ReduceAccuracy(AccuracyDigits);
        //        //linearSystem.RhsVector[i] = ScientificDouble.GetDouble(s);
        //        linearSystem.RhsVector[i] = Double.Parse(String.Format("{0:" + stringFormat + "}", linearSystem.RhsVector[i]));
        //    }
        //}

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } 
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            public SkylineSolver BuildSolver(IStructuralModel_v2 model)
            {
                return new SkylineSolver(model, FactorizationPivotTolerance, DofOrderer);
            }
        }
    }
}
