using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.FEM.Loading.Providers;
using ISAAR.MSolve.FEM.Providers;
using ISSAR.MSolve.Discretization.Loads;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace ISAAR.MSolve.Problems
{
    public class ProblemConvectionDiffusion2 : IConvectionDiffusionIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix> capacity, diffusionConductivityFreeFree, massTransportConductivityFreeFree, stabilizingConductivity,
            firstSpaceDerivativeXMatrix, firstSpaceDerivativeYMatrix, firstSpaceDerivativeZMatrix, 
            secondSpaceDerivativeXMatrix, secondSpaceDerivativeYMatrix, secondSpaceDerivativeZMatrix;
        private Dictionary<int, IVector> stabilizingDomainLoad;
        private Dictionary<int, IMatrixView> diffusionConductivityFreeConstr, diffusionConductivityConstrFree, diffusionConductivityConstrConstr,
            massTransportConductivityFreeConstr, massTransportConductivityConstrFree, massTransportConductivityConstrConstr;
        private readonly Model model;
        private readonly ISolver solver;
        private IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ElementMassTransportConductivityProvider massTransportConductivityProvider = new ElementMassTransportConductivityProvider();
        private ElementDiffusionConductivityProvider diffusionConductivityProvider = new ElementDiffusionConductivityProvider();
        private ElementStructuralStiffnessProvider conductivityProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider capacityProvider = new ElementStructuralMassProvider();
        private ElementStructuralDampingProvider stabilizingConductivityProvider = new ElementStructuralDampingProvider();
        private ElementFirstSpaceDerivativeXProvider firstSpaceDerivativeXProvider = new ElementFirstSpaceDerivativeXProvider();
        private ElementFirstSpaceDerivativeYProvider firstSpaceDerivativeYProvider = new ElementFirstSpaceDerivativeYProvider();
        private ElementFirstSpaceDerivativeZProvider firstSpaceDerivativeZProvider = new ElementFirstSpaceDerivativeZProvider();
        private ElementSecondSpaceDerivativeXProvider secondSpaceDerivativeXProvider = new ElementSecondSpaceDerivativeXProvider();
        private ElementSecondSpaceDerivativeYProvider secondSpaceDerivativeYProvider = new ElementSecondSpaceDerivativeYProvider();
        private ElementSecondSpaceDerivativeZProvider secondSpaceDerivativeZProvider = new ElementSecondSpaceDerivativeZProvider();

        public ProblemConvectionDiffusion2(Model model, ISolver solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(conductivityProvider);
        }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        private IDictionary<int, IMatrix> Capacity
        {
            get
            {
                if (capacity == null) BuildCapacity();
                return capacity;
            }
        }

        private IDictionary<int, IMatrix> DiffusionConductivity
        {
            get
            {
                if (diffusionConductivityFreeFree == null) BuildDiffusionConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return diffusionConductivityFreeFree;
            }
        }
        private IDictionary<int, IMatrix> MassTransportConductivity
        {
            get
            {
                if (massTransportConductivityFreeFree == null) BuildMassTransportConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return massTransportConductivityFreeFree;
            }
        }

        private IDictionary<int, IMatrix> StabilizingConductivity
        {
            get
            {
                if (stabilizingConductivity == null) BuildStabilizingConductivity();
                return stabilizingConductivity;
            }
        }
        private IDictionary<int, IVector> StabilizingDomainLoad
        {
            get
            {
                if (stabilizingDomainLoad == null) stabilizingDomainLoad = model.BuildGlobalStabilizingBodyLoads(solver.DistributeNodalLoads);
                return stabilizingDomainLoad;
            }
        }
        private IDictionary<int, IMatrix> SecondDerivativeXMatrix
        {
            get
            {
                if (secondSpaceDerivativeXMatrix == null) BuildSecondSpaceDerivativeX();
                return secondSpaceDerivativeXMatrix;
            }
        }
        private IDictionary<int, IMatrix> FirstDerivativeXMatrix
        {
            get
            {
                if (firstSpaceDerivativeXMatrix == null) BuildFirstSpaceDerivativeX();
                return firstSpaceDerivativeXMatrix;
            }
        }
        private IDictionary<int, IMatrix> FirstDerivativeYMatrix
        {
            get
            {
                if (firstSpaceDerivativeYMatrix == null) BuildFirstSpaceDerivativeY();
                return firstSpaceDerivativeYMatrix;
            }
        }
        private IDictionary<int, IMatrix> FirstDerivativeZMatrix
        {
            get
            {
                if (firstSpaceDerivativeZMatrix == null) BuildFirstSpaceDerivativeZ();
                return firstSpaceDerivativeZMatrix;
            }
        }
        private IDictionary<int, IMatrix> SecondDerivativeYMatrix
        {
            get
            {
                if (secondSpaceDerivativeYMatrix == null) BuildSecondSpaceDerivativeY();
                return secondSpaceDerivativeYMatrix;
            }
        }
        private IDictionary<int, IMatrix> SecondDerivativeZMatrix
        {
            get
            {
                if (secondSpaceDerivativeZMatrix == null) BuildSecondSpaceDerivativeZ();
                return secondSpaceDerivativeZMatrix;
            }
        }

        private void BuildDiffusionConductivityFreeFree() => diffusionConductivityFreeFree = solver.BuildGlobalMatrices(diffusionConductivityProvider);

        private void BuildDiffusionConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(diffusionConductivityProvider);

            diffusionConductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            diffusionConductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            diffusionConductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            diffusionConductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                diffusionConductivityFreeFree.Add(id, matrices[id].Cff);
                diffusionConductivityFreeConstr.Add(id, matrices[id].Cfc);
                diffusionConductivityConstrFree.Add(id, matrices[id].Ccf);
                diffusionConductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }
        private void BuildMassTransportConductivityFreeFree() => massTransportConductivityFreeFree = solver.BuildGlobalMatrices(massTransportConductivityProvider);

        private void BuildMassTransportConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(massTransportConductivityProvider);

            massTransportConductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            massTransportConductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            massTransportConductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            massTransportConductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                massTransportConductivityFreeFree.Add(id, matrices[id].Cff);
                massTransportConductivityFreeConstr.Add(id, matrices[id].Cfc);
                massTransportConductivityConstrFree.Add(id, matrices[id].Ccf);
                massTransportConductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }

        private void RebuildDiffusionConductivityFreeFree()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (subdomain.StiffnessModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) diffusionConductivityFreeFree = solver.BuildGlobalMatrices(diffusionConductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }
        private void RebuildMassTransportConductivityFreeFree()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (subdomain.StiffnessModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) massTransportConductivityFreeFree = solver.BuildGlobalMatrices(massTransportConductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildFirstSpaceDerivativeX() => firstSpaceDerivativeXMatrix = solver. BuildGlobalMatrices(firstSpaceDerivativeXProvider);
        private void BuildFirstSpaceDerivativeY() => firstSpaceDerivativeYMatrix = solver.BuildGlobalMatrices(firstSpaceDerivativeYProvider);
        private void BuildFirstSpaceDerivativeZ() => firstSpaceDerivativeZMatrix = solver.BuildGlobalMatrices(firstSpaceDerivativeZProvider);
        private void BuildSecondSpaceDerivativeX() => secondSpaceDerivativeXMatrix = solver.BuildGlobalMatrices(secondSpaceDerivativeXProvider);
        private void BuildSecondSpaceDerivativeY() => secondSpaceDerivativeYMatrix = solver.BuildGlobalMatrices(secondSpaceDerivativeYProvider);
        private void BuildSecondSpaceDerivativeZ() => secondSpaceDerivativeZMatrix = solver.BuildGlobalMatrices(secondSpaceDerivativeZProvider);
        private void BuildCapacity() => capacity = solver.BuildGlobalMatrices(capacityProvider);
        private void BuildStabilizingConductivity() => stabilizingConductivity = solver.BuildGlobalMatrices(stabilizingConductivityProvider);
        //private void BuildStabilizingDomainLoad() => stabilizingDomainLoad = solver.DistributeNodalLoads(stabilizingDomainLoadProvider.StabilizingLoad);

        #region IAnalyzerProvider Members
        public void ClearMatrices()
        {
            capacity = null;
            diffusionConductivityFreeFree = null;
            diffusionConductivityFreeConstr = null;
            diffusionConductivityConstrFree = null;
            diffusionConductivityConstrConstr = null;
            massTransportConductivityFreeFree = null;
            massTransportConductivityFreeConstr = null;
            massTransportConductivityConstrFree = null;
            massTransportConductivityConstrConstr = null;
        }

        public void Reset()
        {
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                foreach (IElement element in subdomain.Elements)
                {
                    ((IFiniteElement)element.ElementType).ClearMaterialState();
                }
            }

            diffusionConductivityFreeFree = null;
            diffusionConductivityConstrFree = null;
            diffusionConductivityConstrConstr = null;
            massTransportConductivityFreeFree = null;
            massTransportConductivityConstrFree = null;
            massTransportConductivityConstrConstr = null;
            capacity = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients,
            ISubdomain subdomain)
        {
            #region Backup
            int id = subdomain.ID;
            var conductivityTimesCoeff = DiffusionConductivity[id].Add(MassTransportConductivity[id]);
            conductivityTimesCoeff.ScaleIntoThis(coefficients.Stiffness);

            var capacityTimesCoeff = Capacity[id].Scale(coefficients.Mass);
            capacityTimesCoeff.AddIntoThis(conductivityTimesCoeff);

            return capacityTimesCoeff.LinearCombination(1, StabilizingConductivity[id], coefficients.Damping);
            #endregion
        }

        public void ProcessRhs(ISubdomain subdomain, IVector rhs)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is does not make sense in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is not needed in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.Forces.Clear(); //TODO: this is also done by model.AssignLoads()

            model.AssignLoads(solver.DistributeNodalLoads); // Time-independent nodal loads
            model.AssignTimeDependentNodalLoads(timeStep, solver.DistributeNodalLoads); // Time-dependent nodal loads

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (ISubdomain subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces.Copy());
            return rhsVectors;
        }
        public IVector GetFirstSpaceXDerivatives(ISubdomain subdomain, IVectorView vector)
        {
            var u = (Vector)vector;
            var defX = ((FirstDerivativeXMatrix[subdomain.ID]).GetRow(0)).MultiplyEntrywise(u);
            return defX;
        }
        public IVector GetFirstSpaceYDerivatives(ISubdomain subdomain, IVectorView vector)
        {
            var u = (Vector)vector;
            var defY = ((FirstDerivativeYMatrix[subdomain.ID]).GetRow(0)).MultiplyEntrywise(u);
            return defY;
        }
        public IVector GetFirstSpaceZDerivatives(ISubdomain subdomain, IVectorView vector)
        {
            var u = (Vector)vector;
            var defZ = ((FirstDerivativeZMatrix[subdomain.ID]).GetRow(0)).MultiplyEntrywise(u);
            return defZ;
        }
        public IMatrix GetSecondSpaceDerivatives(ISubdomain subdomain, IVectorView vector)
        {
            var ndofs = model.Subdomains[subdomain.ID].Nodes.Count;
            var v1 = new double[] { 1, 0, 0 };
            var v2 = new double[] { 0, 1, 0 };
            var v3 = new double[] { 0, 0, 1 };
            var i1 = Vector.CreateFromArray(v1);
            var i2 = Vector.CreateFromArray(v2);
            var i3 = Vector.CreateFromArray(v3);
            var secondSpaceDerivatives = Matrix.CreateZero(ndofs, 1);
            var ddX = (Vector)SecondDerivativeXMatrix[subdomain.ID].Multiply(vector);
            var ddY = (Vector)SecondDerivativeYMatrix[subdomain.ID].Multiply(vector);
            var ddZ = (Vector)SecondDerivativeZMatrix[subdomain.ID].Multiply(vector);
            var ddMX = ddX.TensorProduct(i1);
            var ddMY = ddY.TensorProduct(i2);
            var ddMZ = ddZ.TensorProduct(i3);
            ddMX.AddIntoThis(ddMY);
            ddMX.AddIntoThis(ddMZ);
            return ddMX;
        }
        public IVector ConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
        {
            throw new NotImplementedException();
        }
        public IVector SecondDerivativeXMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.SecondDerivativeXMatrix[subdomain.ID].Multiply(vector);
        public IVector SecondDerivativeYMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.SecondDerivativeYMatrix[subdomain.ID].Multiply(vector);
        public IVector SecondDerivativeZMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.SecondDerivativeZMatrix[subdomain.ID].Multiply(vector);
        public IVector CapacityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.Capacity[subdomain.ID].Multiply(vector);
        public IVector DiffusionConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.DiffusionConductivity[subdomain.ID].Multiply(vector);
        public IVector MassTransportConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.MassTransportConductivity[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector StabilizingConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.StabilizingConductivity[subdomain.ID].Multiply(vector);
        public IVector StabilizingRhs(ISubdomain subdomain)
            => this.StabilizingDomainLoad[subdomain.ID];

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain subdomain)
        {
            if (diffusionConductivityFreeFree == null) BuildDiffusionConductivityFreeFree();
            if (massTransportConductivityFreeFree == null) BuildMassTransportConductivityFreeFree();
            return diffusionConductivityFreeFree[subdomain.ID].LinearCombination(1,massTransportConductivityFreeFree[subdomain.ID],1);
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain subdomain)
        {
            int id = subdomain.ID;
            if ((diffusionConductivityFreeFree == null) || (diffusionConductivityFreeConstr == null)
                || (diffusionConductivityConstrFree == null) || (diffusionConductivityConstrConstr == null))
            {
                BuildDiffusionConductivitySubmatrices();
            }
            if ((massTransportConductivityFreeFree == null) || (massTransportConductivityFreeConstr == null)
                || (massTransportConductivityConstrFree == null) || (massTransportConductivityConstrConstr == null))
            {
                BuildMassTransportConductivitySubmatrices();
            }
            var cff = diffusionConductivityFreeFree[id].LinearCombination(1, massTransportConductivityFreeFree[id], 1);
            var cfc = diffusionConductivityFreeConstr[id].LinearCombination(1, massTransportConductivityFreeConstr[id], 1);
            var ccf = diffusionConductivityConstrFree[id].LinearCombination(1, massTransportConductivityConstrFree[id], 1);
            var ccc = diffusionConductivityConstrConstr[id].LinearCombination(1, massTransportConductivityConstrConstr[id], 1);
            return (cff, cfc, ccf, ccc);
        }
        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs) { }

        #endregion
    }
}
