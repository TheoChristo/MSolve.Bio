﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public abstract class Beam2DCorotationalAbstract_v2 : IFiniteElement_v2
    {
        protected static readonly int NATURAL_DEFORMATION_COUNT = 3;
        protected static readonly int FREEDOM_DEGREE_COUNT = 6;
        protected static readonly int AXIS_COUNT = 1;
        protected static readonly int NODE_COUNT = 2;

        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.RotZ };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };

        protected IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();
        protected readonly IFiniteElementMaterial material;
        protected readonly IList<Node_v2> nodes;
        protected readonly double density;
        protected BeamSection2D beamSection;
        protected readonly double initialLength;
        protected double currentLength;
        protected Matrix currentRotationMatrix;
        protected Vector naturalDeformations;
        protected Vector beamAxisX;
        protected Vector beamAxisY;
        protected Vector beamAxisZ;

        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }

        protected Beam2DCorotationalAbstract_v2(IList<Node_v2> nodes, IFiniteElementMaterial material, double density, 
            BeamSection2D beamSection)
        {
            this.nodes = nodes;
            this.material = material;
            this.density = density;
            this.beamSection = beamSection;
            this.initialLength = Math.Sqrt(Math.Pow(nodes[0].X - nodes[1].X, 2) + Math.Pow(nodes[0].Y - nodes[1].Y, 2));
            this.currentLength = this.initialLength;
            this.currentRotationMatrix = Matrix.CreateZero(AXIS_COUNT, AXIS_COUNT);
            this.naturalDeformations = Vector.CreateZero(NATURAL_DEFORMATION_COUNT);
            this.beamAxisX = Vector.CreateZero(AXIS_COUNT);
            this.beamAxisY = Vector.CreateZero(AXIS_COUNT);
        }

        public int ID => 100;
        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
        public bool MaterialModified => material.Modified;
        public IElementDofEnumerator_v2 DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public abstract void SaveGeometryState();
        public abstract void UpdateState(double[] incrementalNodeDisplacements);

        private Matrix CalculateBlockRotationMatrix()
        {
            var blockRotationMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, FREEDOM_DEGREE_COUNT);
            int totalBlocks = 2;
            int blockSize = 3;
            double R11 = this.currentRotationMatrix[0, 0];
            double R12 = this.currentRotationMatrix[0, 1];
            double R13 = this.currentRotationMatrix[0, 2];
            double R21 = this.currentRotationMatrix[1, 0];
            double R22 = this.currentRotationMatrix[1, 1];
            double R23 = this.currentRotationMatrix[1, 2];
            double R31 = this.currentRotationMatrix[2, 0];
            double R32 = this.currentRotationMatrix[2, 1];
            double R33 = this.currentRotationMatrix[2, 2];

            for (int block = 0; block < totalBlocks; block++)
            {
                int cellDistanceCount = block * blockSize;

                blockRotationMatrix[cellDistanceCount, cellDistanceCount] = R11;
                blockRotationMatrix[cellDistanceCount, cellDistanceCount + 1] = R12;
                blockRotationMatrix[cellDistanceCount, cellDistanceCount + 2] = R13;

                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount] = R21;
                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 1] = R22;
                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 2] = R23;

                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount] = R31;
                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 1] = R32;
                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 2] = R33;
            }

            return blockRotationMatrix;
        }

        private Matrix CalculateConstitutiveStiffness()
        {
            var constitutiveStiffness = SymmetricMatrix.CreateZero(FREEDOM_DEGREE_COUNT);
            double E = this.material.YoungModulus;
            double G = E / (2d * (1d + this.material.PoissonRatio));
            double I = this.beamSection.Inertia;
            double A = this.beamSection.Area;
            double L = this.currentLength;
            double LSqruared = L * L;
            double LCubed = L * L * L;
            double phi = (12.0 * E * I) / (LSqruared * G * A);
            double psi = 1.0 / (1.0 + phi);
            double EAOverL = (E * A) / L;
            double EIOverL = (E * I) / L;

            constitutiveStiffness[0, 0] = EAOverL;
            constitutiveStiffness[0, 3] = -EAOverL;

            constitutiveStiffness[1, 1] = 12.0 * psi * E * I / LCubed;
            constitutiveStiffness[1, 2] = 6.0 * psi * E * I / LSqruared;
            constitutiveStiffness[1, 4] = -12.0 * psi * E * I / LCubed;
            constitutiveStiffness[1, 5] = 6.0 * psi * E * I / LSqruared;

            constitutiveStiffness[2, 2] = (3.0 * psi + 1.0) * EIOverL;
            constitutiveStiffness[2, 4] = -6.0 * psi * E * I / LSqruared;
            constitutiveStiffness[2, 5] = (3.0 * psi - 1.0) * EIOverL;

            constitutiveStiffness[3, 3] = EAOverL;

            constitutiveStiffness[4, 4] = 12.0 * psi * E * I / LCubed;
            constitutiveStiffness[4, 5] = -6.0 * psi * E * I / LSqruared;

            constitutiveStiffness[5, 5] = (3.0 * psi + 1.0) * EIOverL;

            return constitutiveStiffness.CopyToFullMatrix();
        }

        private Vector CalculateForcesInGlobalSystem()
        {
            Vector forcesNatural = this.CalculateForcesInNaturalSystem();
            Matrix transformationMatrix = this.CalculateNaturalToLocalTranformMatrix();
            Vector forcesLocal = transformationMatrix * forcesNatural;            
            Matrix blockRotationMatrix = this.CalculateBlockRotationMatrix();
            Vector forcesGlobal = blockRotationMatrix * forcesLocal;
            return forcesGlobal;
        }

        private Vector CalculateForcesInNaturalSystem()
        {
            var forcesNatural = Vector.CreateZero(NATURAL_DEFORMATION_COUNT);
            double E = this.material.YoungModulus;
            double G = E / (2d * (1d + this.material.PoissonRatio));
            double I = this.beamSection.Inertia;
            double A = this.beamSection.Area;
            double L = this.currentLength;
            double phi = (12.0 * E * I) / (L * L * G * A);
            double psi = 1.0 / (1.0 + phi);
            double axialForceNatural = 0d;

            forcesNatural[NaturalDeformationMode2D.EXTENSION] =
                (E * A * this.naturalDeformations[NaturalDeformationMode2D.EXTENSION]) / L;

            axialForceNatural = forcesNatural[NaturalDeformationMode2D.EXTENSION];

            forcesNatural[NaturalDeformationMode2D.SYMMETRIC_BENDING] =
                ((E * I / L) + (L * axialForceNatural / 12.0)) * this.naturalDeformations[NaturalDeformationMode2D.SYMMETRIC_BENDING];
           
            forcesNatural[NaturalDeformationMode2D.ANTISYMMETRIC_BENDING] =
                ((3.0 * psi * E * I/L) + (L*axialForceNatural/20.0)) * this.naturalDeformations[NaturalDeformationMode2D.ANTISYMMETRIC_BENDING];

            return forcesNatural;
        }

        private Matrix CalculateGeometricStiffness()
        {
            double L = this.currentLength;
            var geometricStiffness = SymmetricMatrix.CreateZero(FREEDOM_DEGREE_COUNT);
            var forcesInNaturalSystem = this.CalculateForcesInNaturalSystem();
            double axialForce = forcesInNaturalSystem[NaturalDeformationMode2D.EXTENSION];
            double shearForce = -2.0*forcesInNaturalSystem[2] / L;

            geometricStiffness[0, 1] = -shearForce / L;
            geometricStiffness[0, 4] = shearForce / L;

            geometricStiffness[1, 1] = 6.0 * axialForce / (5.0 * L);
            geometricStiffness[1, 2] = axialForce / 10.0;
            geometricStiffness[1, 3] = shearForce / L;
            geometricStiffness[1, 4] = -6.0 * axialForce / (5.0 * L);
            geometricStiffness[1, 5] = axialForce / 10.0;

            geometricStiffness[2, 2] = 2.0 * axialForce * L / 15.0;
            geometricStiffness[2, 4] = -axialForce / 10.0;
            geometricStiffness[2, 5] = -axialForce * L / 30.0;

            geometricStiffness[3, 4] = -shearForce / L;
            geometricStiffness[4, 4] = 6.0 * axialForce / (5.0 * L);
            geometricStiffness[4, 5] = -axialForce / 10.0;
            geometricStiffness[5, 5] = 2.0 * axialForce * L / 15.0;
            
            return geometricStiffness.CopyToFullMatrix();
        }

        private Matrix CalculateLocalStiffnessMatrix()
        {
            Matrix constitutivePart = this.CalculateConstitutiveStiffness();
            Matrix geometricPart = this.CalculateGeometricStiffness();
            constitutivePart.AddIntoThis(geometricPart);
            return constitutivePart;
        }

        private Matrix CalculateNaturalToLocalTranformMatrix()
        {
            var transformMatrix = Matrix.CreateZero(FREEDOM_DEGREE_COUNT, NATURAL_DEFORMATION_COUNT);
            double L = this.currentLength;

            transformMatrix[0, 0] = -1.0;
            transformMatrix[1, 2] = +2.0 / L;
            transformMatrix[2, 1] = -1.0;
            transformMatrix[2, 2] = +1.0;
            transformMatrix[3, 0] = +1.0;
            transformMatrix[4, 2] = -2.0 / L;
            transformMatrix[5, 1] = +1.0;
            transformMatrix[5, 2] = +1.0;

            return transformMatrix;
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element) => dofTypes;

        public IMatrix StiffnessMatrix(IElement_v2 element)
        {
            Matrix rotationMatrixBlock = this.CalculateBlockRotationMatrix();
            Matrix localStiffnessMatrix = this.CalculateLocalStiffnessMatrix();
            Matrix s = rotationMatrixBlock.MultiplyRight(localStiffnessMatrix).MultiplyRight(rotationMatrixBlock, false, true);
            return s;
        }

        public IMatrix MassMatrix(IElement_v2 element)
        {
            throw new NotImplementedException();
        }

        public IMatrix DampingMatrix(IElement_v2 element)
        {
            throw new NotImplementedException();
            //var m = MassMatrix(element);
            //var lc = m as ILinearlyCombinable;
            //lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            //return m;
        }

        public void ResetMaterialModified() => this.material.ResetModified();

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            UpdateState(localdDisplacements);
            return new Tuple<double[], double[]>(new double[FREEDOM_DEGREE_COUNT], new double[FREEDOM_DEGREE_COUNT]);
        }

        public double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            var internalForces = this.CalculateForcesInGlobalSystem();
            return internalForces.ToRawArray();
        }

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads)
        {
            var accelerations = new double[6];
            IMatrix massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            return massMatrix.Multiply(accelerations);
        }

        public void SaveMaterialState()
        {
            SaveGeometryState();
            material.SaveState();
        }

        public void ClearMaterialState() => material.ClearState();

        public void ClearMaterialStresses() => material.ClearStresses();
    }
}
