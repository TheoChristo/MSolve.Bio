
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using System;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISSAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Loading.BodyLoads;
using System.Reflection;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Tests
{
    public class HyperElasticGrowthOxy
    {
        private const int subdomainID = 0;
        private static readonly double[] loxc = new double[] { .07/24/3600, 1.0 / 24 / 3600 }; //1/s
        private static readonly double[] Aox = new double[] { 2200.0/24/3600, 2200.0 / 24 / 3600 }; //mol/(m^3*s)
        private static readonly double[] Dox = new double[] { 1.78e-9, 1.79e-9}; //m^2/s
        private static readonly double[] kox = new double[] { .00464, .00464 }; //mol/m^3
        private static readonly double[] Koxc = new double[] { 0.0083, 0.0083 }; //mol/m^3
        private static double cvox = 0.2; //mol/m^3
        private static double Lwv = 5e-6; //m
        private static double[][] conv0 = new double[][] { new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } };
        //private static double fox = -((Aox * c_ox) / (kox + c_ox * cvox)) * 0.3;
        private static SkylineSolver.Builder builder = new SkylineSolver.Builder();
        //private static DenseMatrixSolver.Builder builder = new DenseMatrixSolver.Builder();
        private static SkylineSolver.Builder structuralBuilder = new SkylineSolver.Builder();
        private static double[] lgNode;
        private static double[] lgElement;
        private static double[] c_oxNode;
        private static double[] c_oxElement;
        //private static IModelReader structuralMR;
        private static Dictionary<int, IVector> Accelerations;
        private static Dictionary<int, IVector> Velocities;
        private static Dictionary<int, IVector> Displacements;
        private static Tuple<Model, IModelReader> oxModel, gModel, structModel;
        [Fact]
        private static void RunTest()
        {
            var DoxDays = new double[Dox.Length];
            for (int i = 0; i < Dox.Length; i++)
            {
                DoxDays[i] = 24 * 3600 * Dox[i];
            }
            oxModel = CreateOxygenTransportModel(DoxDays, conv0, new double[] { Dox[0] / Lwv * 7e3 * 24 * 3600, Dox[1] / Lwv * 7e3 * 24 * 3600 }, c_oxElement);
            gModel = CreateGrowthModel(0, new double[] { 0, 0, 0 }, 0, lgElement);
            var models = new[] { oxModel.Item1, gModel.Item1 };
            var modelReaders = new[] { oxModel.Item2, gModel.Item2 };            
            //var modelTuple3 = CreateStructuralModel(10e4, 0, new DynamicMaterial(.001, 0, 0, true), 0, new double[] { 0, 0, 0 });
            //IVectorView[] solutions = SolveModels(models, modelReaders);
            IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);

            #region paraview commands
            //string path3 = Path.Combine(Directory.GetCurrentDirectory(), "CD_Structural_9HexaOutput.vtu");
            //var numberOfPoints = models[0].Nodes.Count;
            //var numberOfCells = models[0].Elements.Count;
            //using (StreamWriter outputFile = new StreamWriter(path3))
            //{
            //    outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
            //    outputFile.WriteLine("  <UnstructuredGrid>");
            //    outputFile.WriteLine($"     <Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
            //    outputFile.WriteLine("          <Points>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format =\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{models[0].Nodes[i].X} {models[0].Nodes[i].Y} {models[0].Nodes[i].Z} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("          </Points>");
            //    outputFile.WriteLine("          <PointData>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{i + 1}");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"solution1\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Solutions[0][0][i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"solution2\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Solutions[1][0][i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Displacements[0][3 * i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("          </PointData>");
            //    outputFile.WriteLine("          <CellData>");
            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        outputFile.WriteLine($"{i + 1}");
            //    }
            //    outputFile.WriteLine("              </DataArray>");
            //    outputFile.WriteLine("          </CellData>");
            //    outputFile.WriteLine("          <Cells>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"connectivity\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        for (int j = 0; j < models[0].Elements[i].Nodes.Count; j++)
            //            outputFile.Write($"{models[0].Elements[i].Nodes[j].ID} ");
            //        outputFile.WriteLine("");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    var offset = 0;
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        offset += models[0].Elements[i].Nodes.Count;
            //        outputFile.WriteLine($"{offset} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name =\"types\" NumberOfComponents =\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        if (models[0].Elements[i].Nodes.Count == 8)
            //            outputFile.WriteLine($"{12} ");
            //        else outputFile.WriteLine($"{9} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");
            //    outputFile.WriteLine("          </Cells>");
            //    outputFile.WriteLine("      </Piece>");
            //    outputFile.WriteLine("  </UnstructuredGrid>");
            //    outputFile.WriteLine("</VTKFile>");
            //}end*/
            #endregion
            Assert.True(CompareResults(solutions[0]));
        }

        private static void UpdateModels(Dictionary<int, IVector>[] prevStepSolutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
            IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            c_oxNode = solversToReplace[0].LinearSystems[0].Solution.CopyToArray();
            lgNode = solversToReplace[1].LinearSystems[0].Solution.CopyToArray();
            if (c_oxElement == null) c_oxElement = new double[modelsToReplace[0].Elements.Count];
            foreach (var e in modelsToReplace[0].Elements)
            {
                c_oxElement[e.ID] = 0;
                for (int i = 0; i < e.Nodes.Count; i++)
                {
                    c_oxElement[e.ID] += c_oxNode[i] / (e.Nodes.Count);
                }
            }
            if (lgElement == null) lgElement = new double[modelsToReplace[0].Elements.Count];
            foreach (var e in modelsToReplace[1].Elements)
            {
                lgElement[e.ID] = 0;
                for (int i = 0; i < e.Nodes.Count; i++)
                {
                    lgElement[e.ID] += lgNode[e.Nodes[i].ID] / (e.Nodes.Count);
                }
            }
            modelsToReplace[0] = CreateOxygenTransportModel(Dox, conv0, new double[] { Dox[0] / Lwv * 7e3 * 24 * 3600, Dox[1] / Lwv * 7e3 * 24 * 3600 }, c_oxElement).Item1;
            modelsToReplace[1] = CreateGrowthModel(0, new double[] { 0, 0, 0 }, 0, lgElement).Item1;
            for (int i = 0; i < modelsToReplace.Length; i++)
            {
                solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
                providersToReplace[i] = new ProblemConvectionDiffusion2((Model)modelsToReplace[i], solversToReplace[i]);
                childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
            }
        }

        private static void ReplaceLambdaGInModel(IStructuralModel model, double[] lg)
        {
            foreach (var e in model.Elements)
            {
                var et = (ContinuumElement3DNonLinearDefGrad)e.ElementType;
                var bindFlags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static;
                FieldInfo field = typeof(ContinuumElement3DNonLinearDefGrad).GetField("lambdag", bindFlags);
                field.SetValue(et, lg[e.ID]);
            }
        }

        private static void UpdateNewmarkModel(Dictionary<int, IVector> accelerations, Dictionary<int, IVector> velocities, Dictionary<int, IVector> displacements, IStructuralModel[] modelsToReplace,
            ISolver[] solversToReplace, IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            foreach (Element e in structModel.Item2.elementDomains[1])
            { 
                lgElement[e.ID] = 1d; 
            }
            ReplaceLambdaGInModel(modelsToReplace[0], lgElement);
            solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            //solversToReplace[0].HandleMatrixWillBeSet();
            //childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-6;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 5;
            childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
        }

        private static bool CompareResults(IVectorView solution)
        {
            var comparer = new ValueComparer(1E-5);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
            int numFreeDofs = 6;
            if (solution.Length != 6) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }


        private static Tuple<Model, IModelReader> CreateOxygenTransportModel(double[] k, double[][] U, double[] L, double[] coxElement)
        {
            ComsolMeshReader2 modelReader;
            Model model;

            if (oxModel == null)
            {
                Console.WriteLine("Creating Oxygen Model");
                string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "meshXXCoarse.mphtxt");
                modelReader = new ComsolMeshReader2(filename, new double[] { 1, 1 }, k, U, L);
                model = modelReader.CreateModelFromFile();
            }
            else
            {
                Console.WriteLine("Updating Oxygen Model...");
                modelReader = (ComsolMeshReader2)oxModel.Item2;
                modelReader = modelReader.UpdateModelReader(k, U, L);
                model = modelReader.UpdateModel();
            }

            var materials = new ConvectionDiffusionMaterial[] { new ConvectionDiffusionMaterial(1, k[0], U[0], L[0]), new ConvectionDiffusionMaterial(1, k[1], U[1], L[1])};
            if (coxElement == null)
            {
                c_oxElement = new double[model.Elements.Count];
                for (int i = 0; i < model.Elements.Count; i++)
                {
                    c_oxElement[i] = 0;/* 0.9673;*/
                }
            }
            double[] fox = new double[model.Elements.Count];
            int[] domainIDs = new int[] { 0, 1 };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    if (domainID == 0)
                    {
                        fox[element.ID] = (Dox[domainID] / Lwv * 7e3 - ((Aox[domainID] * c_oxElement[element.ID]) / (kox[domainID] + c_oxElement[element.ID] * cvox)) * 0.3) * (24d * 3600d);
                    }
                    else
                    {
                        fox[element.ID] = Dox[domainID] / Lwv * 7e3 * (24d * 3600d);
                    }
                    var nodes = (IReadOnlyList<Node>)element.Nodes;
                    var domainLoad = new ConvectionDiffusionDomainLoad(materials[domainID], fox[element.ID], ThermalDof.Temperature);
                    var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
                    var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
                    model.BodyLoads.Add(bodyLoadElement);
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateGrowthModel(double k, double[] U, double L, double[] lgr)
        {
            ComsolMeshReader3 modelReader;
            Model model;

            if (gModel == null)
            {
                Console.WriteLine("Creating Growth Model");
                string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "meshXXCoarse.mphtxt");
                int[] modelDomains = new int[] { 0 };
                int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
                modelReader = new ComsolMeshReader3(filename, new double[] { 1, 1 }, new double[] { k, k }, new double[][] { U, U }, new double[] { L, 0 });
                model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
            }
            else
            {
                Console.WriteLine("Updating Growth Model...");
                modelReader = (ComsolMeshReader3)gModel.Item2;
                modelReader = modelReader.UpdateModelReader(new double[] { 1, 1 }, new double[] { k, k }, new double[][] { U, U }, new double[] { L, 0 });
                model = modelReader.UpdateModel(structModel.Item1,Displacements, false);
            }

            if (lgr == null)
            {
                lgr = new double[model.Elements.Count];
                for (int i = 0; i < model.Elements.Count; i++)
                {
                    lgr[i] = 1;
                }
            }

            var materialODE = new ConvectionDiffusionMaterial(1, k, U, L);
            double[] Grox = new double[model.Elements.Count];
            double[] fg = new double[model.Elements.Count];
            int[] domainIDs = new int[] { 0, };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    Grox[element.ID] = (loxc[domainID] * cvox * c_oxElement[element.ID]) / (cvox * c_oxElement[element.ID] + Koxc[domainID]);
                    fg[element.ID] = 24d * 3600d * Grox[element.ID] * lgr[element.ID] / 3d;
                    var nodes = (IReadOnlyList<Node>)element.Nodes;
                    var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, fg[element.ID], ThermalDof.Temperature);
                    var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
                    var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
                    model.BodyLoads.Add(bodyLoadElement);
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateStructuralModel(double[] MuLame, double[] PoissonV, IDynamicMaterial[] commonDynamicMaterialProperties, 
            double b, double[] l, double[] lambdag)
        {
            double[] C1 = new double[MuLame.Length];
            double[] C2 = new double[MuLame.Length];
            double[] bulkModulus = new double[MuLame.Length];
            for (int i = 0; i < MuLame.Length; i++)
            {
                //poissonV[i] = 0.2;
                C1[i] = MuLame[i] / 2;
                C2[i] = 0;
                bulkModulus[i] = 2 * MuLame[i] * (1 + PoissonV[i]) / (3 * (1 - 2 * PoissonV[i]));
            }

            ComsolMeshReader1 modelReader;
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "meshXXCoarse.mphtxt");
            if (lambdag == null)
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties);
            else
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties, lambdag);
            Model model = modelReader.CreateModelFromFile();
            
            //Boundary Conditions
            var lx = l[0];
            var ly = l[1];
            var lz = l[2];
            var distributedLoad = new DistributedLoad(lx, ly, lz);


            int[] boundaryIDs = new int[] { 0, 3 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    node.Constraints.Add(new Constraint()
                    {
                        Amount = b,
                        DOF = StructuralDof.TranslationX
                    });
                }
            }
            boundaryIDs = new int[] { 1, 4 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    node.Constraints.Add(new Constraint()
                    {
                        Amount = b,
                        DOF = StructuralDof.TranslationY
                    });
                }
            }
            boundaryIDs = new int[] { 3, 7 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    node.Constraints.Add(new Constraint()
                    {
                        Amount = b,
                        DOF = StructuralDof.TranslationZ
                    });
                }
            }
            boundaryIDs = new int[] { 6 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    model.Loads.Add(new Load() { Node = node, DOF = StructuralDof.TranslationZ, Amount = 1e-4 });
                }
            }
            boundaryIDs = new int[] { 8 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    model.Loads.Add(new Load() { Node = node, DOF = StructuralDof.TranslationY, Amount = 1e-4 });
                }
            }
            boundaryIDs = new int[] { 9 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    model.Loads.Add(new Load() { Node = node, DOF = StructuralDof.TranslationX, Amount = 1e-4 });
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }

        private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
        {
            Vector[] initialValues = new Vector[models.Length];
            var value0 = new Dictionary<int, double[]>();
            for (int i = 0; i < models.Length; i++)
            {
                double[] v0 = new double[models[i].Nodes.Count];
                value0.Add(i, v0);
            }
            foreach (Node node in models[0].Nodes)
            {
                value0[0][node.ID] = 0; /* 0.96733;*/
            }
            foreach (Node node in models[1].Nodes)
            {
                value0[1][node.ID] = 1;
            }

            SkylineSolver[] solvers = new SkylineSolver[models.Length];
            //DenseMatrixSolver[] solvers = new DenseMatrixSolver[models.Length];
            IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                initialValues[i] = Vector.CreateFromArray(value0[i]);
                //var builder = new DenseMatrixSolver.Builder();
                //builder.IsMatrixPositiveDefinite = false;
                solvers[i] = builder.BuildSolver(models[i]);
                providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            const double timestep = 1;
            const double time = 30;
            var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
                providers, childAnalyzers, timestep, time, initialTemperature: initialValues);
            parentAnalyzer.Initialize();
            double[] muLame = new double[] { 6e4, 2.1e4 };
            double[] poissonV = new double[] { .45, .2};
            IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true), new DynamicMaterial(.001, 0, 0, true) };
            structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);//.Item1; // new Model();
            var structuralModel = structModel.Item1;
            var structuralSolver = structuralBuilder.BuildSolver(structuralModel);
            var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
            //var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
            var increments = 2;
            var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
            structuralChildAnalyzerBuilder.ResidualTolerance = 1E-5;
            structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = 5;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
            var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
                structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);
            structuralParentAnalyzer.Initialize();

            for (int i = 0; i < time / timestep; i++)
            {
                parentAnalyzer.SolveTimestep(i);
                structuralParentAnalyzer.SolveTimestep(i);
            }

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }
    }
}