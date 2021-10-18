using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using System;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Tests
{
	public class largeTumModel_CCellTransport_only
	{
		//Special inits
		private static double Grox = 7.348597449303836E-07;
		//=============
		private const double timestep = 1;
		private const double time = 30;
		private const int subdomainID = 0;
		private static readonly double[] loxc = new double[] { .07 / 24d / 3600d, 1.0 / 24d / 3600d }; //1/s
		private static readonly double[] Dox = new double[] { 1.78e-9, 1.79e-9 }; //m^2/s
		private static readonly double[] Koxc = new double[] { 0.0083, 0.0083 }; //mol/m^3
		private static readonly double[] Dcell = new double[] { 5.4e-3, 1.8e-4 }; //m^2/s
		private static double cvox = 0.2; //mol/m^3
		private static int solverSymmetric = 0, solverNonSymmetric = 0;
		private static bool reordering = false;
		private static ISolverBuilder builder;
		private static DenseMatrixSolver.Builder asymBuilder;
		private static double[] tumcNode;
		private static double[] tumcElement;
		private static double[] lgElement;
		private static Dictionary<int, double[]> CancerTransportU;
		private static Dictionary<int, double> CancerTransportL;
		private static Dictionary<int, IVector> Displacements;
		private static Tuple<Model, IModelReader> ctModel, structModel;
		private static string inputFile = "mesh.mphtxt";
		private static double MultiModelAnalyzerTolerance = 5e-3;
		static largeTumModel_CCellTransport_only()
		{
			if (solverNonSymmetric == 0)
			{
				asymBuilder = new DenseMatrixSolver.Builder();
			}
			//else
			//{
			//	asymBuilder = new CSparseLUSolver.Builder()
			//	{
			//		DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering())
			//	};
			//}
			if (solverSymmetric == 0)
			{
				builder = new SkylineSolver.Builder();
			}
			else if (solverSymmetric == 1)
			{
				IDofReorderingStrategy reorderingStrategy;
				if (reordering)
				{
					reorderingStrategy = AmdReordering.CreateWithCSparseAmd();
				}
				else
				{
					reorderingStrategy = new NullReordering();
				}
				builder = new CSparseCholeskySolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
			}
			else
			{
				IDofReorderingStrategy reorderingStrategy;
				if (reordering)
				{
					reorderingStrategy = AmdReordering.CreateWithSuiteSparseAmd();
				}
				else
				{
					reorderingStrategy = new NullReordering();
				}
				builder = new SuiteSparseSolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
			}
		}

		[Fact]
		private static void RunTest()
		{
			var path1 = Path.Combine(Directory.GetCurrentDirectory(), $"solutionNorms");
			if (!Directory.Exists(path1))
			{
				Directory.CreateDirectory(path1);
			}
			var path2 = Path.Combine(path1, $"solutionNorm.txt");
			Discretization.Logging.GlobalLogger.OpenOutputFile(path2);

			var DcellDays = new double[Dcell.Length];
			for (int i = 0; i < Dox.Length; i++)
			{
				DcellDays[i] = 24 * 3600 * Dcell[i];
			}

			double[] muLame = new double[] { 6e4, 2.1e4 };
			double[] poissonV = new double[] { .45, .2 };
			IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true), new DynamicMaterial(.001, 0, 0, true) };
			structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);//.Item1; // new Model();

			ctModel = CreateCancerTransportModel(DcellDays[0]);
			var models = new[] { ctModel.Item1 };
			var modelReaders = new[] { ctModel.Item2 };
			IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);
			Discretization.Logging.GlobalLogger.CloseCurrentOutputFile();
			Assert.True(true);
		}
		private static void UpdateModels(Dictionary<int, IVector>[] prevStepSolutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
			IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
		{
			//Tumc
			tumcNode = solversToReplace[0].LinearSystems[0].Solution.CopyToArray();
			if (tumcElement == null) tumcElement = new double[modelsToReplace[0].Elements.Count];
			foreach (var e in modelsToReplace[0].Elements)
			{
				tumcElement[e.ID] = 0; 
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					tumcElement[e.ID] += tumcNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			modelsToReplace[0] = CreateCancerTransportModel(Dcell[0]).Item1;
			asymBuilder.IsMatrixPositiveDefinite = false; //only for models 0, 1, 2
			solversToReplace[0] = asymBuilder.BuildSolver(modelsToReplace[0]); //only for models 0, 1, 2
			providersToReplace[0] = new ProblemConvectionDiffusion2((Model)modelsToReplace[0], solversToReplace[0]);
			childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
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
			string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
			ComsolMeshReader1 modelReader;
			if (lambdag == null)
			{
				modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties);
			}
			else
			{
				modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties, lambdag);
			}
			Model model = modelReader.CreateModelFromFile();
			return new Tuple<Model, IModelReader>(model, modelReader);
		}

		private static void TumorCellsCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			CancerTransportU = u;

			if (tumcElement == null)
			{
				tumcElement = new double[structModel.Item2.elementDomains[0].Count];
				for (int i = 0; i < tumcElement.Count(); i++)
					tumcElement[i] = 0.96;
			}

			foreach (Element e in structModel.Item2.elementDomains[0])
				l[e.ID] = -24d * 3600d * Grox;// * tumcElement[e.ID]; //Rtumc = Grox * Tumc

            CancerTransportL = l;

			//Conservative Flux Convection Coefficient (a)
			//((-mtox * d(cox, x)) + (-WvsTc * mtox * d(Cs, x) * Cs0 )) * stop
			//((-8e-3 * d(cox, x)) + (  -1   * 8e-3 * d(Cs, x) * 1e-3)) *  1
			//cox = cvox * c_ox = 0.2 * c_ox.
			//c_ox considered 0 => cox = 0 everywhere
			//Cs considered 0 everywhere
		}
		private static Tuple<Model, IModelReader> CreateCancerTransportModel(double k)
		{
			ComsolMeshReader4 modelReader;
			Model model;
			if (ctModel == null)
			{
				Console.WriteLine("Creating Cancer Transport Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { k }, TumorCellsCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Cancer Transport Model...");
				modelReader = (ComsolMeshReader4)ctModel.Item2;
				TumorCellsCoefficientsCalculation(CancerTransportU, CancerTransportL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { k }, CancerTransportU, CancerTransportL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}

		private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
		{
			Vector[] initialValues = new Vector[models.Length];
			var value0 = new Dictionary<int, double[]>();
			for (int i = 0; i < models.Length; i++)
			{
				double[] v0;
				v0 = new double[models[i].Nodes.Count];
				value0.Add(i, v0);
				foreach (Node node in models[i].Nodes)
				{
					value0[0][node.ID] = 0.96;//TumC
				}
			}

			IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
			IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
			var solvers = new ISolver[models.Length];
			for (int i = 0; i < models.Length; i++)
			{
				initialValues[i] = Vector.CreateFromArray(value0[i]);
				asymBuilder.IsMatrixPositiveDefinite = false; //only for 0,1,2
				solvers[i] = asymBuilder.BuildSolver(models[i]); //only for 0,1,2
				providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
				childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
			}

            var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
                providers, childAnalyzers, timestep, time, initialTemperature: initialValues, tolerance: MultiModelAnalyzerTolerance);
            parentAnalyzer.Initialize();

			for (int i = 0; i < time / timestep; i++)
			{
				parentAnalyzer.SolveTimestep(i);
			}

			return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
		}
	}
}