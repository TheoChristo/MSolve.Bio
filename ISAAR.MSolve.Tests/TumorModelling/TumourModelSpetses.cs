using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using System;
using Xunit;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.Tests
{
	public class largeTumModel_Full_spetses
	{
		private const double timestep = 1;
		private const double time = 30;
        //Solver
		private static ISolverBuilder builder, asymBuilder, structuralBuilder;

		private static string inputFile = "Code_new2Mesh.mphtxt";

		static largeTumModel_Full_spetses()
        {
            //Solver config
            asymBuilder = new DenseMatrixSolver.Builder();
            builder = new SkylineSolver.Builder();
            structuralBuilder = new SkylineSolver.Builder();

            //Output config
            var path1 = Path.Combine(Directory.GetCurrentDirectory(), $"solutionNorms");
            if (!Directory.Exists(path1)) Directory.CreateDirectory(path1);
            Discretization.Logging.GlobalLogger.OpenOutputFile(Path.Combine(path1, $"solutionNorm.txt"));


            //Init models





            //IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);

            //Output
            Discretization.Logging.GlobalLogger.CloseCurrentOutputFile();
            //Compare Results 
        }

		[Fact]
		private static void RunTest()
        {


        }

        private static void CreateStructuralModel()
        {
            //double[] muLame = new double[] { 6e4, 2.1e4 };
            double[] poissonV = new double[] { .45, .2 };
            IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true), new DynamicMaterial(.001, 0, 0, true) };
            structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);
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
        }

        private static bool CompareResults(IVectorView solution)
        {
            return true;
        }

        private static void Paraview(int timeStep)
        {
            //var path0 = Path.Combine(Directory.GetCurrentDirectory(), $"paraviewOutput");
            //var path3 = Path.Combine(path0, $"results{timeStep}.vtu");
            //var numberOfPoints = structModel.Item1.Nodes.Count;
            //var numberOfCells = structModel.Item1.Elements.Count;
            //using (StreamWriter outputFile = new StreamWriter(path3))
            //{
            //    outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
            //    outputFile.WriteLine("  <UnstructuredGrid>");
            //    outputFile.WriteLine($"     <Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
            //    outputFile.WriteLine("          <Points>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format =\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{structModel.Item1.Nodes[i].X} {structModel.Item1.Nodes[i].Y} {structModel.Item1.Nodes[i].Z} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("          </Points>");
            //    outputFile.WriteLine("          <PointData>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{i + 1}");
            //    outputFile.WriteLine("              </DataArray>");


            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"totalDisplacement\" NumberOfComponents=\"1\" format =\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //    {
            //        double dist = Math.Sqrt(Math.Pow(oxModel.Item1.Nodes[i].X - structModel.Item1.Nodes[i].X, 2) +
            //            Math.Pow(oxModel.Item1.Nodes[i].Y - structModel.Item1.Nodes[i].Y, 2) +
            //            Math.Pow(oxModel.Item1.Nodes[i].Z - structModel.Item1.Nodes[i].Z, 2));
            //        outputFile.WriteLine($"{dist} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"oxygen\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{c_oxNode[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"growth factor\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < lgNode.Length; i++)
            //        outputFile.WriteLine($"{lgNode[i]} ");
            //    for (int i = lgNode.Length; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"1 ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"cancer cells\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < tumcNode.Length; i++)
            //        outputFile.WriteLine($"{tumcNode[i]} ");
            //    for (int i = tumcNode.Length; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"0 ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{pNode[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Cs\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < CsNode.Length; i++)
            //        outputFile.WriteLine($"{CsNode[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"SvD\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < SvDNode.Length; i++)
            //        outputFile.WriteLine($"{SvDNode[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Cv\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < CvNode.Length; i++)
            //        outputFile.WriteLine($"{CvNode[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Ang1\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < a1Node.Length; i++)
            //        outputFile.WriteLine($"{a1Node[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Ang2\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < a2Node.Length; i++)
            //        outputFile.WriteLine($"{a2Node[i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"phis\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < phisNode.Length; i++)
            //        outputFile.WriteLine($"{phisNode[i]} ");
            //    for (int i = phisNode.Length; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"0.3 ");
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
            //        for (int j = 0; j < structModel.Item1.Elements[i].Nodes.Count; j++)
            //            outputFile.Write($"{structModel.Item1.Elements[i].Nodes[j].ID} ");
            //        outputFile.WriteLine("");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    var offset = 0;
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        offset += structModel.Item1.Elements[i].Nodes.Count;
            //        outputFile.WriteLine($"{offset} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name =\"types\" NumberOfComponents =\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        if (structModel.Item1.Elements[i].Nodes.Count == 4)
            //            outputFile.WriteLine($"{10} ");
            //        else outputFile.WriteLine($"{5} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");
            //    outputFile.WriteLine("          </Cells>");
            //    outputFile.WriteLine("      </Piece>");
            //    outputFile.WriteLine("  </UnstructuredGrid>");
            //    outputFile.WriteLine("</VTKFile>");
            //}
        }

		private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
        {
            

            return new Vector[12];
        }
    }
}