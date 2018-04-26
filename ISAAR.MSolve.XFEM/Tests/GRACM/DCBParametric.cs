﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Solvers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class DCBParametric
    {
        public static void Main()
        {
            //SingleTest();
            CompareSolvers();
            //ParametricCrackLength(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
            //ParametricMesh(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
            //GridSearch(jIntegralRadiusOverElementSize, fractureToughness, maxIterations);
        }

        private static void SingleTest()
        {
            double growthLength = 0.3;
            double elementSize = 0.08;
            var benchmark = new DCB(elementSize, growthLength);
            benchmark.UniformMesh = false;
            benchmark.UseLSM = true;
            benchmark.InitializeModel();
            Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                elementSize, benchmark.Model.Elements.Count, growthLength);

            //Print
            //VTKWriter writer = new VTKWriter(benchmark.model);
            //writer.InitializeFile("dcb_transfinite");
            //writer.CloseCurrentFile();

            var solver = new SkylineSolver();
            IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);
            Console.WriteLine("Crack path:");
            foreach (var point in crackPath)
            {
                Console.WriteLine("{0} {1}", point.X, point.Y);
            }
        }

        private static void CompareSolvers()
        {
            double growthLength = 0.3;
            double elementSize = 0.08;
            int repetitions = 10;

            var solvers = new Dictionary<string, ISolver>
            {
                {  "Skyline", new SkylineSolver()},
                {  "Jacobi Preconditioned CG", new PCGSolver(1, 1e-8)}
            };

            var solverTimes = new Dictionary<string, long>
            {
                {  "Skyline", 0},
                {  "Jacobi Preconditioned CG", 0}
            };

            for (int t = 0; t < repetitions; ++t)
            {
                Console.WriteLine($"Repetition: {t}");
                foreach (var solverName in solvers.Keys)
                {
                    var benchmark = new DCB(elementSize, growthLength);
                    benchmark.UniformMesh = false;
                    benchmark.UseLSM = true;
                    //TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip elements are 
                    //      found. It happens at iteration 10.
                    benchmark.MaxIterations = 10;
                    benchmark.InitializeModel();
                    //Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                    //    elementSize, benchmark.Model.Elements.Count, growthLength);

                    //Print
                    //VTKWriter writer = new VTKWriter(benchmark.model);
                    //writer.InitializeFile("dcb_transfinite");
                    //writer.CloseCurrentFile();
                    solvers[solverName].Logger.Clear();
                    IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solvers[solverName]);
                    long totalTime = solvers[solverName].Logger.CalcTotalTime();
                    Console.WriteLine($"Solver {solverName}: total time = {totalTime} ms.");
                    solverTimes[solverName] += totalTime;
                }
                Console.WriteLine();
            }

            foreach (var solverName in solverTimes.Keys)
            {
                Console.WriteLine($"Solver {solverName}: Total time = {solverTimes[solverName] / repetitions} ms");
            }

        }

        private static void ParametricCrackLength(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations)
        {
            double fineElementSize = 0.038;
            double[] growthLengths = new double[] { 0.05, 0.1, 0.2, 0.4 };
            Console.WriteLine("------------------------------------ Parametric Growth Length ------------------------------------");
            for (int i = 0; i < growthLengths.Length; ++i)
            {
                var benchmark = new DCB(fineElementSize, growthLengths[i]);
                benchmark.UniformMesh = false;
                benchmark.UseLSM = true;
                benchmark.InitializeModel();
                Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                    fineElementSize, benchmark.Model.Elements.Count, growthLengths[i]);

                try
                {
                    var solver = new SkylineSolver();
                    IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);
                    Console.WriteLine("Crack path:");
                    foreach (var point in crackPath)
                    {
                        Console.WriteLine("{0} {1}", point.X, point.Y);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    //throw e;
                }
                Console.WriteLine();
            }
        }

        private static void ParametricMesh(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations)
        {
            double propagationLength = 0.3;
            //double[] fineElementSizes = new double[] { 0.046 };
            double[] fineElementSizes = new double[] { 0.046, 0.1, 0.15, 0.2, 0.25 };
            Console.WriteLine("------------------------------------ Parametric Mesh ------------------------------------");
            for (int i = 0; i < fineElementSizes.Length; ++i)
            {
                var benchmark = new DCB(fineElementSizes[i], propagationLength);
                benchmark.UniformMesh = false;
                benchmark.UseLSM = true;
                benchmark.InitializeModel();
                Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                    fineElementSizes[i], benchmark.Model.Elements.Count, propagationLength);

                try
                {
                    var solver = new SkylineSolver();
                    IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);
                    Console.WriteLine("Crack path:");
                    foreach (var point in crackPath)
                    {
                        Console.WriteLine("{0} {1}", point.X, point.Y);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
                Console.WriteLine();
            }
        }

        private static void GridSearch(double jIntegralRadiusOverElementSize, double fractureToughness,
            int maxIterations)
        {
            double[] growthLengths = new double[] { 0.2, 0.25, 0.3 };
            double[] fineElementSizes = new double[] { 0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049 };
            //double[] growthLengths = new double[] { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
            //double[] fineElementSizes = new double[] { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 };
            //double[] fineElementSizes = new double[] {0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049 };

            for (int j = 0; j < growthLengths.Length; ++j)
            {
                for (int i = 0; i < fineElementSizes.Length; ++i)
                {
                    var benchmark = new DCB(fineElementSizes[i], growthLengths[j]);
                    benchmark.UniformMesh = false;
                    benchmark.UseLSM = true;
                    benchmark.InitializeModel();
                    Console.WriteLine("------------------ Fine mesh size = {0}, Elements = {1} , Growth length = {2} ------------------",
                        fineElementSizes[i], benchmark.Model.Elements.Count, growthLengths[j]);
                    try
                    {
                        var solver = new SkylineSolver();
                        IReadOnlyList<ICartesianPoint2D> crackPath = benchmark.Analyze(solver);
                        Console.WriteLine("Crack path:");
                        foreach (var point in crackPath)
                        {
                            Console.WriteLine("{0} {1}", point.X, point.Y);
                        }
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine(e.Message);
                    }
                    Console.WriteLine();
                }
            }
        }
    }
}
