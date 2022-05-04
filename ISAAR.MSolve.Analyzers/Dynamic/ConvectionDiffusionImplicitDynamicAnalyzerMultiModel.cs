using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Optimization: I could avoid initialization and GC of some vectors by reusing existing ones.
//TODO: Use a base class for implicit time integration methods (perhaps to together with explicit)
namespace ISAAR.MSolve.Analyzers.Dynamic
{
    public class ConvectionDiffusionImplicitDynamicAnalyzerMultiModel : INonLinearParentAnalyzer //TODO: why is this non linear
    {
        #region BDF
        private readonly int BDFOrder = 5;
        private int currentStep = 0;
        #endregion
        private readonly int maxStaggeredSteps;
        private readonly double timeStep, totalTime, tolerance;
        private IStructuralModel[] models;
        private readonly IReadOnlyDictionary<int, ILinearSystem>[] linearSystems;
        private readonly ISolver[] solvers;
        private readonly IChildAnalyzer[] childAnalyzers;
        private readonly IConvectionDiffusionIntegrationProvider[] providers;
        private readonly IVector[] initialTemperature;
        private readonly NewmarkDynamicAnalyzer structuralParentAnalyzer;
        private readonly Dictionary<int, IVector>[] rhs;
        private readonly Dictionary<int, IVector>[] stabilizingRhs;//TODO: has to be implemented, pertains to domain loads
        private readonly Dictionary<int, IVector>[] rhsPrevious;//TODO: at the moment domain loads are not implemented in this
        private readonly Dictionary<int, IVector>[] temperature;
        private readonly Dictionary<int, IVector>[,] temperatureFromPreviousStaggeredStep;
        private readonly Dictionary<int, IVector>[] capacityTimesTemperature;
        private readonly Dictionary<int, IVector>[] stabilizingConductivityTimesTemperature;
        private readonly Action<Dictionary<int, IVector>[], IStructuralModel[], ISolver[], IConvectionDiffusionIntegrationProvider[], IChildAnalyzer[]> CreateNewModel;
        #region debug
        private int i, stagSt;
        #endregion

        public ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(Action<Dictionary<int, IVector>[], IStructuralModel[], ISolver[],
            IConvectionDiffusionIntegrationProvider[], IChildAnalyzer[]> modelCreator, IStructuralModel[] models, ISolver[] solvers, IConvectionDiffusionIntegrationProvider[] providers,
            IChildAnalyzer[] childAnalyzers, double timeStep, double totalTime, int maxStaggeredSteps = 100, double tolerance = 1e-3, NewmarkDynamicAnalyzer structuralParentAnalyzer = null,
            IVector[] initialTemperature = null)
            {
            this.CreateNewModel = modelCreator;
            this.maxStaggeredSteps = maxStaggeredSteps;
            this.tolerance = tolerance;
            this.models = models;
            this.linearSystems = new IReadOnlyDictionary<int, ILinearSystem>[solvers.Length];
            rhs = new Dictionary<int, IVector>[solvers.Length];
            stabilizingRhs = new Dictionary<int, IVector>[solvers.Length];
            rhsPrevious = new Dictionary<int, IVector>[solvers.Length];
            temperature = new Dictionary<int, IVector>[solvers.Length];
            temperatureFromPreviousStaggeredStep = new Dictionary<int, IVector>[BDFOrder-1, solvers.Length];
            capacityTimesTemperature = new Dictionary<int, IVector>[solvers.Length];
            stabilizingConductivityTimesTemperature = new Dictionary<int, IVector>[solvers.Length];
            for (int i = 0; i < solvers.Length; i++)
            {
                this.linearSystems[i] = solvers[i].LinearSystems;
                rhs[i] = new Dictionary<int, IVector>();
                stabilizingRhs[i] = new Dictionary<int, IVector>();//TODO: has to be implemented, pertains to domain loads
                rhsPrevious[i] = new Dictionary<int, IVector>();//TODO: at the moment domain loads are not implemented in this
                temperature[i] = new Dictionary<int, IVector>();
                for (int st = 0; st < BDFOrder - 1; st++)
                    temperatureFromPreviousStaggeredStep[st,i] = new Dictionary<int, IVector>();
                capacityTimesTemperature[i] = new Dictionary<int, IVector>();
                stabilizingConductivityTimesTemperature[i] = new Dictionary<int, IVector>();
            }
            this.solvers = solvers;
            //solver.PreventFromOverwrittingSystemMatrices(); //TODO: If the scheme is purely implicit we can overwrite the matrix.
            this.providers = providers;
            this.ChildAnalyzer = childAnalyzers[0];
            this.childAnalyzers = childAnalyzers;
            this.structuralParentAnalyzer = structuralParentAnalyzer;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;
            this.initialTemperature = initialTemperature;
        }

        public Dictionary<int, IAnalyzerLog[]> Logs => null; //TODO: this can't be right
        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get; }
            = new Dictionary<int, ImplicitIntegrationAnalyzerLog>();

        public IChildAnalyzer ChildAnalyzer { get; }

        public void BuildMatrices()
        {
            double massFactor;
            int bdfOrderInternal = Math.Min(currentStep+1, BDFOrder);
            switch (bdfOrderInternal)
            {
                case 1: massFactor = 1; break;
                case 2: massFactor = 3.0 / 2.0; break;
                case 3: massFactor = 11.0 / 6.0; break;
                case 4: massFactor = 25.0 / 12.0; break;
                case 5: massFactor = 137.0 / 60.0; break;
                default: throw new ArgumentException("Wrong BDF Order");

            }
            var coeffs = new ImplicitIntegrationCoefficients
            {
                #region backup
                //Mass = 1 / timeStep / timeStep,
                //Damping = -1,
                //Stiffness = 1 / timeStep
                #endregion
                #region BDF
                Mass = massFactor / timeStep,
                Damping = -1,
                Stiffness = 1
                #endregion
            };
            for (int i = 0; i < linearSystems.Length; i++)
            {
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.Matrix = providers[i].LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
                }
            }
        }

        public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
        {
            #region old code
            //// u[id]: old solution
            //// v[id]: current solution
            //// vv: old acceleration
            //// v2: current acceleration
            //// v1: current velocity
            ////double vv = v2[id].Data[j];
            ////v2[id].Data[j] = a0 * (v[id].Data[j] - u[id].Data[j]) - a2 * v1[id].Data[j] - a3 * vv;
            ////v1[id].Data[j] += a6 * vv + a7 * v2[id].Data[j];

            //int id = subdomain.ID;
            //Vector<double> currentAcceleration = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> currentVelocity = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    currentAcceleration.Data[j] = a0 * (currentSolution[j] - v[id].Data[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    currentVelocity.Data[j] = v1[id].Data[j] + a6 * v2[id].Data[j] + a7 * currentAcceleration.Data[j];
            //    uu.Data[j] = a0 * currentSolution[j] + a2 * currentVelocity.Data[j] + a3 * currentAcceleration.Data[j];
            //    uc.Data[j] = a1 * currentSolution[j] + a4 * currentVelocity.Data[j] + a5 * currentAcceleration.Data[j];
            //}

            //Vector<double> tempResult = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> result = new Vector<double>(subdomain.Solution.Length);
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);

            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            //return result.Data;

            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //int id = subdomain.ID;
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    uu.Data[j] = -a0 * (v[id].Data[j] - currentSolution[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    uc.Data[j] = -a1 * (v[id].Data[j] - currentSolution[j]) - a4 * v1[id].Data[j] - a5 * v2[id].Data[j];
            //}
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);
            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            ////CalculateRhsImplicit(subdomain, result.Data, false);
            ////result.Scale(-1d);
            #endregion

            // result = M * u
            return providers[0].ConductivityMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            //if (isFirstAnalysis)
            //{
            //    // The order in which the next initializations happen is very important.
            //    model.ConnectDataStructures();
            //    solver.OrderDofs(false);
            //    foreach (ILinearSystem linearSystem in linearSystems.Values)
            //    {
            //        linearSystem.Reset(); // Necessary to define the linear system's size 
            //        linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
            //    }
            //}
            //else
            //{
            //    foreach (ILinearSystem linearSystem in linearSystems.Values)
            //    {
            //        //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
            //        //      modified. Why should the linear system be?
            //        linearSystem.Reset();
            //    }
            //}

            ////TODO: Perhaps this should be called by the child analyzer
            //BuildMatrices();

            //// Loads must be created after building the matrices.
            ////TODO: Some loads may not have to be recalculated each time the stiffness changes.
            //model.AssignLoads(solver.DistributeNodalLoads);
            //foreach (ILinearSystem linearSystem in linearSystems.Values)
            //{
            //    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            //}

            ////InitializeCoefficients();
            //InitializeInternalVectors();
            ////InitializeMatrices();
            //InitializeRhs();

            //if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            //ChildAnalyzer.Initialize(isFirstAnalysis);
        }
        private void SystemReset()
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                models[i].ConnectDataStructures();
                solvers[i].OrderDofs(false);
                //foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                //{
                //    linearSystem.Reset(); // Necessary to define the linear system's size 
                //    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                //}
            }
        }
        private void InitializeInternal()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            for (int i = 0; i < linearSystems.Length; i++)
            {
                models[i].ConnectDataStructures();
                solvers[i].OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices();

            for (int i = 0; i < linearSystems.Length; i++)
            {
                models[i].AssignLoads(solvers[i].DistributeNodalLoads);
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                }

            }

            InitializeInternalVectors();
            for (int i = 0; i < linearSystems.Length; i++)
            {
                InitializeRhs(i);
                childAnalyzers[i].Initialize(true);
            }
        }

        public void SolveTimestep(int t)
        {
            currentStep = t;
            DateTime start = DateTime.Now;
            Debug.WriteLine("Implicit Integration step: {0}", t);
            ISAAR.MSolve.Discretization.Logging.GlobalLogger.WriteLine("Implicit Integration step: " + t);

            int staggeredStep = 0;
            var temperatureNorm = 0d;
            var previousTemperatureNorm = 0d;
            var error = 1d;
            do
            {
                ISAAR.MSolve.Discretization.Logging.GlobalLogger.WriteLine("Staggered step " + staggeredStep);
                previousTemperatureNorm = temperatureNorm;
                #region backup
                //for (int i = 0; i < linearSystems.Length; i++)
                //{
                //    temperatureFromPreviousStaggeredStep[1, i].Clear();
                //    foreach (var sb in temperatureFromPreviousStaggeredStep[0, i].Keys)
                //        temperatureFromPreviousStaggeredStep[1, i].Add(sb, temperatureFromPreviousStaggeredStep[0, i][sb]);

                //    temperatureFromPreviousStaggeredStep[0, i].Clear();
                //    foreach (var linearSystem in linearSystems[i].Values)
                //    {
                //        if (linearSystem.Solution != null)
                //        {
                //            #region debug
                //            ISAAR.MSolve.Discretization.Logging.GlobalLogger.WriteLine($"linear system {i}: solution norm = {linearSystem.Solution.Norm2()}");
                //            #endregion
                //            temperatureFromPreviousStaggeredStep[0, i].Add(linearSystem.Subdomain.ID, linearSystem.Solution.Copy());
                //        }
                //        this.linearSystems[i] = solvers[i].LinearSystems;
                //    }
                //}
                #endregion
                #region theo
                for (int i = 0; i < linearSystems.Length; i++)
                {
                    foreach (var linearSystem in linearSystems[i].Values)
                    {
                        if (linearSystem.Solution != null)
                        {
                            #region debug
                            ISAAR.MSolve.Discretization.Logging.GlobalLogger.WriteLine($"linear system {i}: solution norm = {linearSystem.Solution.Norm2()}");
                            #endregion
                        }
                        this.linearSystems[i] = solvers[i].LinearSystems;
                    }
                }
                #endregion

                InitializeInternal();
                for (int i = 0; i < linearSystems.Length; i++)
                {
                    IDictionary<int, IVector> rhsVectors = providers[i].GetRhsFromHistoryLoad(t);
                    foreach (var l in linearSystems[i].Values) l.RhsVector = rhsVectors[l.Subdomain.ID];
                    InitializeRhs(i);
                    CalculateRhsImplicit(i);
                    childAnalyzers[i].Solve();
                }

                CreateNewModel(temperature, models, solvers, providers, childAnalyzers);

                temperatureNorm = 0;
                if (structuralParentAnalyzer != null)
                {
                    temperatureNorm = 0;
                    structuralParentAnalyzer.SolveTimestep(t, staggeredStep);
                    foreach (var linearSystem in structuralParentAnalyzer.linearSystems.Values)
                    {
                        temperatureNorm += linearSystem.Solution.Norm2();
                    }
                }

                for (int i = 0; i < linearSystems.Length; i++)
                {
                    foreach (var linearSystem in linearSystems[i].Values)
                    {
                        temperatureNorm += linearSystem.Solution.Norm2();
                    }
                }
                error = temperatureNorm != 0 ? Math.Abs(temperatureNorm - previousTemperatureNorm) / temperatureNorm : 0;
                Debug.WriteLine("Staggered step: {0} - error {1}", staggeredStep, error);
                staggeredStep++;
            }
            while (staggeredStep < maxStaggeredSteps && error > tolerance);

            DateTime end = DateTime.Now;
            SystemReset();
            UpdateTemperature(t);
            UpdateResultStorages(start, end);
            Debug.WriteLine("-------------");
        }

        public void Solve()
        {
            int numTimeSteps = (int)(totalTime / timeStep);
            for (int t = 0; t < numTimeSteps; ++t)
            {
                SolveTimestep(t);
            }
        }

        private void CalculateRhsImplicit(int modelNo)
        {
            foreach (ILinearSystem linearSystem in linearSystems[modelNo].Values)
            {
                linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, modelNo, true);
            }
        }

        private IVector CalculateRhsImplicit(ILinearSystem linearSystem, int modelNo, bool addRhs)
        {
            #region backup
            ////TODO: what is the meaning of addRhs? Do we need this when solving dynamic thermal equations?
            ////result = a0*conductuvity*temperature + a2* rhs - StabilizingRhs
            //double a0 = 1 / Math.Pow(timeStep, 2);
            ////double a1 = 1 / (2 * timeStep);
            //double a2 = 1 / timeStep;
            //int id = linearSystem.Subdomain.ID;
            //capacityTimesTemperature[modelNo][id] = providers[modelNo].CapacityMatrixVectorProduct(linearSystem.Subdomain, temperature[modelNo][id]);
            //capacityTimesTemperature[modelNo][id].ScaleIntoThis(a0);
            //rhs[modelNo][id].ScaleIntoThis(a2);
            //var rhsResult = capacityTimesTemperature[modelNo][id].Subtract(stabilizingRhs[modelNo][id]);
            //rhsResult.AddIntoThis(rhs[modelNo][id]);

            //return rhsResult;
            #endregion

            #region BDF
            //rhs_hat = rhs + (temperatureTerm) * MM (CapacityMatrix) * a2
            double a2 = 1 / timeStep;
            int id = linearSystem.Subdomain.ID;
            int bdfOrderInternal = Math.Min(currentStep+1, BDFOrder);

            double[] rhsFactors = new double[bdfOrderInternal];
            switch (bdfOrderInternal)
            {
                case 1:
                    rhsFactors[0] = 1; // T_n
                    break;
                case 2:
                    rhsFactors[0] = 4.0 / 2.0; // T_n == n+1
                    rhsFactors[1] = -1.0 / 2.0; // T_n-1 == n
                    break;
                case 3:
                    rhsFactors[0] = 18.0 / 6.0; // T_n == n+2
                    rhsFactors[1] = -9.0 / 6.0; // T_n-1 == n+1
                    rhsFactors[2] = 2.0 / 6.0; // T_n-2 == n
                    break;
                case 4:
                    rhsFactors[0] =  48.0 / 12.0; // T_n == n+3
                    rhsFactors[1] = -36.0 / 12.0; // T_n-1 == n+2
                    rhsFactors[2] =  16.0 / 12.0; // T_n-2 == n+1
                    rhsFactors[3] = -3.0 / 12.0; // T_n-3 == n
                    break;
                case 5:
                    rhsFactors[0] = 300.0 / 60.0; // T_n == n+4
                    rhsFactors[1] = -300.0 / 60.0; // T_n-1 == n+3
                    rhsFactors[2] = 200.0 / 60.0; // T_n-2 == n+2
                    rhsFactors[3] = -75.0 / 60.0; // T_n-4 == n+1
                    rhsFactors[4] = 12.0 / 60.0; // T_n-5 == n
                    break;
                default: throw new ArgumentException("Wrong BDF Order");
            }

            var temperatureTerm = temperature[modelNo][id].Scale(rhsFactors[0]);
            for (int bdfTerm = 1; bdfTerm < bdfOrderInternal; bdfTerm++)
            {
                if (temperatureFromPreviousStaggeredStep[bdfTerm - 1, modelNo].ContainsKey(id))
                    temperatureTerm.AddIntoThis(temperatureFromPreviousStaggeredStep[bdfTerm - 1, modelNo][id].Scale(rhsFactors[bdfTerm]));
                else
                {//this time step has not yet occured, search for the last initialized
                    for (int jj = bdfTerm-2; jj>=-1; jj--)
                    {
                        if (jj >= 0)
                        {
                            if (temperatureFromPreviousStaggeredStep[jj, modelNo].ContainsKey(id))
                            {
                                temperatureTerm.AddIntoThis(temperatureFromPreviousStaggeredStep[jj, modelNo][id].Scale(rhsFactors[bdfTerm]));
                                break;
                            }
                        }
                        else
                        {
                            temperatureTerm.AddIntoThis(temperature[modelNo][id].Scale(rhsFactors[bdfTerm]));
                            break;
                        }
                    }
                }
            }

            // MM * (temperatureFactor)
            capacityTimesTemperature[modelNo][id] = providers[modelNo].CapacityMatrixVectorProduct(linearSystem.Subdomain, temperatureTerm);
            // MM * (temperatureFactor) * a2
            capacityTimesTemperature[modelNo][id].ScaleIntoThis(a2);
            // MM * (temperatureFactor) * a2 - StabilizingRhs
            var rhsResult = capacityTimesTemperature[modelNo][id].Subtract(stabilizingRhs[modelNo][id]);
            // F +  MM * (temperatureFactor) * a2 - StabilizingRhs
            rhsResult.AddIntoThis(rhs[modelNo][id]);
            return rhsResult;
            #endregion
        }

        private void InitializeInternalVectors()
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                //temperature[i].Clear();
                stabilizingRhs[i].Clear();
                rhs[i].Clear();
                rhsPrevious[i].Clear();
                stabilizingConductivityTimesTemperature[i].Clear();
                capacityTimesTemperature[i].Clear();

                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    capacityTimesTemperature[i].Add(id, linearSystem.CreateZeroVector());
                    stabilizingConductivityTimesTemperature[i].Add(id, linearSystem.CreateZeroVector());
                    //temperature.Add(id, linearSystem.CreateZeroVector());
                    stabilizingRhs[i].Add(id, linearSystem.CreateZeroVector());
                    rhs[i].Add(id, linearSystem.CreateZeroVector());
                    rhsPrevious[i].Add(id, linearSystem.CreateZeroVector());

                    if (temperature[i].ContainsKey(id) == false)
                        temperature[i].Add(id, initialTemperature[i].Copy());

                    //// Account for initial conditions coming from a previous solution. 
                    ////TODO: This doesn't work as intended. The solver (previously the LinearSystem) initializes the solution to zero.
                    //if (linearSystem.Solution != null) temperature[i].Add(id, linearSystem.Solution.Copy());
                    //else temperature[i].Add(id, initialTemperature[i].Copy());
                }
            }
        }

        private void InitializeRhs(int modelNo)
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 0,  //never used
                Stiffness = 0
            };
            foreach (ILinearSystem linearSystem in linearSystems[modelNo].Values)
            {
                providers[modelNo].ProcessRhs(linearSystem.Subdomain, linearSystem.RhsVector);
                rhs[modelNo][linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy(); //TODO: copying the vectors is wasteful.
                stabilizingRhs[modelNo][linearSystem.Subdomain.ID] = providers[modelNo].StabilizingRhs(linearSystem.Subdomain);
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    if (ResultStorages.ContainsKey(id))
                        if (ResultStorages[id] != null)
                            foreach (var l in ChildAnalyzer.Logs[id])
                                ResultStorages[id].StoreResults(start, end, l);
                }
            }
        }

        private void UpdateTemperature(int timeStep)
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                #region theo
                //for (int jj = 1; jj < Math.Min(currentStep + 1, BDFOrder - 1); jj++)
                for (int jj = Math.Min(currentStep, BDFOrder - 2); jj >= 1 ; jj--)
                {
                    temperatureFromPreviousStaggeredStep[jj, i].Clear();
                    foreach (var sb in temperatureFromPreviousStaggeredStep[jj-1, i].Keys)
                        temperatureFromPreviousStaggeredStep[jj, i].Add(sb, temperatureFromPreviousStaggeredStep[jj-1, i][sb].Copy());
                }

                temperatureFromPreviousStaggeredStep[0, i].Clear();
                #endregion
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    #region theo
                    temperatureFromPreviousStaggeredStep[0, i].Add(id, temperature[i][id].Copy());
                    temperature[i][id].CopyFrom(linearSystem.Solution);
                    #endregion
                    //temperature[i][id].AddIntoThis(linearSystem.Solution);
                    if ((timeStep + 1) % 1 == 0)
                    {
                        var path0 = Path.Combine(Directory.GetCurrentDirectory(), $"MSolveOutput");
                        var path = Path.Combine(path0, $"temperature{i}-{timeStep}.txt");
                        var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainVertical };
                        writer.WriteToFile(temperature[i][id], path);
                        //writer.WriteToFile(temperature[id][0], path1);

                        //File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                        //File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                    }
                }
            }
            





        }
    }
}
