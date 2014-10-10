using CNNWB.CNN;
using CNNWB.Common;
using CNNWB.Data;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace CNNWB.ConsoleApp.Test
{
    class Program
    {
        static DataProvider _dp;
        static TrainingRate _data { get; set; }

        static void Main(string[] args)
        {
            // Load MINST train data

            _storageDirectory = @"D:\tmp\cnnwb_storage";

            _dp = new DataProvider(_storageDirectory);
            _dp.LoadDataset(DataProviderSets.MNIST);

            // Create Network LeNet5
            NeuralNetwork cnn = InitializeDefaultNeuralNetwork(_dp);

            // Train

            cnn.RaiseNetworkProgressEvent += cnn_RaiseNetworkProgressEvent;

            CancellationTokenSource cts = new CancellationTokenSource();
            CancellationToken token = cts.Token;

            cnn.StartTraining(token).Wait();

            Console.WriteLine("-----------------------------------------");
            cnn.SaveWeights(Path.Combine(_storageDirectory, "mnist-weights.bin"));

            Console.ReadLine();

        }

        static void cnn_RaiseNetworkProgressEvent(object sender, EventArgs e)
        {
            var NeuralNetwork = sender as NeuralNetwork;

            //Console.WriteLine();
            //Console.WriteLine("Task Duration: {0}", NeuralNetwork.GetTaskDuration());
            //Console.WriteLine("Sample Rate: {0}", NeuralNetwork.SampleRate.ToString("N2"));
            Console.Write(".");

            switch (NeuralNetwork.OperationState)
            {
                case NetworkStates.CalculatingHessian:
                    //Console.WriteLine();
                    //Console.WriteLine("Calculating Pseudo-Hessian...");
                    //Console.WriteLine("Progress Value: {0}", (double)(NeuralNetwork.SampleIndex + 1) / NeuralNetwork.HessianSamples);
                    //Console.WriteLine("Progress Maximum, HessianSamples: {0}", NeuralNetwork.HessianSamples);
                    //Console.WriteLine("Progress Value 2, SampleIndex: {0}", NeuralNetwork.SampleIndex);

                    break;

                case NetworkStates.Training:
                    //Console.WriteLine();
                    //Console.WriteLine("Training...");

                    //Console.WriteLine("Task Progress Value, SampleIndex: {0}",
                    //        (double)(NeuralNetwork.SampleIndex + 1) / _dp.TrainingSamplesCount);
                    //Console.WriteLine("ProgressBar.Maximum, TrainingSamplesCount: {0}", _dp.TrainingSamplesCount);
                    //Console.WriteLine("ProgressBar.Value, TrainingSamplesCount: {0}", NeuralNetwork.SampleIndex);

                    break;

                case NetworkStates.Testing:
                    Console.WriteLine();
                    Console.WriteLine("Testing...");
                    //if (testingPVM.UseTrainingSet)
                    //{
                    //    taskBarItemInfo.ProgressValue = (double)(pageViewModel.NeuralNetwork.SampleIndex + 1) / pageViewModel.DataProvider.TrainingSamplesCount;
                    //    MainView.ProgressBar.Maximum = pageViewModel.DataProvider.TrainingSamplesCount;
                    //    if (pageViewModel.NeuralNetwork.SampleIndex < pageViewModel.DataProvider.TrainingSamplesCount)
                    //        testingPVM.SampleImage = pageViewModel.DataProvider.TrainingSamples[pageViewModel.NeuralNetwork.SampleIndex].GetBitmapSource(pageViewModel.DataProvider);
                    //}
                    //else
                    //{
                    //    taskBarItemInfo.ProgressValue = (double)(pageViewModel.NeuralNetwork.SampleIndex + 1) / pageViewModel.DataProvider.TestingSamplesCount;
                    //    MainView.ProgressBar.Maximum = pageViewModel.DataProvider.TestingSamplesCount;
                    //    if (pageViewModel.NeuralNetwork.SampleIndex < pageViewModel.DataProvider.TestingSamplesCount)
                    //        testingPVM.SampleImage = pageViewModel.DataProvider.TestingSamples[pageViewModel.NeuralNetwork.SampleIndex].GetBitmapSource(pageViewModel.DataProvider);
                    //}
                    //MainView.ProgressBar.Value = pageViewModel.NeuralNetwork.SampleIndex;
                    //progressText.Length = 0;
                    //{
                    //    double testAccuracy = (double)((pageViewModel.NeuralNetwork.SampleIndex + 1) - pageViewModel.NeuralNetwork.TestErrors) * (100D / (pageViewModel.NeuralNetwork.SampleIndex + 1));
                    //    progressText.AppendFormat(CultureInfo.CurrentUICulture, "\nSample Index:\t{0:G}\nAverage Loss:\t{1:N10}\nTest Error %:\t{2:N4}\nTest Errors:\t{3:G}\nAccuracy %:\t{4:N4}", pageViewModel.NeuralNetwork.SampleIndex, pageViewModel.NeuralNetwork.AvgTestLoss, 100D - testAccuracy, pageViewModel.NeuralNetwork.TestErrors, testAccuracy);
                    //    testingPVM.ProgressText = progressText.ToString();
                    //    if (testingPVM.DataProvider.CurrentDataSet == DataProviderSets.CIFAR10)
                    //        testingPVM.SampleLabel = ((CIFAR10Classes)pageViewModel.NeuralNetwork.CurrentSample.Label).ToString();
                    //    else
                    //        testingPVM.SampleLabel = pageViewModel.NeuralNetwork.CurrentSample.Label.ToString();
                    //}
                    break;

                case NetworkStates.CalculatingTestError:
                    //Console.WriteLine();
                    //Console.WriteLine("CalculatingTestError...");
                    //taskBarItemInfo.ProgressValue = (double)(pageViewModel.NeuralNetwork.SampleIndex + 1) / pageViewModel.DataProvider.TestingSamplesCount;
                    //MainView.ProgressBar.Maximum = pageViewModel.DataProvider.TestingSamplesCount;
                    //MainView.ProgressBar.Value = pageViewModel.NeuralNetwork.SampleIndex;
                    //trainingPVM.SampleImage = null;
                    //progressText.Length = 0;
                    //progressText.Append("Calculating Test Error...");
                    //{
                    //    double testAccuracy = ((pageViewModel.NeuralNetwork.SampleIndex + 1) - pageViewModel.NeuralNetwork.TestErrors) * (100D / (pageViewModel.NeuralNetwork.SampleIndex + 1));
                    //    progressText.AppendFormat(CultureInfo.CurrentUICulture, "\n\nSample Index:\t{0:G}\nAverage Loss:\t{1:N10}\nTest Error %:\t{2:N4}\nTest Errors:\t{3:G}\nAccuracy %:\t{4:N4}", pageViewModel.NeuralNetwork.SampleIndex, pageViewModel.NeuralNetwork.AvgTestLoss, 100D - testAccuracy, pageViewModel.NeuralNetwork.TestErrors, testAccuracy);
                    //}
                    //trainingPVM.ProgressText = progressText.ToString();
                    //trainingPVM.SampleLabel = String.Empty;
                    break;

                case NetworkStates.SavingWeights:
                    Console.WriteLine();
                    Console.WriteLine("SavingWeights...");
                    //MainView.SampleRate.Text = String.Empty;
                    //taskBarItemInfo.ProgressState = System.Windows.Shell.TaskbarItemProgressState.None;
                    //MainView.ProgressBar.Value = 0;
                    //progressText.Length = 0;
                    //progressText.Append("Saving Weights...");
                    //trainingPVM.ProgressText = progressText.ToString();
                    //trainingPVM.SampleLabel = String.Empty;
                    break;

                case NetworkStates.NewEpoch:
                    Console.WriteLine();
                    Console.WriteLine("NewEpoch...");

                    TimeSpan span = NeuralNetwork.TaskDuration.Elapsed.Subtract(NeuralNetwork.EpochDuration);

                    NeuralNetwork.EpochDuration = NeuralNetwork.TaskDuration.Elapsed;
                    double testAccuracy = (_dp.TestingSamplesCount - NeuralNetwork.TestErrors) * (100D / _dp.TestingSamplesCount);

                    Console.WriteLine("TrainingResult: Epoch: {0}, DistortionPercentage: {1}, Rate: {2}, AvgTrainLoss: {3}, TrainErrors: {4}, TrainErrorsPercent: {5}, AvgTestLoss: {6}, TestErrors: {7}, 100 - trainAccuracy: {8}, trainAccuracy: {9}, Time: {10}",
                                        NeuralNetwork.CurrentEpoch,
                                        NeuralNetwork.TrainingRate.Distorted ? NeuralNetwork.TrainingRate.DistortionPercentage : 0,
                                        NeuralNetwork.TrainingRate.Rate,
                                        NeuralNetwork.AvgTrainLoss,
                                        NeuralNetwork.TrainErrors,
                                        (double)NeuralNetwork.TrainErrors * (100D / _dp.TrainingSamplesCount),
                                        NeuralNetwork.AvgTestLoss,
                                        NeuralNetwork.TestErrors,
                                        100D - testAccuracy,
                                        testAccuracy,
                                        new TimeSpan(span.Hours, span.Minutes, span.Seconds));

                    NeuralNetwork.OperationState = NetworkStates.Idle;

                    break;
            }
        }


        private static NeuralNetwork InitializeDefaultNeuralNetwork(DataProvider dp)
        {
            using (CNNDataSet.TrainingRatesDataTable table = new CNNDataSet.TrainingRatesDataTable())
            {
                table.ReadXml(@"D:\prj\cnnwb\CNNWB.Data\TrainingSchemes\LeCun2.scheme-xml");
                CNNDataSet.TrainingRatesRow row = table.Rows[0] as CNNDataSet.TrainingRatesRow;
                _data = new TrainingRate(row.Rate, row.Epochs, row.MinimumRate, row.WeightDecayFactor, row.Momentum, row.BatchSize, row.InitialAvgLoss, row.DecayFactor, row.DecayAfterEpochs, row.WeightSaveTreshold, row.Distorted, row.DistortionPercentage, row.SeverityFactor, row.MaxScaling, row.MaxRotation, row.ElasticSigma, row.ElasticScaling);
            }

            //NeuralNetwork network = new NeuralNetwork(_dp, "LeNet-5", 10, 0.8D, LossFunctions.MeanSquareError,
            //                            DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.MaxDegreeOfParallelism = 2;
            
            //network.AddGlobalTrainingRate(_data, true);

            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 6, 28, 28, 5, 5);
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 6, 14, 14, 2, 2);

            //bool[] maps = new bool[6 * 16] 
            //{
            // true, false,false,false,true, true, true, false,false,true, true, true, true, false,true, true,
            // true, true, false,false,false,true, true, true, false,false,true, true, true, true, false,true,
            // true, true, true, false,false,false,true, true, true, false,false,true, false,true, true, true,
            // false,true, true, true, false,false,true, true, true, true, false,false,true, false,true, true,
            // false,false,true, true, true, false,false,true, true, true, true, false,true, true, false,true,
            // false,false,false,true, true, true, false,false,true, true, true, true, false,true, true, true
            //};

            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 16, 10, 10, 5, 5, mappings: new Mappings(maps));
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 16, 5, 5, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 120, 1, 1, 5, 5);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);

            //network.InitializeWeights();

            NeuralNetwork network = new NeuralNetwork(_dp, "Simard-6", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            network.AddLayer(LayerTypes.Input, 1, 32, 32);
            network.AddLayer(LayerTypes.ConvolutionalSubsampling, ActivationFunctions.Tanh, 6, 14, 14, 5, 5);
            network.AddLayer(LayerTypes.ConvolutionalSubsampling, ActivationFunctions.Tanh, 50, 5, 5, 5, 5);
            network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 100);
            network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);

            network.MaxDegreeOfParallelism = 3;

            network.AddGlobalTrainingRate(_data, true);

            network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "Simard-16", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.1D);
            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.ConvolutionalSubsampling, ActivationFunctions.Tanh, 16, 14, 14, 5, 5);
            //network.AddLayer(LayerTypes.ConvolutionalSubsampling, ActivationFunctions.Tanh, 64, 5, 5, 5, 5);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 196);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "LeNet-5", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 6, 28, 28, 5, 5, 1, 1);
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 6, 14, 14, 2, 2, 2, 2);
            //bool[] maps = new bool[6 * 16]
            //{
            //    true,  false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,
            //    true,  true,  false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,
            //    true,  true,  true,  false, false, false, true,  true,  true,  false, false, true,  false, true,  true,  true,
            //    false, true,  true,  true,  false, false, true,  true,  true,  true,  false, false, true,  false, true,  true,
            //    false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,  false, true,
            //    false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 16, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 16, 5, 5, 2, 2, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 120, 1, 1, 5, 5, 1, 1);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "LeNet-5", 10, 1D, LossFunctions.CrossEntropy, DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 6, 28, 28, 5, 5);
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Ident, 6, 14, 14, 2, 2, 2, 2);
            //bool[] maps = new bool[6 * 16]
            //{
            //    true,  false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,
            //    true,  true,  false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,
            //    true,  true,  true,  false, false, false, true,  true,  true,  false, false, true,  false, true,  true,  true,
            //    false, true,  true,  true,  false, false, true,  true,  true,  true,  false, false, true,  false, true,  true,
            //    false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,  false, true,
            //    false, false, false, true,  true,  true,  false, false, true,  true,  true,  true,  false, true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 16, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Ident, 16, 5, 5, 2, 2, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 120, 1, 1, 5, 5);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftMax, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "MyNet-16", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 16, 28, 28, 5, 5);
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 16, 14, 14, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 64, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(16, 64, 66, 1));
            //network.AddLayer(LayerTypes.AvgPooling, ActivationFunctions.Tanh, 64, 5, 5, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 256, 1, 1, 5, 5, 1, 1, 0, 0, new Mappings(64, 256, 66, 2));
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-A", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 28, 28, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 14, 14, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(64, 64, 66, 1));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 5, 5, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 64, 3, 3, 3, 3, 1, 1, 0, 0, new Mappings(64, 64, 66, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 386, 1, 1, 3, 3, 1, 1, 0, 0, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftSign, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-Z2", 10, 1D, LossFunctions.CrossEntropy, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};

            //bool[] maps = new bool[3 * 96]
            //{
            //    true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false,
            //    false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false,
            //    false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true,  false, false, true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 28, 28, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 14, 14, 3, 3, 2, 2);
            ////network.AddLayer(LayerTypes.LocalResponseNormalizationCM, ActivationFunctions.None, 64, 14, 14, 3, 3);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(64, 64, 66, 1));
            ////network.AddLayer(LayerTypes.LocalResponseNormalizationCM, ActivationFunctions.None, 64, 10, 10, 3, 3);
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 5, 5, 3, 3, 2, 2);
            ////network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 64, 1, 1, 5, 5, 1, 1, 0, 0, new Mappings(64, 64, 66, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.Logistic, 384, 1, 1, 5, 5, 1, 1, 0, 0, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftMax, 10);
            //network.InitializeWeights();
            //network.LoadWeights(StorageDirectory + @"\CNN-CIFAR10-Z2 (2259 errors).weights-bin");

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-B2", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 28, 28, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 14, 14, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(64, 64, 66, 1));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 5, 5, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.Tanh, 512, 1, 1, 5, 5, 1, 1, 0, 0, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-Z3", 10, 1D, LossFunctions.CrossEntropy, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 32, 32, 5, 5, 1, 1, 2, 2, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 16, 16, 3, 3, 2, 2);
            ////network.AddLayer(LayerTypes.LocalResponseNormalizationCM, ActivationFunctions.None, 64, 16, 16, 5, 5);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 16, 16, 5, 5, 1, 1, 2, 2, new Mappings(64, 64, 66, 1));
            ////network.AddLayer(LayerTypes.LocalResponseNormalizationCM, ActivationFunctions.None, 64, 16, 16, 5, 5);
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 8, 8, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 32, 8, 8, 3, 3, 1, 1, 1, 1, new Mappings(64, 32, 50, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.SoftSign, 32, 8, 8, 3, 3, 1, 1, 1, 1, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftMax, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-C", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 28, 28, 5, 5, 1, 1, 0, 0, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 14, 14, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 10, 10, 5, 5, 1, 1, 0, 0, new Mappings(64, 64, 66, 1));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 5, 5, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 32, 5, 5, 3, 3, 1, 1, 1, 1, new Mappings(64, 32, 50, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 32, 5, 5, 3, 3, 1, 1, 1, 1, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-Z4", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D, 2000, false);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            ////bool[] maps = new bool[3 * 64]
            ////{
            ////    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, false, true,  false, true,  false, true,  false, true,
            ////    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, true, true, true, true, true, true,  false, true,  false, true,  false, true,  false,
            ////    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, true, true, true, true, true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, true, true, true, true, true, false, true,  false, true,  false, true,  false, true
            ////};
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            ////bool[] maps = new bool[3 * 48]    //Z3
            ////{
            ////    true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
            ////    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
            ////    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true 
            ////};

            ////bool[] maps = new bool[3 * 48]    //Z3
            ////{
            ////    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false,
            ////    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            ////    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, true, true, true, true, true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true
            ////};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 32, 32, 5, 5, 1, 1, 2, 2, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 16, 16, 3, 3, 2, 2);
            ////network.AddLayer(LayerTypes.LocalContrastNormalization, ActivationFunctions.None, 64, 16, 16, 5, 5);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 16, 16, 5, 5, 1, 1, 2, 2, new Mappings(64, 64, 66, 1));
            ////network.AddLayer(LayerTypes.LocalContrastNormalization, ActivationFunctions.None, 64, 16, 16, 5, 5);
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 8, 8, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 32, 8, 8, 3, 3, 1, 1, 1, 1, new Mappings(64, 32, 50, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 32, 8, 8, 3, 3, 1, 1, 1, 1, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftSign, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-D", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 64]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true, 
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 32, 32, 5, 5, 1, 1, 2, 2, new Mappings(maps));
            //network.AddLayer(LayerTypes.MaxPoolingWeightless, ActivationFunctions.Ident, 64, 16, 16, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 16, 16, 3, 3, 1, 1, 1, 1, new Mappings(64, 64, 50, 1));
            //network.AddLayer(LayerTypes.AvgPoolingWeightless, ActivationFunctions.Ident, 64, 8, 8, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.BReLU, 32, 8, 8, 3, 3, 1, 1, 1, 1, new Mappings(64, 32, 50, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.Tanh, 32, 8, 8, 3, 3, 1, 1, 1, 1, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            ////network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-CIFAR10-G", 10, 0.8D, LossFunctions.MeanSquareError, DataProviderSets.CIFAR10, TrainingStrategy.SGD, 0.02D);
            //network.AddLayer(LayerTypes.Input, 3, 32, 32);
            //bool[] maps = new bool[3 * 48]
            //{
            //    true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false,
            //    false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false, false, false, false, false, false, true,  true,  true, true, true, true, true, true, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            //    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, true, true, true, true, true, false, false, false, false, false, false, false, false, true,  true,  true,  true,  true,  true,  true,  true
            //};
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 48, 32, 32, 5, 5, 1, 1, 2, 2, new Mappings(maps));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 48, 16, 16, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 64, 16, 16, 5, 5, 1, 1, 2, 2, new Mappings(48, 64, 66, 1));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 64, 8, 8, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 24, 8, 8, 3, 3, 1, 1, 1, 1, new Mappings(64, 24, 50, 2));
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.Tanh, 24, 8, 8, 3, 3, 1, 1, 1, 1, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
            //network.InitializeWeights();

            //NeuralNetwork network = new NeuralNetwork(DataProvider, "CNN-MNIST-A", 10, 1D, LossFunctions.CrossEntropy, DataProviderSets.MNIST, TrainingStrategy.SGD, 0.02D);
            //network.AddLayer(LayerTypes.Input, 1, 32, 32);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 32, 24, 24, 9, 9, 1, 1, 0, 0);
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 32, 12, 12, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Convolutional, ActivationFunctions.ReLU, 32, 8, 8, 5, 5, 1, 1, 0, 0, new Mappings(32, 32, 66));
            //network.AddLayer(LayerTypes.StochasticPooling, ActivationFunctions.Ident, 32, 4, 4, 3, 3, 2, 2);
            //network.AddLayer(LayerTypes.Local, ActivationFunctions.ReLU, 256, 1, 1, 4, 4, 1, 1, 0, 0, 50);
            //network.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.SoftMax, 10);
            //network.InitializeWeights();

            //network.RaiseNetworkProgressEvent += new EventHandler<EventArgs>(NetworkProgressEvent);
            //network.RaiseAddUnrecognizedTestSampleEvent += new EventHandler<AddUnrecognizedTestSampleEventArgs>(AddUnrecognizedTestSampleEvent);

            return network;
        }

        public static string _storageDirectory { get; set; }
    }
}
