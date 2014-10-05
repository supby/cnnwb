using CNNWB.Common;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Reflection;
using System.Resources;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Timers;
using System.Windows.Media.Imaging;
using TaskDialogInterop;

namespace CNNWB.Model
{
	// warning: all enums must start from zero and always be incremental by excactly one otherwise the CNNDataSet constraints are violated!
	public enum TaskState
	{
		Paused = 0,
		Running = 1,
		Stopped = 2
	}

	public enum TrainingStrategy
	{
		MiniBatchSGD = 0,                       // Mini-Batch Stochastic Gradient Descent with weight decay (can be used for both objectives: mse and ce)
        MiniBatchSGDLevenbergMarquardt = 1,     // Mini-Batch Stochastic Gradient Descent with Levenberg-Marquardt second order method (mostly only usable with Mean Square Error objective function)
        MiniBatchSGDLevenbergMarquardtModA = 2, // Mini-Batch Stochastic Gradient Descent with Levenberg-Marquardt Modification A (this allows to use the LM method with the cross-entropy objective function)
        MiniBatchSGDM = 3,                      // Mini-Batch Stochastic Gradient Descent with Momentum
		SGD = 4,                                // Stochastic Gradient Descent with weight decay (or without by setting the weight decay parameter to zero)
        SGDLevenbergMarquardt = 5,              // Stochastic Gradient Descent with Levenberg-Marquardt second order method (mostly only usable with Mean Square Error objective function)
        SGDLevenbergMarquardtModA = 6,          // Stochastic Gradient Descent with Levenberg-Marquardt second order method Modification A (for cross-entropy objective)
        SGDM = 7                                // Stochastic Gradient Descent with Momentum
	}

	public enum LossFunctions
	{
		CrossEntropy = 0,
		MeanSquareError = 1
	}

	public enum LayerTypes
	{
		AvgPooling = 0,
		AvgPoolingWeightless = 1,
		Convolutional = 2,                  
		ConvolutionalSubsampling = 3,           // Patrice Simard's layertype
		FullyConnected = 4,
		Input = 5,
		L2Pooling = 6,
		Local = 7,
		LocalContrastNormalization = 8,         // same map
		LocalResponseNormalization = 9,
		LocalResponseNormalizationCM = 10,      // across maps
		MaxPooling = 11,
		MaxPoolingWeightless = 12,
		RBF = 13,                               // as a possible final output layer for the MNIST dataset
		StochasticPooling = 14
	}

	public enum ActivationFunctions
	{
		Abs = 0,
		AbsTanh = 1,
		BReLU1 = 2,                             // Bounded Rectifier Linear ]0-1]
        BReLU6 = 3,                             //    "        "        "   ]0-6]
		Gaussian = 4,
		Ident = 5,
		Logistic = 6,                           // aka Sigmoid function
		None = 7,
		Ramp = 8,
        ReLU = 9,                               // Rectified linear function ]0-∞]
		STanh = 10,                             // Symmetric Tanh
		SoftMax = 11,
		SoftReLU = 12,                          // aka SoftPlus
		SoftSign = 13,
		Square = 14,
		SquareRoot = 15,
		Tanh = 16,
	}

	public enum NetworkStates
	{
		Idle = 0,
		NewEpoch = 1,
		CalculatingHessian = 2,
		Testing = 3,
		Training = 4,
		CalculatingTestError = 5,
		SavingWeights = 6
	}

	public sealed class Mappings
	{
		public bool[] Mapping;

		public Mappings(bool[] mapping)
		{
			if ((mapping == null))
				throw new ArgumentException("Invalid Mappings parameter(s)");

			Mapping = mapping;
		}

		public Mappings(int previousLayerMapCount, int currentLayerMapCount, int density, int randomSeed = 0)
		{
			if ((previousLayerMapCount < 1) || (currentLayerMapCount < 1))
				throw new ArgumentException("Invalid Mappings parameter(s)");

			Mapping = new bool[previousLayerMapCount * currentLayerMapCount];
			Random random = new Random(randomSeed);
					   
			for (int channel = 0; channel < previousLayerMapCount; channel++)
				for (int map = 0; map < currentLayerMapCount; map++)
					Mapping[(channel * currentLayerMapCount) + map] = (random.Next(100) < density);
		}

		public bool IsMapped(int map, int previousLayerMapCount, int currentLayerMapCount)
		{
			return (Mapping[map + (previousLayerMapCount * currentLayerMapCount)]);
		}
	}

	[StructLayout(LayoutKind.Sequential, Pack = 4)]
	public struct Connection
	{
		public int ToNeuronIndex;
		public int ToWeightIndex;

		public Connection(int toNeuronIndex, int toWeightIndex)
		{
			ToNeuronIndex = toNeuronIndex;
			ToWeightIndex = toWeightIndex;
		}
	}

	[StructLayout(LayoutKind.Sequential, Pack = 8)]
	public struct Weight
	{
		public double Value;
		public double DiagonalHessian;
		public double D1Err;
	}

	[StructLayout(LayoutKind.Sequential, Pack = 8)]
	public struct Neuron
	{
		public double Output;
		public double D1ErrX;
	}

	public sealed class TestingParameters
	{
		public bool UseTrainingSamples { get; set; }  
		public bool Distorted { get; set; }
		public int DistortionPercentage { get; set; }
		public double SeverityFactor { get; set; }
		public double MaxScaling { get; set; }
		public double MaxRotation { get; set; }
		public double ElasticSigma { get; set; }
		public double ElasticScaling { get; set; }

		public TestingParameters()
		{
			UseTrainingSamples = false;
			Distorted = Properties.Settings.Default.Distorted;
			DistortionPercentage = Properties.Settings.Default.DistortionPercentage;
			SeverityFactor = Properties.Settings.Default.SeverityFactor;
			MaxScaling = Properties.Settings.Default.MaxScaling;
			MaxRotation = Properties.Settings.Default.MaxRotation;
			ElasticSigma = Properties.Settings.Default.ElasticSigma;
			ElasticScaling = Properties.Settings.Default.ElasticScaling;
		}

		public TestingParameters(bool useTrainingSamples, bool distorted, int distortionPercentage, double severityFactor, double maxScaling, double maxRotation, double elasticSigma, double elasticScaling)
		{
			UseTrainingSamples = useTrainingSamples;
			Distorted = distorted;
			DistortionPercentage = distortionPercentage;
			SeverityFactor = severityFactor;
			MaxScaling = maxScaling;
			MaxRotation = maxRotation;
			ElasticSigma = elasticSigma;
			ElasticScaling = elasticScaling;
		}
	}

	public sealed class AddUnrecognizedTestSampleEventArgs : EventArgs
	{
		public int SampleIndex { get; private set; }
		public int WrongValue { get; private set; }
		public ImageData Sample {get; private set;}

		public AddUnrecognizedTestSampleEventArgs(int sampleIndex, int wrongValue, ImageData sample)
		{
			SampleIndex = sampleIndex;
			WrongValue = wrongValue;
			Sample = sample;
		}
	}

	public sealed class TrainingRate
	{
		public int Epochs { get; set; }
		public double Rate { get; set; }
		public double MinimumRate { get; set; }
		public double WeightDecayFactor { get; set; }
        public double Momentum { get; set; }
		public int BatchSize { get; set; }
		public double InitialAvgLoss { get; set; }
		public int DecayAfterEpochs { get; set; }
		public double DecayFactor { get; set; }
		public double WeightSaveTreshold { get; set; }
		public bool Distorted {get; set;}
		public int DistortionPercentage { get; set; }
		public double SeverityFactor { get; set; }
		public double MaxScaling { get; set; }
		public double MaxRotation { get; set; }
		public double ElasticSigma { get; set; }
		public double ElasticScaling { get; set; }

		public TrainingRate()
		{
			Epochs = Properties.Settings.Default.Epochs;
			Rate = Properties.Settings.Default.Rate;
			MinimumRate = Properties.Settings.Default.MinimumRate;
			WeightDecayFactor = Properties.Settings.Default.WeightDecayFactor;
            Momentum = Properties.Settings.Default.Momentum;
			BatchSize = Properties.Settings.Default.BatchSize;
			InitialAvgLoss = Properties.Settings.Default.InitialAvgLoss;
			DecayAfterEpochs = Properties.Settings.Default.DecayAfterEpochs;
			DecayFactor = Properties.Settings.Default.DecayFactor;
			WeightSaveTreshold = Properties.Settings.Default.WeightSaveTreshold;
			Distorted = Properties.Settings.Default.Distorted;
			DistortionPercentage = Properties.Settings.Default.DistortionPercentage;
			SeverityFactor = Properties.Settings.Default.SeverityFactor;
			MaxScaling = Properties.Settings.Default.MaxScaling;
			MaxRotation = Properties.Settings.Default.MaxRotation;
			ElasticSigma = Properties.Settings.Default.ElasticSigma;
			ElasticScaling = Properties.Settings.Default.ElasticScaling;
		}

		public TrainingRate(double rate, int epochs, double minRate, double weightDecayFactor, double momentum, int batchSize, double initialAvgLoss, double decayFactor, int decayAfterEpochs, double weightSaveTreshold)
		{
			Epochs = epochs;
			Rate = rate;
			MinimumRate = minRate;
			WeightDecayFactor = weightDecayFactor;
            Momentum = momentum;
			BatchSize = batchSize;
			InitialAvgLoss = initialAvgLoss;
			DecayAfterEpochs = decayAfterEpochs;
			DecayFactor = decayFactor;
			WeightSaveTreshold = weightSaveTreshold;
			Distorted = false;
			DistortionPercentage = Properties.Settings.Default.DistortionPercentage;
			SeverityFactor = Properties.Settings.Default.SeverityFactor;
			MaxScaling = Properties.Settings.Default.MaxScaling;
			MaxRotation = Properties.Settings.Default.MaxRotation;
			ElasticSigma = Properties.Settings.Default.ElasticSigma;
			ElasticScaling = Properties.Settings.Default.ElasticScaling;
		}

		public TrainingRate(bool useDistortions)
		{
			Rate = Properties.Settings.Default.Rate;
			Epochs = Properties.Settings.Default.Epochs;
			WeightDecayFactor = Properties.Settings.Default.WeightDecayFactor;
            Momentum = Properties.Settings.Default.Momentum;
			BatchSize = Properties.Settings.Default.BatchSize;
			DecayAfterEpochs = Properties.Settings.Default.DecayAfterEpochs;
			DecayFactor = Properties.Settings.Default.DecayFactor;
			MinimumRate = Properties.Settings.Default.MinimumRate;
			InitialAvgLoss = Properties.Settings.Default.InitialAvgLoss;
			WeightSaveTreshold = Properties.Settings.Default.WeightSaveTreshold;
			Distorted = useDistortions;
			DistortionPercentage = Properties.Settings.Default.DistortionPercentage;
			SeverityFactor = Properties.Settings.Default.SeverityFactor;
			MaxScaling = Properties.Settings.Default.MaxScaling;
			MaxRotation = Properties.Settings.Default.MaxRotation;
			ElasticSigma = Properties.Settings.Default.ElasticSigma;
			ElasticScaling = Properties.Settings.Default.ElasticScaling;
		}

		public TrainingRate(double rate, int epochs, double minRate, double weightDecayFactor, double momentum, int batchSize, double initialAvgLoss, double decayFactor, int decayAfterEpochs, double weightSaveTreshold, bool distorted, int distortionPercentage, double severityFactor, double maxScaling, double maxRotation, double elasticSigma, double elasticScaling)
		{
			Rate = rate;
			Epochs = epochs;
			MinimumRate = minRate;
            WeightDecayFactor = weightDecayFactor;
            Momentum = momentum;
			BatchSize = batchSize;
			InitialAvgLoss = initialAvgLoss;
			DecayFactor = decayFactor;
			DecayAfterEpochs = decayAfterEpochs;
			WeightSaveTreshold = weightSaveTreshold;
			Distorted = distorted;
			DistortionPercentage = distortionPercentage;
			SeverityFactor = severityFactor;
			MaxScaling = maxScaling;
			MaxRotation = maxRotation;
			ElasticSigma = elasticSigma;
			ElasticScaling = elasticScaling;
		}
	}

	public class NeuralNetwork: IDisposable 
	{
		public delegate void LossFunctionActionDelegate(int correctClass);
		public delegate double GetSampleLossDelegate(int correctClass);
		public delegate int RecognizedDelegate();
		
		public DataProvider DataProvider { get; private set; }
		public int NetworkIndex { get; private set; }
		public string Name { get; private set; }
		public int ClassCount { get; private set; }
		public double TrainToValue { get; private set; }
		public LossFunctions LossFunction { get; private set; }
		public DataProviderSets DataProviderSet { get; private set; }
		public TrainingStrategy TrainingStrategy { get; set; }
        public bool SubstractMean { get; private set; }
		public double dMicron { get; private set; }
		public bool DropOutUsed { get; set; }
		public Layer LastLayer { get; private set; }
		public NetworkStates OperationState { get; set; }
		public NetworkStates OldOperationState { get; private set; }
		public LossFunctionActionDelegate LossFunctionAction;
		public GetSampleLossDelegate GetSampleLoss;
		public RecognizedDelegate Recognized;
		
		private string description = String.Empty;
		public string Description
		{
			get
			{
				description = String.Empty;
				description += "Name CNN:     \t" + Name + "\r\n";
				description += "Dataset Used: \t" + DataProviderSet.ToString() + "\r\n";
				description += "Loss Function:\t" + LossFunction.ToString() + "\r\n";
				foreach (Layer layer in this.Layers)
				{
					description += layer.Name;
					if (layer.NextLayer != null)
						description += "\r\n";
				}
				return description;
			}
		}

		public int ParallelTreshold = Properties.Settings.Default.ParallelTreshold;
		public TrainingRate TrainingRate { get; set; }
		public Layer[] Layers;
		public List<TrainingRate> TrainingRates { get; private set; }
		public List<List<Byte>> RbfWeights {get; private set;}
		public ThreadSafeRandom RandomGenerator;
		public ParallelOptions ParallelOption;
		
		private int maxDegreeOfParallelism = Environment.ProcessorCount;
		public int MaxDegreeOfParallelism
		{
			get
			{
				return maxDegreeOfParallelism;
			}
			set
			{
				if (value == maxDegreeOfParallelism)
					return;

				maxDegreeOfParallelism = value;
				ParallelOption.MaxDegreeOfParallelism = maxDegreeOfParallelism;
			}
		}

		public TestingParameters TestParameters { get; set; }
		public ImageData CurrentSample { get; set; }
		public TaskState CurrentTaskState = TaskState.Stopped;
		public Stopwatch TaskDuration;
		public Stopwatch SampleSpeedTimer;
		public TimeSpan EpochDuration;
		public int TotalEpochs;
		public int CurrentEpoch;
		public int SampleIndex;
		public int TrainErrors;
		public int TestErrors;
		public int HessianSamples = Properties.Settings.Default.HessianSamplesUsed;
		public double AvgTrainLoss;
		public double AvgTestLoss;
		public double SampleRate;
		public event EventHandler<EventArgs> RaiseNetworkProgressEvent = delegate { };
		public event EventHandler<AddUnrecognizedTestSampleEventArgs> RaiseAddUnrecognizedTestSampleEvent = delegate { };
        public double Spread { get; private set; }
        private double min = -1D;
		public double Min
		{
			get { return min; }
			set
			{
				min = value;
				Spread = Math.Abs(max) - min;
			} 
		}
        private double max = 1D;
        public double Max 
		{
			get { return max; }
			set
			{
				max = value;
				Spread = Math.Abs(max) - Min;
			} 
		}

		private StringBuilder timeStringBuilder;
		private static Thread workerThread;
		private System.Timers.Timer workerTimer;
		private object[] emptyParams;
		
		protected virtual void OnRaiseProgressEvent(EventArgs e)
		{
			EventHandler<EventArgs> handler = null;
			lock (this)
			{
				handler = RaiseNetworkProgressEvent;
			}
			if (handler != null)
				foreach (EventHandler<EventArgs> _handler in handler.GetInvocationList())
					System.Windows.Application.Current.Dispatcher.Invoke(_handler, new object[] { this, e });
		}

		protected virtual void OnRaiseTestErrorEvent(AddUnrecognizedTestSampleEventArgs e)
		{
			EventHandler<AddUnrecognizedTestSampleEventArgs> handler = null;
			lock (this)
			{
				handler = RaiseAddUnrecognizedTestSampleEvent;
			}

			if (handler != null)
				foreach (EventHandler<AddUnrecognizedTestSampleEventArgs> _handler in handler.GetInvocationList())
					System.Windows.Application.Current.Dispatcher.Invoke(_handler, new object[] { this, e });  
		}

		public NeuralNetwork(DataProvider dataProvider, string name = "Neural Network", int classCount = 10, double trainTo = 0.8D, LossFunctions lossFunction = LossFunctions.MeanSquareError, DataProviderSets dataProviderSet = DataProviderSets.MNIST, TrainingStrategy trainingStrategy = Model.TrainingStrategy.SGDLevenbergMarquardt, double dmicron = 0.02D, int hessianSamples = 2000, bool substractMean = false, double min= -1.0, double max = 1.0)
		{
            if (dataProvider == null)
                throw new ArgumentException("DataProvider is null");
            if ((dataProviderSet == DataProviderSets.MNIST) && (substractMean))
                throw new ArgumentException("SubstractMean parameter is true (should be false)");

            DataProvider = dataProvider;
			NetworkIndex = 0;
			Name = name.Trim();
			ClassCount = classCount;
			TrainToValue = trainTo;
			LossFunction = lossFunction;
			switch (LossFunction)
			{
				case LossFunctions.MeanSquareError:
					LossFunctionAction = MeanSquareErrorLossFunction;
					GetSampleLoss = GetSampleLossMSE;
					Recognized = ArgMax;
					break;

				case LossFunctions.CrossEntropy:
					LossFunctionAction = CrossEntropyLossFunction;
					GetSampleLoss = GetSampleLossCrossEntropy;
					Recognized = ArgMax;  
					break;
			}
		   
			DataProviderSet = dataProviderSet;
			TrainingStrategy = trainingStrategy;
            SubstractMean = substractMean;
			dMicron = dmicron;
            Min = min;
            Max = max;
			DropOutUsed = false;
			
			Layers = new Layer[0];
			RbfWeights = GetRBFLeNet5WeightSamples();
		   
			RandomGenerator = new ThreadSafeRandom();
		   
			ParallelOption = new ParallelOptions();
			ParallelOption.TaskScheduler = null;
			ParallelOption.MaxDegreeOfParallelism = maxDegreeOfParallelism;

			SampleRate = 0D;
			TaskDuration = new Stopwatch();
			SampleSpeedTimer = new Stopwatch();
			timeStringBuilder = new StringBuilder();
			emptyParams = new object[] { this, EventArgs.Empty };
		}

		public void AddLayer(LayerTypes layerType, int mapCount, int mapWidth, int mapHeight)
		{
			Layer newLayer = new Layer(this, layerType, mapCount, mapWidth, mapHeight);

			Array.Resize(ref Layers, Layers.Length + 1);
			Layers[Layers.Length - 1] = newLayer;
		}

		public void AddLayer(LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight, int receptiveFieldWidth, int receptiveFieldHeight, int strideX = 1, int strideY = 1, int padX = 0, int padY = 0, Mappings mappings = null)
		{
			Layer newLayer = new Layer(this, layerType, activationFunction, mapCount, mapWidth, mapHeight, receptiveFieldWidth, receptiveFieldHeight, strideX, strideY, padX, padY, mappings);

			Array.Resize(ref Layers, Layers.Length + 1);
			Layers[Layers.Length - 1] = newLayer;
		}

		public void AddLayer(LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight, int receptiveFieldWidth, int receptiveFieldHeight, int strideX, int strideY, int padX, int padY, int dropOutPercentage)
		{
			Layer newLayer = new Layer(this, layerType, activationFunction, mapCount, mapWidth, mapHeight, receptiveFieldWidth, receptiveFieldHeight, strideX, strideY, padX, padY, true, dropOutPercentage);

			Array.Resize(ref Layers, Layers.Length + 1);
			Layers[Layers.Length - 1] = newLayer;
		}

		public void AddLayer(LayerTypes layerType, ActivationFunctions activationFunction, int neuronCount)
		{
			Layer newLayer = new Layer(this, layerType, activationFunction, neuronCount);
			
			Array.Resize(ref Layers, Layers.Length + 1);
			Layers[Layers.Length - 1] = newLayer;
		}

		~NeuralNetwork()
		{
			// In case the client forgets to call
			// Dispose , destructor will be invoked for
			Dispose(false);
		}
			  
		protected virtual void Dispose(bool disposing)
		{
			if (disposing)
			{
				// Free managed objects.
                //if (Layers != null)
                //    Array.Clear(Layers, 0, Layers.Length);
                //if (TrainingRates != null)
                //    TrainingRates.Clear();
                //if (RbfWeights != null)
                //    RbfWeights.Clear();
				if (RandomGenerator != null)
					RandomGenerator.Dispose();
				if (workerTimer != null)
					workerTimer.Dispose();
			   
                //if (TaskDuration != null)
                //    TaskDuration = null;
                //Layers = null;
                //TrainingRates = null;
                //RbfWeights = null;
			}
			// Free unmanaged objects
		}

		public void Dispose()
		{
			Dispose(true);
			// Ensure that the destructor is not called
			GC.SuppressFinalize(this);
		}
		
		public override string ToString()
		{
			return Description;
		}

		public void ChangeTrainingStrategy()
		{
			foreach (Layer layer in Layers)
				layer.ChangeTrainingStrategy();

            GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
            GC.WaitForPendingFinalizers();
            GC.Collect();
		}

		public void SaveDefinition(string fileName)
		{
			if (fileName.Contains("definition-xml"))
			{
				using (CNNDataSet ds = new CNNDataSet())
				{
					ds.RemotingFormat = SerializationFormat.Xml;
					ds.SchemaSerializationMode = SchemaSerializationMode.IncludeSchema;

					CNNDataSet.DataProviderSetsRow dataProviderSetRow;
					CNNDataSet.LossFunctionsRow lossFunctionRow;
					CNNDataSet.TrainingStrategiesRow trainingStrategyRow;
					CNNDataSet.LayerTypesRow layerTypeRow;
					CNNDataSet.ActivationFunctionsRow activationTypeRow;

					ds.DataProviderSets.BeginLoadData();
					foreach (string dataProviderSetName in Enum.GetNames(typeof(DataProviderSets)))
						ds.DataProviderSets.AddDataProviderSetsRow(dataProviderSetName);
					ds.DataProviderSets.EndLoadData();

					ds.LossFunctions.BeginLoadData();
					foreach (string lossFunctionName in Enum.GetNames(typeof(LossFunctions)))
						ds.LossFunctions.AddLossFunctionsRow(lossFunctionName);
					ds.LossFunctions.EndLoadData();

					ds.TrainingStrategies.BeginLoadData();
					foreach (string trainingStrategyName in Enum.GetNames(typeof(TrainingStrategy)))
						ds.TrainingStrategies.AddTrainingStrategiesRow(trainingStrategyName);
					ds.TrainingStrategies.EndLoadData();

					ds.LayerTypes.BeginLoadData();
					foreach (string layerTypeName in Enum.GetNames(typeof(LayerTypes)))
						ds.LayerTypes.AddLayerTypesRow(layerTypeName);
					ds.LayerTypes.EndLoadData();

					ds.ActivationFunctions.BeginLoadData();
					foreach (string activationTypeName in Enum.GetNames(typeof(ActivationFunctions)))
						ds.ActivationFunctions.AddActivationFunctionsRow(activationTypeName);
					ds.ActivationFunctions.EndLoadData();

					dataProviderSetRow = ds.DataProviderSets.FindByDataProviderSet((int)DataProviderSet);
					lossFunctionRow = ds.LossFunctions.FindByLossFunction((int)LossFunction);
					trainingStrategyRow = ds.TrainingStrategies.FindByTrainingStrategy((int)TrainingStrategy);

					ds.NeuralNetworks.BeginLoadData();
					CNNDataSet.NeuralNetworksRow networkRow = ds.NeuralNetworks.AddNeuralNetworksRow((byte)NetworkIndex, Name, Description, ClassCount, dataProviderSetRow, trainingStrategyRow, dMicron, TrainToValue, lossFunctionRow, DropOutUsed, SubstractMean, HessianSamples, Min, Max);
					ds.NeuralNetworks.EndLoadData();

					ds.Layers.BeginLoadData();
					foreach (Layer layer in Layers)
					{
						layerTypeRow = ds.LayerTypes.FindByLayerType((int)layer.LayerType);
						activationTypeRow = ds.ActivationFunctions.FindByActivationFunction((int)layer.ActivationFunctionId);

						ds.Layers.AddLayersRow(networkRow, (byte)layer.LayerIndex, layerTypeRow, activationTypeRow, layer.NeuronCount, layer.UseMapInfo, layer.MapCount, layer.MapWidth, layer.MapHeight, layer.ReceptiveFieldWidth, layer.ReceptiveFieldHeight, layer.StrideX, layer.StrideY, layer.PadX, layer.PadY, layer.IsFullyMapped, layer.LockedWeights, layer.UseDropOut, layer.DropOutPercentage);

						if (layer.Mappings != null)
						{
							ds.Mappings.BeginLoadData();
							int previousMapIndex = 0;
							int currentMapIndex = 0;
							foreach (bool mapped in layer.Mappings.Mapping)
							{
								ds.Mappings.AddMappingsRow((byte)NetworkIndex, (byte)layer.LayerIndex, previousMapIndex, currentMapIndex, mapped);

								currentMapIndex++;
								if (currentMapIndex >= layer.MapCount)
								{
									currentMapIndex = 0;
									previousMapIndex++;
								}
							}
							ds.Mappings.EndLoadData();
						}

						//if (includeWeights)
						//{
						//    if (layer.Weights != null)
						//    {
						//        ds.Weights.BeginLoadData();
						//        int weightIndex = 0;
						//        foreach (Weight weight in layer.Weights)
						//            ds.Weights.AddWeightsRow((byte)NetworkIndex, (byte)layer.LayerIndex, weightIndex++, weight.Value);
						//        ds.Weights.EndLoadData();
						//    }
						//}
					}
					ds.Layers.EndLoadData();

					ds.WriteXml(fileName);
				   
					//BinaryFormatter fmt = new BinaryFormatter();
					//fmt.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
					//using (FileStream outFile = File.Create(fileName))
					//{
					//    fmt.Serialize(outFile, ds);
					//}
				}

				GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
				GC.WaitForPendingFinalizers();
				GC.Collect();
			}
		}
		
		public void SaveWeights(string fileName)
		{
			int totalWeightsCount = 0;
			for (int l=1; l < Layers.Length; l++)
				if (Layers[l].HasWeights)
					totalWeightsCount += Layers[l].WeightCount;

			int indexSize = sizeof(byte);
			int weightSize = sizeof(double);
			int recordSize = indexSize + weightSize;
			int fileSize = totalWeightsCount * recordSize;
			byte[] info = new byte[fileSize];
			int layerWeightsOffset = 0;
			for (int l = 1; l < Layers.Length; l++)
			{
				if (Layers[l].HasWeights)
				{
					if (Layers[l].UseWeightPartitioner)
					{
						Parallel.For(0, Layers[l].WeightCount, ParallelOption, i =>
						{
							int idx = layerWeightsOffset + (i * recordSize);
							info[idx] = (byte)l;
							idx += indexSize;
							byte[] temp = new byte[weightSize];
							temp = BitConverter.GetBytes(Layers[l].Weights[i].Value);
							for (int j = 0; j < weightSize; j++)
								info[idx + j] = temp[j];
						});
					}
					else
					{
						int idx;
						for (int i = 0; i < Layers[l].WeightCount; i++)
						{
							idx = layerWeightsOffset + (i * recordSize);
							info[idx] = (byte)l;
							idx += indexSize;
							byte[] temp = new byte[weightSize];
							temp = BitConverter.GetBytes(Layers[l].Weights[i].Value);
							for (int j = 0; j < weightSize; j++)
								info[idx + j] = temp[j];
						}
					}

					layerWeightsOffset += Layers[l].WeightCount * recordSize;
				}
			}

			using (FileStream outFile = File.Create(fileName, fileSize, FileOptions.RandomAccess))
			{
				outFile.Write(info, 0, fileSize);
				outFile.Flush();
			}
			info = null;

			GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
			GC.WaitForPendingFinalizers();
			GC.Collect();
		}

		public void SaveFullyZipped(string fileName)
		{
			string tempPath = Path.GetTempPath();
			string definitionFileName = Name + ".definition-xml";
			string weightsFileName = Name + ".weights-bin";

			SaveDefinition(tempPath + @"\" + definitionFileName);
			SaveWeights(tempPath + @"\" + weightsFileName);
			using (ZipArchive archive = ZipFile.Open(fileName, ZipArchiveMode.Update))
			{
				archive.CreateEntryFromFile(tempPath + @"\" + definitionFileName, definitionFileName + "-gz", CompressionLevel.Optimal);
				archive.CreateEntryFromFile(tempPath + @"\" + weightsFileName, weightsFileName + "-gz", CompressionLevel.Optimal); 
			}
			File.Delete(tempPath + @"\" + definitionFileName);
			File.Delete(tempPath + @"\" + weightsFileName);
		}

		public static NeuralNetwork LoadDefinition(string fileName, DataProvider dataProvider)
		{
			NeuralNetwork network = null;

			if (fileName.Contains("definition-xml"))
			{
				CNNDataSet ds = new CNNDataSet();
				ds.RemotingFormat = SerializationFormat.Xml;
				ds.ReadXml(fileName, XmlReadMode.Auto);


				if (ds != null && ds.NeuralNetworks.Count > 0)
				{
					CNNDataSet.NeuralNetworksRow networkRow = ds.NeuralNetworks.Rows[0] as CNNDataSet.NeuralNetworksRow;

					network = new NeuralNetwork(dataProvider, networkRow.Name, networkRow.ClassCount, networkRow.TrainToValue, (LossFunctions)networkRow.LossFunction, (DataProviderSets)networkRow.DataProviderSet, (TrainingStrategy)networkRow.TrainingStrategy, networkRow.DMicron, networkRow.HessianSamples, networkRow.SubstractMean);
                    network.Min = networkRow.Min;
                    network.Max = networkRow.Max;

					Layer layer = null;
					Layer previousLayer = null;
					CNNDataSet.LayersRow[] layerRows = networkRow.GetLayersRows();
					foreach (CNNDataSet.LayersRow layerRow in layerRows)
					{
						bool[] mapping = new bool[0];
						CNNDataSet.MappingsRow[] mappingRows = layerRow.GetMappingsRows();
						foreach (CNNDataSet.MappingsRow mappingRow in mappingRows)
						{
							Array.Resize(ref mapping, mapping.Length + 1);
							mapping[mapping.Length - 1] = mappingRow.IsMapped;
						}

						Mappings mappings = null;
						if (mapping.Length > 0)
							mappings = new Mappings(mapping);

						LayerTypes l = (LayerTypes)layerRow.LayerType;
						ActivationFunctions f = (ActivationFunctions)layerRow.ActivationFunction;

						layer = new Layer(network, layerRow.LayerIndex, l, f, layerRow.NeuronCount, layerRow.UseMapInfo, layerRow.MapCount, layerRow.MapWidth, layerRow.MapHeight, layerRow.IsFullyMapped, layerRow.ReceptiveFieldWidth, layerRow.ReceptiveFieldHeight, layerRow.StrideX, layerRow.StrideY, layerRow.PadX, layerRow.PadY, previousLayer, mappings, layerRow.LockedWeights, layerRow.UseDropOut, layerRow.DropOutPercentage);
						//if (includeWeights)
						//{
						//    CNNDataSet.WeightsRow[] rows = layerRow.GetWeightsRows();
						//    if (rows.Count() > 0)
						//        foreach (CNNDataSet.WeightsRow weightRow in rows)
						//            layer.Weights[weightRow.WeightIndex].Value = weightRow.Value;
						//}
						Array.Resize(ref network.Layers, network.Layers.Length + 1);
						network.Layers[network.Layers.Length - 1] = layer;
						
						previousLayer = layer;
					}

					network.LastLayer = layer;

					network.InitializeWeights();

					layer = null;
					previousLayer = null;
					layerRows = null;
				}
				else
				{
					TaskDialogOptions config = new TaskDialogOptions();
					config.MainIcon = VistaTaskDialogIcon.Information;
					config.MainInstruction = "Invalid data format.";
					config.Title = "Information";
					config.CommandButtons = new string[] { "&OK" };
					config.AllowDialogCancellation = false;
					TaskDialog.Show(config);
					network = null;
				}

				if (ds != null)
				{
					ds.Dispose();
					ds = null;
				}

				GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
				GC.WaitForPendingFinalizers();
				GC.Collect();
			}

			return network;
		}

		public bool LoadWeights(string fileName)
		{
			byte[] buffer = File.ReadAllBytes(fileName);

			int indexSize = sizeof(byte);
			int weightSize = sizeof(double);
			int recordSize = indexSize + weightSize;
			int totalWeightCount = buffer.Length / recordSize;
			int checkWeightCount = 0;
			for (int l = 1; l < Layers.Length; l++)
				if (Layers[l].HasWeights)
					checkWeightCount += Layers[l].WeightCount;

            if (totalWeightCount == checkWeightCount)
            {
                byte oldLayerIdx = 0;
                int weightIdx = 0;
                for (int index = 0; index < buffer.Length; index += recordSize)
                {
                    if (buffer[index] != oldLayerIdx)
                    {
                        weightIdx = 0;
                        oldLayerIdx = buffer[index];
                    }
                    byte[] temp = new byte[weightSize];
                    for (int j = 0; j < weightSize; j++)
                        temp[j] = buffer[index + j + indexSize];
                    Layers[buffer[index]].Weights[weightIdx++].Value = BitConverter.ToDouble(temp, 0);
                }
            }
            else
                return false;

			buffer = null;
			GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
			GC.WaitForPendingFinalizers();
			GC.Collect();

			return true;
		}

		public bool LoadWeightsXmlBin(string fileName)  // for importing and converting old style weight files
		{
			CNNDataSet ds = null;
			System.Runtime.Serialization.Formatters.Binary.BinaryFormatter fmt = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
			fmt.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
			
			try
			{
				using (FileStream inFile = File.OpenRead(fileName))
				{
					ds = (CNNDataSet)fmt.Deserialize(inFile);
				}
			}
			catch (Exception)
			{
				return false;
			}

			if (ds != null)
			{
				CNNDataSet.NeuralNetworksRow networkRow = (CNNDataSet.NeuralNetworksRow)ds.NeuralNetworks.Rows[0];
				if (Description.Contains(networkRow.Description.Trim()))  // check if we're dealing with an identical nertwork definition as the current one
				{
					CNNDataSet.LayersRow[] layerRows = networkRow.GetLayersRows();
					foreach (CNNDataSet.LayersRow row in layerRows)
					{
						CNNDataSet.WeightsRow[] weightRows = row.GetWeightsRows();
						foreach (CNNDataSet.WeightsRow weight in weightRows)
							Layers[row.LayerIndex].Weights[weight.WeightIndex].Value = weight.Value;
					}
				}
			}
			else
				return false;

			return true;
		}

		public static NeuralNetwork LoadFullyZipped(string fileName, DataProvider dataProvider)
		{
			string tempPath = Path.GetTempPath();
			string definitionFileName = String.Empty;
			string weightsFileName = String.Empty;

			bool correctFormat = false;
			using (ZipArchive archive = ZipFile.OpenRead(fileName))
			{
				if ((archive.Entries[0].Name.Contains(".definition-xml-gz")) && (archive.Entries[1].Name.Contains(".weights-bin-gz")))
				{
					definitionFileName = tempPath + @"\" + archive.Entries[0].Name.Replace("-gz", String.Empty);
					weightsFileName = tempPath + @"\" + archive.Entries[1].Name.Replace("-gz", String.Empty);
					archive.Entries[0].ExtractToFile(definitionFileName, true);
					archive.Entries[1].ExtractToFile(weightsFileName, true);
					correctFormat = true;
				}
				//archive.ExtractToDirectory(tempPath);
			}

			NeuralNetwork network = null;
			if (correctFormat)
			{
				network = LoadDefinition(definitionFileName, dataProvider);
				if (network != null)
					network.LoadWeights(weightsFileName);

				File.Delete(definitionFileName);
				File.Delete(weightsFileName);
			}

			return network;
		}

		public void AddGlobalTrainingRate(TrainingRate rate, Boolean clear = true)
		{
			 if (TrainingRates == null)
				 TrainingRates = new List<TrainingRate>();
			 else
				if (clear)
					TrainingRates.Clear();

			if (rate.DecayFactor == 1D)
				TrainingRates.Add(rate);
			else
			{
				// Decaying Training Rate
				if (rate.Epochs < rate.DecayAfterEpochs)
					rate.DecayAfterEpochs = rate.Epochs;

				int totIteration = rate.Epochs / rate.DecayAfterEpochs;
				double newRating = rate.Rate;
				
				for (int i = 0; i < totIteration; i++)
				{
					TrainingRates.Add(new TrainingRate(newRating, rate.DecayAfterEpochs, rate.MinimumRate, rate.WeightDecayFactor, rate.Momentum, rate.BatchSize, rate.InitialAvgLoss, 1, 1, rate.WeightSaveTreshold, rate.Distorted, rate.DistortionPercentage, rate.SeverityFactor, rate.MaxScaling, rate.MaxRotation, rate.ElasticSigma, rate.ElasticScaling));
					if (newRating * rate.DecayFactor > rate.MinimumRate)
						newRating *= rate.DecayFactor;
					else
						newRating = rate.MinimumRate;
				}

				if ((totIteration * rate.DecayAfterEpochs) < rate.Epochs)
					TrainingRates.Add(new TrainingRate(newRating, rate.Epochs - (totIteration * rate.DecayAfterEpochs), rate.MinimumRate, rate.WeightDecayFactor, rate.Momentum, rate.BatchSize, rate.InitialAvgLoss, 1, 1, rate.WeightSaveTreshold, rate.Distorted, rate.DistortionPercentage, rate.SeverityFactor, rate.MaxScaling, rate.MaxRotation, rate.ElasticSigma, rate.ElasticScaling));
			}
		}

		public void InitializeWeights()
		{
			for (int i=1; i < Layers.Length; i++)
			{
				if (Layers[i].HasWeights)
					Layers[i].InitializeWeights();

				// set the last layer in the network class!
				if (Layers[i].NextLayer == null)
					LastLayer = Layers[i];
			}
		}

		public BitmapSource BitmapToSource(System.Drawing.Image image)
		{
			MemoryStream stream = new MemoryStream();
			image.Save(stream, System.Drawing.Imaging.ImageFormat.Png);
			PngBitmapDecoder decoder = new PngBitmapDecoder(stream, BitmapCreateOptions.PreservePixelFormat, BitmapCacheOption.Default);

			BitmapSource destination = decoder.Frames[0];
			if (destination.CanFreeze)
				destination.Freeze();

			return destination;
		}

		public List<List<byte>> GetRBFLeNet5WeightSamples()
		{
			ResourceManager rm = new ResourceManager("CNNWB.Properties.Resources", Assembly.GetExecutingAssembly());
			CultureInfo ci = Thread.CurrentThread.CurrentCulture;

			List<List<byte>> rbfImages = new List<List<byte>>(10);
			int width = 7;
			int height = 12;

			for (int i = 0; i < 10; ++i)
			{
				System.Drawing.Image image = rm.GetObject("Image" + i.ToString(), ci) as System.Drawing.Image;
				BitmapSource bitmapSource = BitmapToSource(image);
				if (bitmapSource.CanFreeze)
					bitmapSource.Freeze();

				int rawStride = (width * bitmapSource.Format.BitsPerPixel + 7) / 8;
				byte[] rawImage = new byte[rawStride * height];
				bitmapSource.CopyPixels(rawImage, rawStride, 0);
				
				List<byte> inputPattern = new List<byte>(rawStride * 12);
				for (int j = 0; j < (rawStride * height); j++)
					inputPattern.Add(rawImage[j]);
				
				rbfImages.Add(inputPattern);
			}

			return rbfImages;
		}

		private int ArgMin()
		{
			int bestIndex = 0;
			double minValue = double.MaxValue;

			for (int i = 0; i < ClassCount; i++)
			{
				if (LastLayer.Neurons[i].Output < minValue)
				{
					minValue = LastLayer.Neurons[i].Output;
					bestIndex = i;
				}
			}

			return bestIndex;
		}

		private int ArgMax()
		{
			int bestIndex = 0;
			double maxValue = double.MinValue;

			for (int i = 0; i < ClassCount; i++)
			{
				if (LastLayer.Neurons[i].Output > maxValue)
				{
					maxValue = LastLayer.Neurons[i].Output;
					bestIndex = i;
				}
			}

			return bestIndex;
		}

		private double GetSampleLossMSE(int correctClass)
		{
			double patternLoss = 0D;

			for (int i = 0; i < ClassCount; i++)
			{
				if (i == correctClass)
					patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
				else
					patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
			}

			return patternLoss * 0.5D;
		}

		private double GetSampleLossCrossEntropy(int correctClass)
		{
			return -0.5D * MathUtil.Log(LastLayer.Neurons[correctClass].Output);
		}
		
		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public string GetTaskDuration()
		{
			timeStringBuilder.Length = 0;
			if (TaskDuration.Elapsed.Days > 0)
				if (TaskDuration.Elapsed.Days == 1)
					timeStringBuilder.AppendFormat("{0:D} day {1:D2}:{2:D2}:{3:D2}", TaskDuration.Elapsed.Days, TaskDuration.Elapsed.Hours, TaskDuration.Elapsed.Minutes, TaskDuration.Elapsed.Seconds);
				else
					timeStringBuilder.AppendFormat("{0:D} days {1:D2}:{2:D2}:{3:D2}", TaskDuration.Elapsed.Days, TaskDuration.Elapsed.Hours, TaskDuration.Elapsed.Minutes, TaskDuration.Elapsed.Seconds);
			else
				timeStringBuilder.AppendFormat("{0:D2}:{1:D2}:{2:D2}", TaskDuration.Elapsed.Hours, TaskDuration.Elapsed.Minutes, TaskDuration.Elapsed.Seconds);

			return timeStringBuilder.ToString();
		}

		private bool CheckStateChange()
		{
			if (CurrentTaskState != TaskState.Running)
			{
				if (CurrentTaskState == TaskState.Paused)
				{
					workerTimer.Stop();
					TaskDuration.Stop();
					SampleSpeedTimer.Stop();
					while (CurrentTaskState == TaskState.Paused)
						Thread.Sleep(100);
					TaskDuration.Start();
					workerTimer.Start();
					SampleSpeedTimer.Start();
				}
				else
				{
					TaskDuration.Stop();
					SampleSpeedTimer.Stop();
					return false;
				}
			}

			return true;
		}

		private void WorkerTimerElapsed(Object data, System.Timers.ElapsedEventArgs eventArgs)
		{
			double secs = SampleSpeedTimer.ElapsedMilliseconds * 0.001D;
			if (secs > 0.5D)
				SampleRate = (SampleIndex + 1) / secs;
			else
				SampleRate = 0D;

			//OnRaiseProgressEvent(EventArgs.Empty);
			System.Windows.Application.Current.Dispatcher.BeginInvoke(RaiseNetworkProgressEvent, emptyParams);
		}

		public void StartTraining()
		{
			if (CurrentTaskState == TaskState.Stopped)
			{
				//CurrentSample = DataProvider.TrainingSamples[0];

				workerThread = new Thread(new ThreadStart(TrainingTask));
				workerThread.IsBackground = true;
                                
				workerTimer = new System.Timers.Timer(1000);
				workerTimer.Elapsed += new ElapsedEventHandler(WorkerTimerElapsed);
				TaskDuration.Reset();
				SampleSpeedTimer.Reset();

				CurrentTaskState = TaskState.Running;
				EpochDuration = TimeSpan.Zero;
				workerThread.Start();
				while (!workerThread.IsAlive) Thread.SpinWait(1);

				TaskDuration.Start();
				SampleSpeedTimer.Start();
				workerTimer.Start();
			}
		}

		public void StopTraining()
		{
			if (CurrentTaskState != TaskState.Stopped)
			{
				CurrentTaskState = TaskState.Stopped;
                
				workerThread.Join();
                
                TaskDuration.Stop();
                SampleSpeedTimer.Stop();
                workerTimer.Stop();
				workerTimer.Elapsed -= new ElapsedEventHandler(WorkerTimerElapsed);
				workerTimer.Dispose();
                
                GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
                GC.WaitForPendingFinalizers();
                GC.Collect();
			}
		}

		public void TrainingTask()
		{
            bool calculatePseudoHessian = (TrainingStrategy == TrainingStrategy.SGDLevenbergMarquardt) || (TrainingStrategy == TrainingStrategy.SGDLevenbergMarquardtModA) || (TrainingStrategy == TrainingStrategy.MiniBatchSGDLevenbergMarquardt) || (TrainingStrategy == TrainingStrategy.MiniBatchSGDLevenbergMarquardtModA);
            double prevAvgLoss;
			string oldSaveWeightsFileName = String.Empty;
			int learningRateIndex = 0;
			int learningRateEpochs = TrainingRates[0].Epochs;
			int bestScore = (int)(((double)DataProvider.TestingSamplesCount / 100D) * (100D - TrainingRates[0].WeightSaveTreshold));
			int sampleSize = DataProvider.SampleSize * DataProvider.SampleChannels;
			int sampleLabel;
			int bestIndex = 0;
			double totLoss = 0D;
			double patternLoss;
			double maxValue;
			TotalEpochs = 0;
			foreach (TrainingRate rate in TrainingRates)
				TotalEpochs += rate.Epochs;

			TrainingRate = TrainingRates[0];
			CurrentEpoch = 0;
			prevAvgLoss = TrainingRate.InitialAvgLoss * 0.1D;
			while ((CurrentTaskState != TaskState.Stopped) && (CurrentEpoch < TotalEpochs))
			{
				if (CurrentEpoch == learningRateEpochs)
				{
					learningRateIndex++;
					TrainingRate = TrainingRates[learningRateIndex];
					learningRateEpochs += TrainingRate.Epochs;
				}

				CurrentEpoch++;

                if (calculatePseudoHessian)
				{
					SampleSpeedTimer.Restart();
					OperationState = NetworkStates.CalculatingHessian;
					
					DataProvider.ScrambleTrainingSamples();
					EraseHessian();

					double[] D2ErrX;
					for (SampleIndex = 0; SampleIndex < HessianSamples; SampleIndex++)
					{
                        if (SubstractMean)
                        {
                            double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                            for (int i = 0; i < sampleSize; i++)
                                Layers[0].Neurons[i].Output = data[i];
                        }
                        else
                            for (int i = 0; i < sampleSize; i++)
                                Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                        
						// fprop
						for (int i = 1; i < Layers.Length; i++)
							Layers[i].CalculateAction();

						// bbprop
						D2ErrX = new double[ClassCount];
						for (int i = 0; i < ClassCount; i++)
							D2ErrX[i] = 1D; // TrainToValue;  Is this also true with cross entropy loss and a softmax output layer ???

						for (int i = Layers.Length - 1; i > 1; i--)
							D2ErrX = Layers[i].BackpropagateSecondDerivativesAction(D2ErrX);

						if (CurrentTaskState != TaskState.Running)
							if (!CheckStateChange())
							{
								DivideHessianBy(SampleIndex + 1);
								break;
							}
					}

					DivideHessianBy(HessianSamples);
				}
			   
				SampleSpeedTimer.Restart();
				OperationState = NetworkStates.Training;
				DataProvider.ScrambleTrainingSamples();
				AvgTrainLoss = 0D;
				totLoss = 0D;
				TrainErrors = 0;
				
				int temp = (int)Math.Floor((double)DataProvider.TrainingSamplesCount / TrainingRate.BatchSize);
				int end = temp * TrainingRate.BatchSize; 
				int finalBatch = DataProvider.TrainingSamplesCount - end;
				
				switch (TrainingStrategy)
				{ 
					case TrainingStrategy.SGD:
					case TrainingStrategy.SGDLevenbergMarquardt:
                    case TrainingStrategy.SGDLevenbergMarquardtModA:
                    case TrainingStrategy.SGDM:
						switch (LossFunction)
						{
							case LossFunctions.MeanSquareError:
								if (TrainingRate.Distorted)
								{
									for (SampleIndex = 0; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
									{
                                        if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
                                        {
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                            {
                                                CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                            }
                                        }
                                        else
                                        {
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                                        }
										sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;
								
										// fprop
										for (int i = 1; i < Layers.Length; i++)
											Layers[i].CalculateAction();

										// Mean Square Error loss
										patternLoss = 0D;
										for (int i = 0; i < ClassCount; i++)
										{
											if (i == sampleLabel)
												patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
											else
												patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
										}
										patternLoss *= 0.5D;

										totLoss += patternLoss;
										AvgTrainLoss = totLoss / (SampleIndex + 1);

										bestIndex = 0;
										maxValue = LastLayer.Neurons[0].Output;
										for (int i = 1; i < ClassCount; i++)
										{
											if (LastLayer.Neurons[i].Output > maxValue)
											{
												maxValue = LastLayer.Neurons[i].Output;
												bestIndex = i;
											}
										}
								
										if (sampleLabel != bestIndex)
											TrainErrors++;

										if (patternLoss > prevAvgLoss)
										{
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

										   
											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
											{
												Layers[i].EraseGradientWeights();
												Layers[i].BackpropagateAction();
												Layers[i].UpdateWeights(1);
											}
										}

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								else
								{
									for (SampleIndex = 0; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
									{
                                        if (SubstractMean)
                                        {
                                            double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                            for (int i = 0; i < sampleSize; i++)
                                                Layers[0].Neurons[i].Output = data[i];
                                        }
                                        else
                                            for (int i = 0; i < sampleSize; i++)
                                                Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;

										sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;
								
										// fprop
										for (int i = 1; i < Layers.Length; i++)
											Layers[i].CalculateAction();

										// Mean Square Error loss
										patternLoss = 0D;
										for (int i = 0; i < ClassCount; i++)
										{
											if (i == sampleLabel)
												patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
											else
												patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
										}
										patternLoss *= 0.5D;

										totLoss += patternLoss;
										AvgTrainLoss = totLoss / (SampleIndex + 1);

										bestIndex = 0;
										maxValue = LastLayer.Neurons[0].Output;
										for (int i = 1; i < ClassCount; i++)
										{
											if (LastLayer.Neurons[i].Output > maxValue)
											{
												maxValue = LastLayer.Neurons[i].Output;
												bestIndex = i;
											}
										}
								
										if (sampleLabel != bestIndex)
											TrainErrors++;

										if (patternLoss > prevAvgLoss)
										{
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
											{
												Layers[i].EraseGradientWeights();
												Layers[i].BackpropagateAction();
												Layers[i].UpdateWeights(1);
											}
										}

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								break;

							case LossFunctions.CrossEntropy:
								if (TrainingRate.Distorted)
								{
									for (SampleIndex = 0; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
									{
                                        if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
                                        {
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                            {
                                                CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                            }
                                        }
                                        else
                                        {
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                                        }					
										sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;
							   
										// fprop
										for (int i = 1; i < Layers.Length; i++)
											Layers[i].CalculateAction();

										// Cross Entropy Loss
										patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

										totLoss += patternLoss;
										AvgTrainLoss = totLoss / (SampleIndex + 1);
								
										bestIndex = 0;
										maxValue = LastLayer.Neurons[0].Output;
										for (int i = 1; i < ClassCount; i++)
										{
											if (LastLayer.Neurons[i].Output > maxValue)
											{
												maxValue = LastLayer.Neurons[i].Output;
												bestIndex = i;
											}
										}

										if (sampleLabel != bestIndex)
											TrainErrors++;

										if (patternLoss > prevAvgLoss)
										{
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
											{
												Layers[i].EraseGradientWeights();
												Layers[i].BackpropagateAction();
												Layers[i].UpdateWeights(1);
											}
										}

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								else
								{
									for (SampleIndex = 0; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
									{
                                        if (SubstractMean)
                                        {
                                            double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                            for (int i = 0; i < sampleSize; i++)
                                                Layers[0].Neurons[i].Output = data[i];
                                        }
                                        else
                                            for (int i = 0; i < sampleSize; i++)
                                                Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                                        
                                        sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;
								
										// fprop
										for (int i = 1; i < Layers.Length; i++)
											Layers[i].CalculateAction();

										// Cross Entropy Loss
										patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

										totLoss += patternLoss;
										AvgTrainLoss = totLoss / (SampleIndex + 1);
								
										bestIndex = 0;
										maxValue = LastLayer.Neurons[0].Output;
										for (int i = 1; i < ClassCount; i++)
										{
											if (LastLayer.Neurons[i].Output > maxValue)
											{
												maxValue = LastLayer.Neurons[i].Output;
												bestIndex = i;
											}
										}

										if (sampleLabel != bestIndex)
											TrainErrors++;

										if (patternLoss > prevAvgLoss)
										{
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;
											
											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
											{
												Layers[i].EraseGradientWeights();
												Layers[i].BackpropagateAction();
												Layers[i].UpdateWeights(1);
											}
										}

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
							break;
						}
						break;

					case TrainingStrategy.MiniBatchSGD:
					case TrainingStrategy.MiniBatchSGDLevenbergMarquardt:
                    case TrainingStrategy.MiniBatchSGDLevenbergMarquardtModA:
                    case TrainingStrategy.MiniBatchSGDM:
						switch (LossFunction)
						{
							case LossFunctions.MeanSquareError:
								if (TrainingRate.Distorted)
								{
									for (SampleIndex = 0; SampleIndex < end; SampleIndex += TrainingRate.BatchSize)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (int batchIndex = 0; batchIndex < TrainingRate.BatchSize; batchIndex++)
										{
											if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
											{
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                {
                                                    CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                                }
											}
											else
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Image[i] / 255D) * Spread) + Min;
                                            }

											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Mean Square Error loss
											patternLoss = 0D;
											for (int i = 0; i < ClassCount; i++)
											{
												if (i == sampleLabel)
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
												else
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
											}
											patternLoss *= 0.5D;

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + batchIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										}

										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(TrainingRate.BatchSize);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}

									if (finalBatch > 0)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (SampleIndex = DataProvider.TrainingSamplesCount - finalBatch; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
										{
                                            if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
                                            {
                                                
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                {
                                                    CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                                }
                                            }
                                            else
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                                            }

											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Mean Square Error loss
											patternLoss = 0D;
											for (int i = 0; i < ClassCount; i++)
											{
												if (i == sampleLabel)
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
												else
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
											}
											patternLoss *= 0.5D;

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
											
										}

										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(finalBatch);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								else
								{
									for (SampleIndex = 0; SampleIndex < end; SampleIndex += TrainingRate.BatchSize)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (int batchIndex = 0; batchIndex < TrainingRate.BatchSize; batchIndex++)
										{
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Image[i] / 255D) * Spread) + Min;
											
                                            sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Mean Square Error loss
											patternLoss = 0D;
											for (int i = 0; i < ClassCount; i++)
											{
												if (i == sampleLabel)
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
												else
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
											}
											patternLoss *= 0.5D;

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + batchIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										   
										}
										
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(TrainingRate.BatchSize);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}

									if (finalBatch > 0)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (SampleIndex = DataProvider.TrainingSamplesCount-finalBatch; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex ++)
										{
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
											
                                            sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Mean Square Error loss
											patternLoss = 0D;
											for (int i = 0; i < ClassCount; i++)
											{
												if (i == sampleLabel)
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
												else
													patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
											}
											patternLoss *= 0.5D;

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										}
									}

									for (int i = Layers.Length - 1; i > 1; i--)
										Layers[i].UpdateWeights(finalBatch);

									if (CurrentTaskState != TaskState.Running)
										if (!CheckStateChange())
											break;
								}
								break;

							case LossFunctions.CrossEntropy:
								if (TrainingRate.Distorted)
								{
									for (SampleIndex = 0; SampleIndex < end; SampleIndex += TrainingRate.BatchSize)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (int batchIndex = 0; batchIndex < TrainingRate.BatchSize; batchIndex++)
										{
                                            if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                {
                                                    CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                                }

                                            }
                                            else
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Image[i] / 255D) * Spread) + Min;
                                            }

											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Cross Entropy Loss
											patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + batchIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;
											
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										}

										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(TrainingRate.BatchSize);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}

									if (finalBatch > 0)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (SampleIndex = DataProvider.TrainingSamplesCount-finalBatch; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex ++)
										{
                                            if (RandomGenerator.NextPercentage() < TrainingRate.DistortionPercentage)
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling).SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                {
                                                    CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
                                                }
                                            }
                                            else
                                            {
                                                if (SubstractMean)
                                                {
                                                    double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = data[i];
                                                }
                                                else
                                                    for (int i = 0; i < sampleSize; i++)
                                                        Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;
                                            }
											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Cross Entropy Loss
											patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;
											   
											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										}
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(finalBatch);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								else
								{
									for (SampleIndex = 0; SampleIndex < end; SampleIndex += TrainingRate.BatchSize)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (int batchIndex = 0; batchIndex < TrainingRate.BatchSize; batchIndex++)
										{
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Image[i] / 255D) * Spread) + Min;

											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex + batchIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Cross Entropy Loss
											patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + batchIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();
										}
										
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(TrainingRate.BatchSize);
										
										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}

									if (finalBatch > 0)
									{
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].EraseGradientWeights();

										for (SampleIndex = DataProvider.TrainingSamplesCount-finalBatch; SampleIndex < DataProvider.TrainingSamplesCount; SampleIndex++)
										{
                                            if (SubstractMean)
                                            {
                                                double[] data = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].SubstractMean(DataProvider);
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = data[i];
                                            }
                                            else
                                                for (int i = 0; i < sampleSize; i++)
                                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Image[i] / 255D) * Spread) + Min;

											sampleLabel = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[SampleIndex]].Label;

											// fprop
											for (int i = 1; i < Layers.Length; i++)
												Layers[i].CalculateAction();

											// Cross Entropy Loss
											patternLoss = -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);

											totLoss += patternLoss;
											AvgTrainLoss = totLoss / (SampleIndex + 1);

											bestIndex = 0;
											maxValue = LastLayer.Neurons[0].Output;
											for (int i = 1; i < ClassCount; i++)
											{
												if (LastLayer.Neurons[i].Output > maxValue)
												{
													maxValue = LastLayer.Neurons[i].Output;
													bestIndex = i;
												}
											}

											if (sampleLabel != bestIndex)
												TrainErrors++;

											for (int i = 0; i < ClassCount; i++)
												LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
											LastLayer.Neurons[sampleLabel].D1ErrX = LastLayer.Neurons[sampleLabel].Output - TrainToValue;

											// bprop
											for (int i = Layers.Length - 1; i > 1; i--)
												Layers[i].BackpropagateAction();

										}
										for (int i = Layers.Length - 1; i > 1; i--)
											Layers[i].UpdateWeights(finalBatch);

										if (CurrentTaskState != TaskState.Running)
											if (!CheckStateChange())
												break;
									}
								}
								break;
						}
						break;
				}
				
				if (CurrentTaskState != TaskState.Running)
					if (!CheckStateChange())
						break;

				// Calculate Test Error 
				OperationState = NetworkStates.CalculatingTestError;
				SampleSpeedTimer.Restart();

				double totTestLoss = 0;
				AvgTestLoss = 0D;
				TestErrors = 0;

				switch (LossFunction)
				{
					case LossFunctions.MeanSquareError:
						for (SampleIndex = 0; SampleIndex < DataProvider.TestingSamplesCount; SampleIndex++)
						{
                            if (SubstractMean)
                            {
                                double[] data = DataProvider.TestingSamples[SampleIndex].SubstractMean(DataProvider);
                                for (int i = 0; i < sampleSize; i++)
                                    Layers[0].Neurons[i].Output = data[i]; 
                            }
                            else
                                for (int i = 0; i < sampleSize; i++)
                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TestingSamples[SampleIndex].Image[i] / 255D) * Spread) + Min;

							sampleLabel = DataProvider.TestingSamples[SampleIndex].Label;
							
							// fprop
							for (int i = 1; i < Layers.Length; i++)
								Layers[i].CalculateAction();
							
							// Mean Square Error loss
							patternLoss = 0D;
							for (int i = 0; i < ClassCount; i++)
							{
								if (i == sampleLabel)
									patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output - TrainToValue);
								else
									patternLoss += MathUtil.Pow2(LastLayer.Neurons[i].Output + TrainToValue);
							}
							patternLoss *= 0.5D;
							totTestLoss += patternLoss;
							AvgTestLoss = totTestLoss / (SampleIndex + 1);

							bestIndex = 0;
							maxValue = LastLayer.Neurons[0].Output;
							for (int i = 1; i < ClassCount; i++)
							{
								if (LastLayer.Neurons[i].Output > maxValue)
								{
									maxValue = LastLayer.Neurons[i].Output;
									bestIndex = i;
								}
							}

							if (sampleLabel != bestIndex)
								TestErrors++;

							if (CurrentTaskState != TaskState.Running)
								if (!CheckStateChange())
									break;
						}
						break;

					case LossFunctions.CrossEntropy:
						for (SampleIndex = 0; SampleIndex < DataProvider.TestingSamplesCount; SampleIndex++)
						{
                            if (SubstractMean)
                            {
                                double[] data = DataProvider.TestingSamples[SampleIndex].SubstractMean(DataProvider);
                                for (int i = 0; i < sampleSize; i++)
                                    Layers[0].Neurons[i].Output = data[i];
                            }
                            else
                                for (int i = 0; i < sampleSize; i++)
                                    Layers[0].Neurons[i].Output = (((double)DataProvider.TestingSamples[SampleIndex].Image[i] / 255D) * Spread) + Min;

							sampleLabel = DataProvider.TestingSamples[SampleIndex].Label;
							
							// fprop
							for (int i = 1; i < Layers.Length; i++)
								Layers[i].CalculateAction();

							// Cross Entropy Loss
							totTestLoss += -0.5D * MathUtil.Log(LastLayer.Neurons[sampleLabel].Output);
							AvgTestLoss = totTestLoss / (SampleIndex + 1);

							bestIndex = 0;
							maxValue = LastLayer.Neurons[0].Output;
							for (int i = 1; i < ClassCount; i++)
							{
								if (LastLayer.Neurons[i].Output > maxValue)
								{
									maxValue = LastLayer.Neurons[i].Output;
									bestIndex = i;
								}
							}

							if (sampleLabel != bestIndex)
								TestErrors++;

							if (CurrentTaskState != TaskState.Running)
								if (!CheckStateChange())
									break;
						}
						break;
				}

				// Save best score
				if (TestErrors <= bestScore)
				{
					bestScore = TestErrors;
					string fileName = DataProvider.StorageDirectory + @"\" + Name + " (epoch " + CurrentEpoch.ToString() + " - " + bestScore.ToString() + " errors).weights-bin";
					OperationState = NetworkStates.SavingWeights;
					SaveWeights(fileName);
					if ((oldSaveWeightsFileName != String.Empty) && (File.Exists(oldSaveWeightsFileName)))
						File.Delete(oldSaveWeightsFileName);
					oldSaveWeightsFileName = fileName;
				}
				
				if (CurrentTaskState != TaskState.Running)
					if (!CheckStateChange())
						break;

				OperationState = NetworkStates.NewEpoch;
				prevAvgLoss = AvgTrainLoss * 0.1D;
			   
				while (OperationState == NetworkStates.NewEpoch)
				{
					Thread.SpinWait(1);
					if (CurrentTaskState == TaskState.Stopped)
						break;
				}
			}
			CurrentTaskState = TaskState.Stopped;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private int SetTrainingInputSample(int index)
		{
			if (TrainingRate.Distorted)
			{
				if (RandomGenerator.Next(100) <= TrainingRate.DistortionPercentage)
					CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[index]].Distorted(DataProvider, TrainingRate.SeverityFactor, TrainingRate.MaxScaling, TrainingRate.MaxRotation, TrainingRate.ElasticSigma, TrainingRate.ElasticScaling);
				else
					CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[index]];
			}
			else
				CurrentSample = DataProvider.TrainingSamples[DataProvider.RandomTrainingSample[index]];

            if (SubstractMean)
            {
                double[] data = CurrentSample.SubstractMean(DataProvider);
                for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                    Layers[0].Neurons[i].Output = data[i];
            }
            else
                for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                    Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
			
			return CurrentSample.Label;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private int SetTestingInputSample(int index)
		{
			CurrentSample = DataProvider.TestingSamples[index];

            if (SubstractMean)
            {
                double[] data = CurrentSample.SubstractMean(DataProvider);
                for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                    Layers[0].Neurons[i].Output = data[i];
            }
            else
               for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;
			
			return CurrentSample.Label;
		}

		public void StartTesting()
		{
			if (CurrentTaskState == TaskState.Stopped)
			{
                
				workerThread = new Thread(new ThreadStart(TestingTask));
				workerThread.IsBackground = true;
                
				workerTimer = new System.Timers.Timer(1000);
				workerTimer.Elapsed += new ElapsedEventHandler(WorkerTimerElapsed);
				TaskDuration.Reset();
				SampleSpeedTimer.Reset();

				CurrentTaskState = TaskState.Running;
				workerThread.Start();
				while (!workerThread.IsAlive) Thread.SpinWait(1);

				TaskDuration.Start();
				SampleSpeedTimer.Start();
				workerTimer.Start();
			}
		}

		public void StopTesting()
		{
			if (CurrentTaskState != TaskState.Stopped)
			{
				CurrentTaskState = TaskState.Stopped;
								
				workerThread.Join(250);
                if (workerThread.IsAlive)
                    workerThread.Abort();

				TaskDuration.Stop();
                SampleSpeedTimer.Stop();
                workerTimer.Stop();
				workerTimer.Elapsed -= new ElapsedEventHandler(WorkerTimerElapsed);
				workerTimer.Dispose();

                GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
                GC.WaitForPendingFinalizers();
                GC.Collect();
			}
		}

		public void TestingTask()
		{
			double totLoss = 0D;
			int bestIndex = 0;
			int totalSamples;
			if (TestParameters.UseTrainingSamples)
				totalSamples = DataProvider.TrainingSamplesCount;
			else
				totalSamples = DataProvider.TestingSamplesCount;
			
			OperationState = NetworkStates.Testing;
			AvgTestLoss = 0D;
			TestErrors = 0;
			SampleSpeedTimer.Restart();
			for (SampleIndex = 0; SampleIndex < totalSamples; SampleIndex++)
			{
				if (TestParameters.UseTrainingSamples)
				{
					if (TestParameters.Distorted)
						CurrentSample = DataProvider.TrainingSamples[SampleIndex].Distorted(DataProvider, TestParameters.SeverityFactor, TestParameters.MaxScaling, TestParameters.MaxRotation, TestParameters.ElasticSigma, TestParameters.ElasticScaling);
					else
						CurrentSample = DataProvider.TrainingSamples[SampleIndex];
				}
				else
				{
				   if (TestParameters.Distorted)
						CurrentSample = DataProvider.TestingSamples[SampleIndex].Distorted(DataProvider, TestParameters.SeverityFactor, TestParameters.MaxScaling, TestParameters.MaxRotation, TestParameters.ElasticSigma, TestParameters.ElasticScaling);
					else
						CurrentSample = DataProvider.TestingSamples[SampleIndex];
				}

                if (SubstractMean)
                {
                    double[] data = CurrentSample.SubstractMean(DataProvider);
                    for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                        Layers[0].Neurons[i].Output = data[i];
                }
                else
                    for (int i = 0; i < DataProvider.SampleSize * DataProvider.SampleChannels; i++)
                        Layers[0].Neurons[i].Output = (((double)CurrentSample.Image[i] / 255D) * Spread) + Min;

				// fprop
				for (int i = 1; i < Layers.Length; i++)
					Layers[i].CalculateAction();

				totLoss += GetSampleLoss(CurrentSample.Label);
				AvgTestLoss = totLoss / (SampleIndex + 1);
				bestIndex = Recognized();

				if (bestIndex != CurrentSample.Label)
				{
					TestErrors++;
					System.Windows.Application.Current.Dispatcher.Invoke(RaiseAddUnrecognizedTestSampleEvent, System.Windows.Threading.DispatcherPriority.Send, new object[] { null, new AddUnrecognizedTestSampleEventArgs(SampleIndex, bestIndex, CurrentSample) });  
					//OnRaiseTestErrorEvent(new AddUnrecognizedTestSampleEventArgs(SampleIndex, bestIndex, CurrentSample));
				}

				if (CurrentTaskState != TaskState.Running)
					if (!CheckStateChange())
						break;
			}
			SampleSpeedTimer.Stop();
			CurrentTaskState = TaskState.Stopped;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public void Calculate()
		{
			for (int i = 1; i < Layers.Length; i++)
				Layers[i].CalculateAction();
		}

		public void MeanSquareErrorLossFunction(int correctClass)
		{
			for (int i = 0; i < ClassCount; i++)
				LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output + TrainToValue;
			LastLayer.Neurons[correctClass].D1ErrX = LastLayer.Neurons[correctClass].Output - TrainToValue;
		}

		public void CrossEntropyLossFunction(int correctClass)
		{
			for (int i = 0; i < ClassCount; i++)
				LastLayer.Neurons[i].D1ErrX = LastLayer.Neurons[i].Output;
			LastLayer.Neurons[correctClass].D1ErrX = LastLayer.Neurons[correctClass].Output - TrainToValue;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public void Backpropagate(int correctClass)
		{
			LossFunctionAction(correctClass);

			for (int i = Layers.Length - 1; i > 1; i--)
				Layers[i].BackpropagateAction();
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public void BackpropagateSecondDerivatives()
		{
			// start the process by calculating the second derivative for the last layer.
			// for the standard MSE Err function (i.e., 0.5*sumof( (desiredOutput)^2 ), this differential is 
			// exactly one

			double[] D2ErrX = new double[ClassCount];
			for (int i = 0; i < ClassCount; i++)
				D2ErrX[i] = 1D; // TrainToValue;  Is this also true with cross entropy loss and a softmax output layer ???

			for (int i = Layers.Length - 1; i > 1; i--)
				D2ErrX = Layers[i].BackpropagateSecondDerivativesAction(D2ErrX);
		}

		public void EraseHessian()
		{
			for (int i=1; i < Layers.Length; i++)
			   if (Layers[i].HasWeights)
				{
					if (Layers[i].UseWeightPartitioner)
						Parallel.For(0, Layers[i].WeightCount, ParallelOption, j => Layers[i].Weights[j].DiagonalHessian = 0D);
					else
						for (int j = 0; j < Layers[i].WeightCount; j++)
							Layers[i].Weights[j].DiagonalHessian = 0D;
				}
		}

		public void DivideHessianBy(double divisor)
		{
			for (int i = 1; i < Layers.Length; i++)
				if (Layers[i].HasWeights)
				{
					if (Layers[i].UseWeightPartitioner)
						Parallel.For(0, Layers[i].WeightCount, ParallelOption, j => Layers[i].Weights[j].DiagonalHessian /= divisor);
					else
						for (int j = 0; j < Layers[i].WeightCount; j++)
							Layers[i].Weights[j].DiagonalHessian /= divisor;
				}
		}
	}
 
	public sealed class Layer
	{
		public NeuralNetwork Network;
		public string Name;
		public LayerTypes LayerType;
		public ActivationFunctions ActivationFunctionId;
		public int LayerIndex;
		public int NeuronCount;
		public int WeightCount;
		public int MapCount;
		public int MapWidth;
		public int MapHeight;
		public int MapSize;
		public int ReceptiveFieldWidth;
		public int ReceptiveFieldHeight;
		public int ReceptiveFieldSize;
		public int StrideX;
		public int StrideY;
		public int PadX;
		public int PadY;
		public int DropOutPercentage;
		public bool UseMapInfo;
		public bool IsFullyMapped;
		public bool HasWeights;
		public bool LockedWeights;
		public bool UseDropOut;
		public bool UseWeightPartitioner;
		public bool UseNeuronPartitioner;
		public double SubsamplingScalingFactor;
		public double DropOutFactor;
		public double InitWeight;
		public double InitBias;
		public Layer PreviousLayer;
		public Layer NextLayer;
		public Func<double, double> ActivationFunction;
		public Func<double,double> DerivativeActivationFunction;
		public Action CalculateAction;
		public Action BackpropagateAction;
		public Action EraseGradientWeights;
		public Action<int> UpdateWeights;
		public Func<double[], double[]> BackpropagateSecondDerivativesAction;
		public Mappings Mappings;
        public double[][] Gaussian2DKernel;
		public int[] NeuronActive;
		public Neuron[] Neurons;
		public Weight[] Weights;
        public readonly Connection[][] Connections;
		
		public Layer(NeuralNetwork network, LayerTypes layerType, int mapCount, int mapWidth, int mapHeight) : this(network, network.Layers.Length, layerType, ActivationFunctions.None, mapCount * mapWidth * mapHeight, true, mapCount, mapWidth, mapHeight, true, 0, 0, 1, 1, 0, 0, ((network.Layers.Length == 0) ? (null) : (network.Layers[network.Layers.Length - 1])), null, false, false, 50, 0.01D, 0) { }
		public Layer(NeuralNetwork network, LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight) : this(network, network.Layers.Length, layerType, activationFunction, mapCount * mapWidth * mapHeight, true, mapCount, mapWidth, mapHeight, true, 0, 0, 1, 1, 0, 0, ((network.Layers.Length == 0) ? (null) : (network.Layers[network.Layers.Length - 1])), null, false, false, 50, 0.01D, 0) { }
		public Layer(NeuralNetwork network, LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight, int receptiveFieldWidth, int receptiveFieldHeight, int strideX, int strideY, int padX, int padY, Mappings mappings, bool useDropOut = false, int dropOutPercentage = 50, double initWeight = 0.01D, double initBias = 0D) : this(network, network.Layers.Length, layerType, activationFunction, mapCount * mapWidth * mapHeight, true, mapCount, mapWidth, mapHeight, mappings == null ? true : false, receptiveFieldWidth, receptiveFieldHeight, strideX, strideY, padX, padY, network.Layers[network.Layers.Length - 1], mappings, false, useDropOut, dropOutPercentage, initWeight, initBias) { }
		public Layer(NeuralNetwork network, LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight, int receptiveFieldWidth, int receptiveFieldHeight, int strideX, int strideY, int padX, int padY, bool useDropOut = false, int dropOutPercentage = 50, double initWeight = 0.01D, double initBias = 0D) : this(network, network.Layers.Length, layerType, activationFunction, mapCount * mapWidth * mapHeight, true, mapCount, mapWidth, mapHeight, true, receptiveFieldWidth, receptiveFieldHeight, strideX, strideY, padX, padY, network.Layers[network.Layers.Length - 1], null, false, useDropOut, dropOutPercentage, initWeight, initBias) { }
		public Layer(NeuralNetwork network, LayerTypes layerType, ActivationFunctions activationFunction, int mapCount, int mapWidth, int mapHeight, int receptiveFieldWidth, int receptiveFieldHeight, int strideX = 1, int strideY = 1) : this(network, network.Layers.Length, layerType, activationFunction, mapCount * mapWidth * mapHeight, true, mapCount, mapWidth, mapHeight, true, receptiveFieldWidth, receptiveFieldHeight, strideX, strideY, 0, 0, network.Layers[network.Layers.Length - 1], null, false, false, 50, 0.01D, 0) { }
		public Layer(NeuralNetwork network, LayerTypes layerType, ActivationFunctions activationFunction, int neuronCount, bool useDropOut = false, int dropOutPercentage = 50, double initWeight = 0.01D, double initBias = 0D) : this(network, network.Layers.Length, layerType, activationFunction, neuronCount, false, 1, 1, 1, true, 0, 0, 0, 0, 0, 0, ((network.Layers.Length == 0) ? (null) : (network.Layers[network.Layers.Length - 1])), null, false, useDropOut, dropOutPercentage, initWeight, initBias) { }
		public Layer(NeuralNetwork network, int layerIndex, LayerTypes layerType, ActivationFunctions activationFunction, int neuronCount, bool useMapInfo, int mapCount, int mapWidth, int mapHeight, bool isFullyMapped, int receptiveFieldWidth, int receptiveFieldHeight, int strideX, int strideY, int padX, int padY, Layer previousLayer, Mappings mappings, bool lockedWeights, bool useDropOut, int dropOutPercentage = 50, double initWeight = 0.01D, double initBias = 0D)
		{
			Network = network;
			LayerIndex = layerIndex;
			LayerType = layerType;
			ActivationFunctionId = activationFunction;
			NeuronCount = neuronCount;
			UseMapInfo = useMapInfo;
			MapCount = mapCount;
			MapWidth = mapWidth;
			MapHeight = mapHeight;
			MapSize = MapWidth * MapHeight;
			IsFullyMapped = isFullyMapped;
			ReceptiveFieldWidth = receptiveFieldWidth;
			ReceptiveFieldHeight = receptiveFieldHeight;
			ReceptiveFieldSize = ReceptiveFieldWidth * ReceptiveFieldHeight;
			StrideX = strideX;
			StrideY = strideY;
			PadX = padX;
			PadY = padY;
			Mappings = mappings;
			PreviousLayer = previousLayer;
			UseDropOut = useDropOut;
			DropOutPercentage = dropOutPercentage;
			DropOutFactor = 1D / (1D - ((double)DropOutPercentage / 100D));
			LockedWeights = lockedWeights;
			InitWeight = initWeight;
			InitBias = initBias;
			
			if (UseDropOut)
				Network.DropOutUsed = true;

			NeuronActive = new int[NeuronCount];
			Neurons = new Neuron[NeuronCount];
			Connections = new Connection[NeuronCount][];
			for (int i = 0; i < NeuronCount; i++)
			{
				NeuronActive[i] = 1;
				Neurons[i].Output = 0D;
				Neurons[i].D1ErrX = 0D;
				Connections[i] = new Connection[0];
			}
			UseNeuronPartitioner = NeuronCount > network.ParallelTreshold;
			
			int totalMappings = 0;
			int maskWidth = 0;
			int maskHeight = 0;
			int maskSize = 0;
			int cMid = ReceptiveFieldWidth / 2;
			int rMid = ReceptiveFieldHeight / 2;
			int[] kernelTemplate;
			int[] maskMatrix;

			switch (LayerType)
			{
				case LayerTypes.Input:
					ActivationFunctionId = ActivationFunctions.None;
					HasWeights = false;
					ActivationFunction = null;
					UseWeightPartitioner = false;
					WeightCount = 0;
					Weights = null;
					CalculateAction = null;
					BackpropagateAction = null;
					BackpropagateSecondDerivativesAction = null;
					break;

				case LayerTypes.Local:
					if (IsFullyMapped)
						totalMappings = PreviousLayer.MapCount * MapCount;
					else
					{
						if (Mappings != null)
						{
							if (Mappings.Mapping.Count() == PreviousLayer.MapCount * MapCount)
								totalMappings = Mappings.Mapping.Count(p => p == true);
							else
								throw new ArgumentException("Invalid mappings definition");
						}
						else
							throw new ArgumentException("Empty mappings definition");
					}

					HasWeights = true;
                    WeightCount = totalMappings * MapSize * (ReceptiveFieldSize + 1);
					Weights = new Weight[WeightCount];
					UseWeightPartitioner = WeightCount > network.ParallelTreshold;

					CalculateAction = CalculateCCF; //CalculateLocalConnected;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateCCFParallel;
					else
						BackpropagateAction = BackpropagateCCFSerial;

					if (UseDropOut)
					{
						CalculateAction = CalculateCCFDropOut; // CalculateLocalConnectedDropOut;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateCCFDropOutParallel;
						else
							BackpropagateAction = BackpropagateCCFDropOutSerial;
					}

					if (UseWeightPartitioner)
					{
						EraseGradientWeights = EraseGradientsWeightsParallel;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
					}
					else
					{
						EraseGradientWeights = EraseGradientsWeightsSerial;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
					}

					ChangeTrainingStrategy();

					maskWidth = PreviousLayer.MapWidth + (2 * PadX);
					maskHeight = PreviousLayer.MapHeight + (2 * PadY);
					maskSize = maskWidth * maskHeight;

					kernelTemplate = new int[ReceptiveFieldSize];
					for (int row=0; row < ReceptiveFieldHeight; row++)
						for (int column = 0; column < ReceptiveFieldWidth; column++)
							kernelTemplate[column + (row * ReceptiveFieldWidth)] = column + (row * maskWidth);
					
					maskMatrix = new int[maskSize * PreviousLayer.MapCount];
					for (int i = 0; i < maskSize * PreviousLayer.MapCount; i++)
						maskMatrix[i] = -1;
					Parallel.For(0, PreviousLayer.MapCount, map =>
					{
						for (int y = PadY; y < PreviousLayer.MapHeight + PadY; y++)
							for (int x = PadX; x < PreviousLayer.MapWidth + PadX; x++)
								maskMatrix[x + (y * maskHeight) + (map * maskSize)] = (x - PadX) + ((y - PadY) * PreviousLayer.MapWidth) + (map * PreviousLayer.MapSize);
					});

                    if (!IsFullyMapped)
                    {
                        int mapping = 0;
                        int[] mappingCount = new int[MapCount * PreviousLayer.MapCount];
                        for (int curMap = 0; curMap < MapCount; curMap++)
                            for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
                            {
                                mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] = mapping;
                                if (Mappings.IsMapped(curMap, prevMap, MapCount))
                                    mapping++;
                            }

                        Parallel.For(0, MapCount, curMap =>
						{
							for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
							{
								int positionPrevMap = prevMap * maskSize;

								if (Mappings.IsMapped(curMap, prevMap, MapCount))
								{
									int iNumWeight = mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] * (ReceptiveFieldSize + 1) * MapSize;

									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
											
											AddBias(ref Connections[position], iNumWeight++);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
								}
							}
						});
					}
					else // Fully mapped
					{
						if (totalMappings > MapCount)
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * maskSize;
									int mapping = prevMap + (curMap * PreviousLayer.MapCount);
									int iNumWeight = mapping * (ReceptiveFieldSize + 1) * MapSize;
									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
									  
											AddBias(ref Connections[position], iNumWeight++);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
								}
							});
						}
						else // PreviousLayer has only one map         // 36*36 phantom input , padXY=2, filterSize=5, results in 32x32 conv layer
						{
							Parallel.For(0, MapCount, curMap =>
							{
								int iNumWeight = curMap * (ReceptiveFieldSize + 1) * MapSize;
								for (int y = 0; y < MapHeight; y++)
									for (int x = 0; x < MapWidth; x++)
									{
										int position = x + (y * MapWidth) + (curMap * MapSize);
										
										AddBias(ref Connections[position], iNumWeight++);

										int pIndex;
										for (int row = 0; row < ReceptiveFieldHeight; row++)
											for (int column = 0; column < ReceptiveFieldWidth; column++)
											{
												pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)];
												if (maskMatrix[pIndex] != -1)
													AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
											}
									}
							});
						}
					}
					break;

				case LayerTypes.Convolutional:
					if (IsFullyMapped)
						totalMappings = PreviousLayer.MapCount * MapCount;
					else
					{
						if (Mappings != null)
						{
							if (Mappings.Mapping.Count() == PreviousLayer.MapCount * MapCount)
								totalMappings = Mappings.Mapping.Count(p => p == true);
							else
								throw new ArgumentException("Invalid mappings definition");
						}
						else
							throw new ArgumentException("Empty mappings definition");
					}

					HasWeights = true;
					WeightCount = (totalMappings * ReceptiveFieldSize) + MapCount;
					Weights = new Weight[WeightCount];
					UseWeightPartitioner = WeightCount > network.ParallelTreshold;
				   
					CalculateAction = CalculateCCF;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateCCFParallel;
					else
						BackpropagateAction = BackpropagateCCFSerial;

					if (UseDropOut)
					{
						CalculateAction = CalculateCCFDropOut;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateCCFDropOutParallel;
						else
							BackpropagateAction = BackpropagateCCFDropOutSerial;
					}

					if (UseWeightPartitioner)
					{
						EraseGradientWeights = EraseGradientsWeightsParallel;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
					}
					else
					{
						EraseGradientWeights = EraseGradientsWeightsSerial;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
					}

					ChangeTrainingStrategy();
					
					maskWidth = PreviousLayer.MapWidth + (2 * PadX);
					maskHeight = PreviousLayer.MapHeight + (2 * PadY);
					maskSize = maskWidth * maskHeight;

					kernelTemplate = new int[ReceptiveFieldSize];
					for (int row=0; row < ReceptiveFieldHeight; row++)
						for (int column = 0; column < ReceptiveFieldWidth; column++)
							kernelTemplate[column + (row * ReceptiveFieldWidth)] = column + (row * maskWidth);
					
					maskMatrix = new int[maskSize * PreviousLayer.MapCount];
					for (int i = 0; i < maskSize * PreviousLayer.MapCount; i++)
						maskMatrix[i] = -1;
					Parallel.For(0, PreviousLayer.MapCount, map =>
					{
						for (int y = PadY; y < PreviousLayer.MapHeight + PadY; y++)
							for (int x = PadX; x < PreviousLayer.MapWidth + PadX; x++)
								maskMatrix[x + (y * maskWidth) + (map * maskSize)] = (x - PadX) + ((y - PadY) * PreviousLayer.MapWidth) + (map * PreviousLayer.MapSize);
					});

					if (!IsFullyMapped)
					{
						int mapping = 0;
						int[] mappingCount = new int[MapCount * PreviousLayer.MapCount];
						for (int curMap = 0; curMap < MapCount; curMap++)
							for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
							{
								mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] = mapping;
								if (Mappings.IsMapped(curMap, prevMap, MapCount))
									mapping++;
							}

						Parallel.For(0, MapCount, curMap =>
						{
							for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
							{
								int positionPrevMap = prevMap * maskSize;

								if (Mappings.IsMapped(curMap, prevMap, MapCount))
									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
											int iNumWeight = (mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] * ReceptiveFieldSize) + MapCount;

											AddBias(ref Connections[position], curMap);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
							}
						});
					}
					else // Fully mapped
					{
						if (totalMappings > MapCount)
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * maskSize;
									int mapping = prevMap + (curMap * PreviousLayer.MapCount);
									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
											int iNumWeight = (mapping * ReceptiveFieldSize) + MapCount;

											AddBias(ref Connections[position], curMap);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
								}
							});
						}
						else // PreviousLayer has only one map         // 36*36 phantom input , padXY=2, filterSize=5, results in 32x32 conv layer
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int y = 0; y < MapHeight; y++)
									for (int x = 0; x < MapWidth; x++)
									{
										int position = x + (y * MapWidth) + (curMap * MapSize);
										int iNumWeight = MapCount + (curMap * ReceptiveFieldSize);

										AddBias(ref Connections[position], curMap);

										int pIndex;
										for (int row = 0; row < ReceptiveFieldHeight; row++)
											for (int column = 0; column < ReceptiveFieldWidth; column++)
											{
												pIndex = x + (y * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)];
												if (maskMatrix[pIndex] != -1)
													AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
											}
									}
							});
						}
					}
					break;

				case LayerTypes.ConvolutionalSubsampling:  // Simard's implementation
					if (IsFullyMapped)
						totalMappings = PreviousLayer.MapCount * MapCount;
					else
					{
						if (Mappings != null)
						{
							if (Mappings.Mapping.Count() == PreviousLayer.MapCount * MapCount)
								totalMappings = Mappings.Mapping.Count(p => p == true);
							else
								throw new ArgumentException("Invalid mappings definition");
						}
						else
							throw new ArgumentException("Empty mappings definition");
					}

					HasWeights = true;
					WeightCount = (totalMappings * ReceptiveFieldSize) + MapCount;
					Weights = new Weight[WeightCount];
					UseWeightPartitioner = WeightCount > network.ParallelTreshold;
				   
					CalculateAction = CalculateCCF;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateCCFParallel;
					else
						BackpropagateAction = BackpropagateCCFSerial;

					if (UseDropOut)
					{
						CalculateAction = CalculateCCFDropOut;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateCCFDropOutParallel;
						else
							BackpropagateAction = BackpropagateCCFDropOutSerial;
					}

					if (UseWeightPartitioner)
					{
						EraseGradientWeights = EraseGradientsWeightsParallel;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
					}
					else
					{
						EraseGradientWeights = EraseGradientsWeightsSerial;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
					}

					ChangeTrainingStrategy();            

					maskWidth = PreviousLayer.MapWidth + (2 * PadX);
					maskHeight = PreviousLayer.MapHeight + (2 * PadY);
					maskSize = maskWidth * maskHeight;

					kernelTemplate = new int[ReceptiveFieldSize];
					for (int row=0; row < ReceptiveFieldHeight; row++)
						for (int column = 0; column < ReceptiveFieldWidth; column++)
							kernelTemplate[column + (row * ReceptiveFieldWidth)] = column + (row * maskWidth);
					
					maskMatrix = new int[maskSize * PreviousLayer.MapCount];
					for (int i = 0; i < maskSize * PreviousLayer.MapCount; i++)
						maskMatrix[i] = -1;
					for (int map = 0; map < PreviousLayer.MapCount; map++)
						for (int y = PadY; y < PreviousLayer.MapHeight + PadY; y++)
							for (int x = PadX; x < PreviousLayer.MapWidth + PadX; x++)
								maskMatrix[x + (y * maskHeight) + (map * maskSize)] = (x - PadX) + ((y - PadY) * PreviousLayer.MapWidth) + (map * PreviousLayer.MapSize);
						
					if (!IsFullyMapped) // not fully mapped
					{
						int mapping = 0;
						int[] mappingCount = new int[MapCount * PreviousLayer.MapCount];
						for (int curMap = 0; curMap < MapCount; curMap++)
							for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
							{
								mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] = mapping;
								if (Mappings.IsMapped(curMap, prevMap, MapCount))
									mapping++;
							}

						Parallel.For(0, MapCount, curMap =>
						{
							for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
							{
								int positionPrevMap = prevMap * maskSize;

								if (Mappings.IsMapped(curMap, prevMap, MapCount))
									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
											int iNumWeight = (mappingCount[prevMap + (curMap * PreviousLayer.MapCount)] * ReceptiveFieldSize) + MapCount;

											AddBias(ref Connections[position], curMap);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = (x * 2) + (y * 2 * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
							}
						});
					}
					else // Fully mapped
					{
						if (totalMappings > MapCount)
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * maskSize;
									int mapping = prevMap + (curMap * PreviousLayer.MapCount);
									for (int y = 0; y < MapHeight; y++)
										for (int x = 0; x < MapWidth; x++)
										{
											int position = x + (y * MapWidth) + (curMap * MapSize);
											int iNumWeight = (mapping * ReceptiveFieldSize) + MapCount;

											AddBias(ref Connections[position], curMap);

											int pIndex;
											for (int row = 0; row < ReceptiveFieldHeight; row++)
												for (int column = 0; column < ReceptiveFieldWidth; column++)
												{
													pIndex = (x * 2) + (y * 2 * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)] + positionPrevMap;
													if (maskMatrix[pIndex] != -1)
														AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
												}
										}
								}
							});
						}
						else // PreviousLayer has only one map         // 36*36 phantom input , padXY=2, filterSize=5, results in 32x32 conv layer
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int y = 0; y < MapHeight; y++)
									for (int x = 0; x < MapWidth; x++)
									{
										int position = x + (y * MapWidth) + (curMap * MapSize);
										int iNumWeight = MapCount + (curMap * ReceptiveFieldSize);

										AddBias(ref Connections[position], curMap);
										
										int pIndex;
										for (int row = 0; row < ReceptiveFieldHeight; row++)
											for (int column = 0; column < ReceptiveFieldWidth; column++)
											{
												pIndex = (x * 2) + (y * 2 * maskWidth) + kernelTemplate[column + (row * ReceptiveFieldWidth)];
												if (maskMatrix[pIndex] != -1)
													AddConnection(ref Connections[position], maskMatrix[pIndex], iNumWeight++);
											}
									}
							});
						}
					}
					break;

				case LayerTypes.AvgPooling:
				case LayerTypes.MaxPooling:
					HasWeights = true;
					WeightCount = MapCount * 2;
					Weights = new Weight[WeightCount];
				   
					if (LayerType == LayerTypes.AvgPooling)
					{
						CalculateAction = CalculateAveragePooling;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateAveragePoolingParallel;
						else
							BackpropagateAction = BackpropagateAveragePoolingSerial;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesAveragePooling;
					}
					else
					{
						CalculateAction = CalculateMaxPooling;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateMaxPoolingParallel;
						else
							BackpropagateAction = BackpropagateMaxPoolingSerial;
						BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesMaxPooling;
					}

					if (UseWeightPartitioner)
						EraseGradientWeights = EraseGradientsWeightsParallel;
					else
						EraseGradientWeights = EraseGradientsWeightsSerial;

					ChangeTrainingStrategy();

					SubsamplingScalingFactor = 1.0D / (StrideX * StrideY);
				   
					if (PreviousLayer.MapCount > 1) //fully symmetrical mapped
					{
						if (ReceptiveFieldSize != (StrideX * StrideY))
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * PreviousLayer.MapSize;

									if (prevMap == curMap)
									{
										for (int y = 0; y < MapHeight; y++)
											for (int x = 0; x < MapWidth; x++)
											{
												int position = x + (y * MapWidth) + (curMap * MapSize);
												int iNumWeight = curMap * 2;
												AddBias(ref Connections[position], iNumWeight++);

												bool outOfBounds = false;
												for (int row = -rMid; row <= rMid; row++)
													for (int col = -cMid; col <= cMid; col++)
													{
														if (row + (y * StrideY) < 0)
															outOfBounds = true;
														if (row + (y * StrideY) >= PreviousLayer.MapHeight)
															outOfBounds = true;
														if (col + (x * StrideX) < 0)
															outOfBounds = true;
														if (col + (x * StrideX) >= PreviousLayer.MapWidth)
															outOfBounds = true;
														if (!outOfBounds)
															AddConnection(ref Connections[position], (col + (x * StrideX)) + ((row + (y * StrideY)) * PreviousLayer.MapWidth) + positionPrevMap, iNumWeight);
														else
															outOfBounds = false;
													}
											}
									}
								}
							});
						}
						else
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * PreviousLayer.MapWidth * PreviousLayer.MapHeight;

									if (prevMap == curMap)
									{
										for (int y = 0; y < MapHeight; y++)
											for (int x = 0; x < MapWidth; x++)
											{
												int position = x + (y * MapWidth) + (curMap * MapSize);
												int iNumWeight = curMap * 2;
												AddBias(ref Connections[position], iNumWeight++);

												for (int row = 0; row < ReceptiveFieldHeight; row++)
													for (int col = 0; col < ReceptiveFieldWidth; col++)
														AddConnection(ref Connections[position], (col + (x * StrideX)) + ((row + (y * StrideY)) * PreviousLayer.MapWidth) + positionPrevMap, iNumWeight);
											}
									}
								}
							});
						}
					}
					break;

				case LayerTypes.AvgPoolingWeightless:
				case LayerTypes.MaxPoolingWeightless:
				case LayerTypes.L2Pooling:
				case LayerTypes.StochasticPooling:
					HasWeights = false;
					WeightCount = 0;
					Weights = null;

					switch (LayerType)
					{
						case LayerTypes.AvgPoolingWeightless:
							CalculateAction = CalculateAveragePoolingWeightless;
							if (PreviousLayer.UseNeuronPartitioner)
								BackpropagateAction = BackpropagateAveragePoolingWeightlessParallel;
							else
								BackpropagateAction = BackpropagateAveragePoolingWeightlessSerial;
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesAveragePoolingWeightless;
							break;

						case LayerTypes.MaxPoolingWeightless:
							CalculateAction = CalculateMaxPoolingWeightless;
							if (PreviousLayer.UseNeuronPartitioner)
								BackpropagateAction = BackpropagateMaxPoolingWeightlessParallel;
							else
								BackpropagateAction = BackpropagateMaxPoolingWeightlessSerial;
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesMaxPoolingWeightless;
							break;

						case LayerTypes.L2Pooling:
							CalculateAction = CalculateL2Pooling;
							if (PreviousLayer.UseNeuronPartitioner)
								BackpropagateAction = BackpropagateL2PoolingParallel;
							else
								BackpropagateAction = BackpropagateL2PoolingSerial;
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesL2Pooling;
							break;

						case LayerTypes.StochasticPooling:
							CalculateAction = CalculateStochasticPooling;
							if (PreviousLayer.UseNeuronPartitioner)
								BackpropagateAction = BackpropagateStochasticPoolingParallel;
							else
								BackpropagateAction = BackpropagateStochasticPoolingSerial;
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesStochasticPooling;
							break;
					}

					EraseGradientWeights = NoErase;
					UpdateWeights = NoUpdate;

					SubsamplingScalingFactor = 1.0D / (StrideX * StrideY);
				   
					if (PreviousLayer.MapCount > 1) //fully symmetrical mapped
					{
						if (ReceptiveFieldSize != (StrideX * StrideY))
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * PreviousLayer.MapSize;

									if (prevMap == curMap)
										for (int y = 0; y < MapHeight; y++)
											for (int x = 0; x < MapWidth; x++)
											{
												int position = x + (y * MapWidth) + (curMap * MapSize);
												bool outOfBounds = false;
												for (int row = -rMid; row <= rMid; row++)
													for (int col = -cMid; col <= cMid; col++)
													{
														if (row + (y * StrideY) < 0)
															outOfBounds = true;
														if (row + (y * StrideY) >= PreviousLayer.MapHeight)
															outOfBounds = true;
														if (col + (x * StrideX) < 0)
															outOfBounds = true;
														if (col + (x * StrideX) >= PreviousLayer.MapWidth)
															outOfBounds = true;
														if (!outOfBounds)
															AddConnection(ref Connections[position], (col + (x * StrideX)) + ((row + (y * StrideY)) * PreviousLayer.MapWidth) + positionPrevMap, 0);
														else
															outOfBounds = false;
													}
											}
								}
							});
						}
						else
						{
							Parallel.For(0, MapCount, curMap =>
							{
								for (int prevMap = 0; prevMap < PreviousLayer.MapCount; prevMap++)
								{
									int positionPrevMap = prevMap * PreviousLayer.MapWidth * PreviousLayer.MapHeight;

									if (prevMap == curMap)
										for (int y = 0; y < MapHeight; y++)
											for (int x = 0; x < MapWidth; x++)
											{
												int position = x + (y * MapWidth) + (curMap * MapSize);
											   
												for (int row = 0; row < ReceptiveFieldHeight; row++)
													for (int col = 0; col < ReceptiveFieldWidth; col++)
														AddConnection(ref Connections[position], (col + (x * StrideX)) + ((row + (y * StrideY)) * PreviousLayer.MapWidth) + positionPrevMap, 0);
											}
								}
							});
						}
					}
					break;

				case LayerTypes.FullyConnected:
					HasWeights = true;
					WeightCount = (PreviousLayer.NeuronCount + 1) * NeuronCount;
					Weights = new Weight[WeightCount];
					UseWeightPartitioner = WeightCount > network.ParallelTreshold;
					
					if (ActivationFunctionId == ActivationFunctions.SoftMax)
						CalculateAction = CalculateSoftMax;
					else
						CalculateAction = CalculateFullyConnected;

					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateCCFParallel;
					else
						BackpropagateAction = BackpropagateCCFSerial;

					if (UseDropOut)
					{
						CalculateAction = CalculateCCFDropOut;
						if (PreviousLayer.UseNeuronPartitioner)
							BackpropagateAction = BackpropagateCCFDropOutParallel;
						else
							BackpropagateAction = BackpropagateCCFDropOutSerial;
					}

					if (UseWeightPartitioner)
					{
						EraseGradientWeights = EraseGradientsWeightsParallel;
                        if ((ActivationFunctionId == ActivationFunctions.SoftMax) && ((Network.TrainingStrategy == TrainingStrategy.SGDLevenbergMarquardtModA) || (Network.TrainingStrategy == TrainingStrategy.MiniBatchSGDLevenbergMarquardtModA)))
                            BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxParallel;
                        else
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
					}
					else
					{
						EraseGradientWeights = EraseGradientsWeightsSerial;
                        if ((ActivationFunctionId == ActivationFunctions.SoftMax) && ((Network.TrainingStrategy == TrainingStrategy.SGDLevenbergMarquardtModA) || (Network.TrainingStrategy == TrainingStrategy.MiniBatchSGDLevenbergMarquardtModA)))
                            BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxSerial;
                        else
							BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
					}

					if (UseWeightPartitioner)
						EraseGradientWeights = EraseGradientsWeightsParallel;
					else
						EraseGradientWeights = EraseGradientsWeightsSerial;

					ChangeTrainingStrategy();
					
					if (UseMapInfo)
					{
						int iNumWeight = 0;
						Parallel.For(0, MapCount, curMap =>
						{
							for (int yc = 0; yc < MapHeight; yc++)
								for (int xc = 0; xc < MapWidth; xc++)
								{
									int position = xc + (yc * MapWidth) + (curMap * MapSize);
									AddBias(ref Connections[position], iNumWeight++);

									for (int prevMaps = 0; prevMaps < PreviousLayer.MapCount; prevMaps++)
										for (int y = 0; y < PreviousLayer.MapHeight; y++)
											for (int x = 0; x < PreviousLayer.MapWidth; x++)
												AddConnection(ref Connections[position], (x + (y * PreviousLayer.MapWidth) + (prevMaps * PreviousLayer.MapSize)), iNumWeight++);

								}
						});
					}
					else
					{
						int iNumWeight = 0;
						for (int y=0; y < NeuronCount; y++)
						{
							AddBias(ref Connections[y], iNumWeight++);
							for (int x = 0; x < PreviousLayer.NeuronCount; x++)
								AddConnection(ref Connections[y], x, iNumWeight++);
						}
					}
					break;

				case LayerTypes.RBF:
					HasWeights = true;
					WeightCount = PreviousLayer.NeuronCount * NeuronCount; // no biasses
					Weights = new Weight[WeightCount];
					
					if (ActivationFunctionId == ActivationFunctions.SoftMax)
						CalculateAction = CalculateSoftMaxRBF;
					else
						CalculateAction = CalculateRBF;
					
					BackpropagateAction = BackpropagateRBF;
					BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesRBF;

					if (UseWeightPartitioner)
						EraseGradientWeights = EraseGradientsWeightsParallel;
					else
						EraseGradientWeights = EraseGradientsWeightsSerial;

					ChangeTrainingStrategy();

					if (UseMapInfo)
					{
						int iNumWeight = 0;
						for (int n = 0; n < NeuronCount; n++)
							for (int prevMaps = 0; prevMaps < PreviousLayer.MapCount; prevMaps++)
								for (int y = 0; y < PreviousLayer.MapHeight; y++)
									for (int x = 0; x < PreviousLayer.MapWidth; x++)
										AddConnection(ref Connections[n], (x + (y * PreviousLayer.MapWidth) + (prevMaps * PreviousLayer.MapSize)), iNumWeight++);
					}
					else
					{
						int iNumWeight = 0;
						for (int y = 0; y < NeuronCount; y++)
							for (int x = 0; x < PreviousLayer.NeuronCount; x++)
								AddConnection(ref Connections[y], x, iNumWeight++);
					}
					break;

				case LayerTypes.LocalResponseNormalization:   // same map
					if (ActivationFunctionId != ActivationFunctions.None)
						throw new Exception("LocalContrastNormaliztion layer cannot have an ActivationFunction, specify ActivationFunctions.None");
					if (MapHeight != PreviousLayer.MapHeight)
						throw new Exception("MapHeight must be equal to the MapHeight of the previous layer");
					if (MapWidth != PreviousLayer.MapWidth)
						throw new Exception("MapWidth must be equal to the MapWidth of the previous layer");

					HasWeights = false;
					WeightCount = 0;
					Weights = null;

					Parallel.For(0, MapCount, map =>
					{
						for (int y = 0; y < MapHeight; y++)
							for (int x = 0; x < MapWidth; x++)
							{
								int position = x + (y * MapWidth) + (map * MapSize);
								bool outOfBounds = false;
								for (int row = -rMid; row <= rMid; row++)
									for (int col = -cMid; col <= cMid; col++)
									{
										if (col + x < 0)
											outOfBounds = true;
										if (col + x >= MapWidth)
											outOfBounds = true;
										if (row + y < 0)
											outOfBounds = true;
										if (row + y >= MapHeight)
											outOfBounds = true;

										if (!outOfBounds)
											AddConnection(ref Connections[position], (col + x) + ((row + y) * MapWidth) + (map * MapSize), -1);
										else
											outOfBounds = false;
									}
							}
					});
   
					CalculateAction = CalculateLocalResponseNormalization;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateLocalResponseNormalizationParallel;
					else
						BackpropagateAction = BackpropagateLocalResponseNormalizationSerial;
					BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativeLocalResponseNormalization;
					
					EraseGradientWeights = NoErase;
					UpdateWeights = NoUpdate;
					break;

				case LayerTypes.LocalResponseNormalizationCM: // across maps
					if (ActivationFunctionId != ActivationFunctions.None)
						throw new Exception("LocalContrastNormaliztionCM layer cannot have an ActivationFunction, specify ActivationFunctions.None");
					if (MapHeight != PreviousLayer.MapHeight)
						throw new Exception("MapHeight must be equal to the MapHeight of the previous layer");
					if (MapWidth != PreviousLayer.MapWidth)
						throw new Exception("MapWidth must be equal to the MapWidth of the previous layer");

					HasWeights = false;
					WeightCount = 0;
					Weights = null;

					int size = 9;
					int a, b;
					for (int map = 0; map < MapCount; map++)
						for (int y = 0; y < MapHeight; y++)
							for (int x = 0; x < MapWidth; x++)
							{
								a = Math.Max(0, map - (size / 2));   // from map
								b = Math.Min(MapCount, map - (size / 2) + size); // to map
								int position = x + (y * MapWidth) + (map * MapSize);
								for (int f = a; f < b; f++)
									AddConnection(ref Connections[position], x + (y * MapWidth) + (f * MapSize), -1);
							}
					
					CalculateAction = CalculateLocalResponseNormalizationCM;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagationLocalResponseNormalizationCMParallel;
					else
						BackpropagateAction = BackpropagationLocalResponseNormalizationCMSerial;
					BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativeLocalResponseNormalizationCM;

					EraseGradientWeights = NoErase;
					UpdateWeights = NoUpdate;
					break;

				case LayerTypes.LocalContrastNormalization:
					if (ActivationFunctionId != ActivationFunctions.None)
						throw new Exception("LocalContrastNormaliztion layer cannot have an ActivationFunction, specify ActivationFunctions.None");
					if (MapHeight != PreviousLayer.MapHeight)
						throw new Exception("MapHeight must be equal to the MapHeight of the previous layer");
					if (MapWidth != PreviousLayer.MapWidth)
						throw new Exception("MapWidth must be equal to the MapWidth of the previous layer");
                    Gaussian2DKernel = MathUtil.CreateGaussian2DKernel(ReceptiveFieldWidth, ReceptiveFieldHeight);
					HasWeights = false;
					WeightCount = 0;
					Weights = null;

					Parallel.For(0, MapCount, map =>
					{
						for (int y = 0; y < MapHeight; y++)
							for (int x = 0; x < MapWidth; x++)
							{
								int position = x + (y * MapWidth) + (map * MapSize);
								bool outOfBounds = false;
								for (int row = -rMid; row <= rMid; row++)
									for (int col = -cMid; col <= cMid; col++)
									{
										if (col + x < 0)
											outOfBounds = true;
										if (col + x >= MapWidth)
											outOfBounds = true;
										if (row + y < 0)
											outOfBounds = true;
										if (row + y >= MapHeight)
											outOfBounds = true;

										if (!outOfBounds)
											AddConnection(ref Connections[position], (col + x) + ((row + y) * MapWidth) + (map * MapSize), -1);
										else
											outOfBounds = false;
									}
							}
					});

					CalculateAction = CalculateLocalContrastNormalization;
					if (PreviousLayer.UseNeuronPartitioner)
						BackpropagateAction = BackpropagateLocalContrastNormalizationParallel;
					else
						BackpropagateAction = BackpropagateLocalContrastNormalizationSerial;
					BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativeLocalContrastNormalization;
					
					EraseGradientWeights = NoErase;
					UpdateWeights = NoUpdate;
					break;
			};

			switch (ActivationFunctionId)
			{
				case ActivationFunctions.Abs:
					ActivationFunction = Abs;
					DerivativeActivationFunction = DAbs;
					break;

				case ActivationFunctions.AbsTanh:
					ActivationFunction = AbsTanh;
					DerivativeActivationFunction = DAbsTanh;
					break;

				case ActivationFunctions.BReLU1:
					ActivationFunction = BReLU1;
					DerivativeActivationFunction = DBReLU1;
					break;

                case ActivationFunctions.BReLU6:
                    ActivationFunction = BReLU6;
                    DerivativeActivationFunction = DBReLU6;
                    break;

				case ActivationFunctions.Gaussian:
					ActivationFunction = Gaussian;
					DerivativeActivationFunction = DGaussian;
					break;

				case ActivationFunctions.Ident:
					ActivationFunction = Ident;
					DerivativeActivationFunction = DIdent;
					break;

				case ActivationFunctions.Logistic:
					ActivationFunction = Logistic;
					DerivativeActivationFunction = DLogistic;
					break;

				case ActivationFunctions.None:
					ActivationFunction = null;
					DerivativeActivationFunction = null;
					break;

				case ActivationFunctions.Tanh:
					ActivationFunction = Tanh;
					DerivativeActivationFunction = DTanh;
					break;

				case ActivationFunctions.Ramp:
					ActivationFunction = Ramp;
					DerivativeActivationFunction = DRamp;
					break; 

				case ActivationFunctions.ReLU:
					ActivationFunction = ReLU;
					DerivativeActivationFunction = DReLU;
					break;

				case ActivationFunctions.SoftMax:
					ActivationFunction = Ident;  
					DerivativeActivationFunction = DSoftMax;
					break;

				case ActivationFunctions.SoftReLU:
					ActivationFunction = SoftReLU;
					DerivativeActivationFunction = DSoftReLU;
					break;

				case ActivationFunctions.SoftSign:
					ActivationFunction = SoftSign;
					DerivativeActivationFunction = DSoftSign;
					break;

				case ActivationFunctions.Square:
					ActivationFunction = Square;
					DerivativeActivationFunction = DSquare;
					break;

				case ActivationFunctions.SquareRoot:
					ActivationFunction = SquareRoot;
					DerivativeActivationFunction = DSquareRoot;
					break;

				case ActivationFunctions.STanh:
					ActivationFunction = STanh;
					DerivativeActivationFunction = DSTanh;
					break;
			};

			int totalConnections = 0;
			for (int i = 0; i < NeuronCount; i++)
				totalConnections += Connections[i].Length;

			Name += "\r\nLayer: " + LayerIndex.ToString();
			Name += "\r\nLayer Type: " + LayerType.ToString() +
					((LayerType == LayerTypes.Input) ? ("") : ("\r\nActivation Function: " + ActivationFunctionId.ToString())) +
					((LayerType == LayerTypes.Local || LayerType == LayerTypes.Convolutional || LayerType == LayerTypes.AvgPooling || LayerType == LayerTypes.MaxPooling || LayerType == LayerTypes.AvgPoolingWeightless || LayerType == LayerTypes.MaxPoolingWeightless || LayerType == LayerTypes.StochasticPooling || LayerType == LayerTypes.L2Pooling) ? ("\r\nReceptive Field: " + ReceptiveFieldWidth.ToString() + "x" + ReceptiveFieldHeight.ToString() + "\r\nStride: " + StrideX.ToString() + "x" + StrideY.ToString() + "\r\nPadding: " + PadX.ToString() + "x" + PadY.ToString()) : "") +
					((UseMapInfo) ? ("\r\nMaps: " + MapCount.ToString() + "x(" + MapWidth.ToString() + "x" + MapHeight.ToString() + ")") : ("")) +
					"\r\nNeurons: " + NeuronCount.ToString() +
					((HasWeights) ? ("\r\nWeights: " + Weights.Count().ToString()) : ("")) +
					((LayerType == LayerTypes.Input) ? ("") : ("\r\nConnections: " + totalConnections.ToString())) +
					((UseDropOut) ? ("\r\nDropOut: " + DropOutPercentage.ToString() + "%") : (""));
			
			if (PreviousLayer != null)
				PreviousLayer.NextLayer = this;
		}

		public void AddConnection(ref Connection[] connections, int neuronIndex, int weightIndex)
		{
			Array.Resize(ref connections, connections.Length + 1);
			connections[connections.Length - 1] = new Connection(neuronIndex, weightIndex);
		}

		public void AddBias(ref Connection[] connections, int weightIndex)
		{
			Array.Resize(ref connections, connections.Length + 1);
			connections[connections.Length - 1] = new Connection(int.MaxValue, weightIndex);
		}

		public void ChangeTrainingStrategy()
		{
            if (HasWeights)
            {
                UpdateWeights = null;

                switch (Network.TrainingStrategy)
                {
                    case TrainingStrategy.SGD:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeighsSGDSerial(batchSize); };
                        break;

                    case TrainingStrategy.SGDLevenbergMarquardt:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDLMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDLMSerial(batchSize); };

                        if ((Network.LossFunction == LossFunctions.CrossEntropy) && (ActivationFunctionId == ActivationFunctions.SoftMax))
                        {
                            if (UseWeightPartitioner)
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
                            else
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
                        }
                        break;

                    case TrainingStrategy.SGDLevenbergMarquardtModA:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDLMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDLMSerial(batchSize); };

                        if ((Network.LossFunction == LossFunctions.CrossEntropy) && (ActivationFunctionId == ActivationFunctions.SoftMax))
                        {
                            if (UseWeightPartitioner)
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxParallel;
                            else
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxSerial;
                        }
                        break;

                    case TrainingStrategy.SGDM:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsSGDMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeighsSGDMSerial(batchSize); };
                        break;

                    case TrainingStrategy.MiniBatchSGD:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDSerial(batchSize); };
                        break;

                    case TrainingStrategy.MiniBatchSGDLevenbergMarquardt:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDLMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDLMSerial(batchSize); };

                        if ((Network.LossFunction == LossFunctions.CrossEntropy) && (ActivationFunctionId == ActivationFunctions.SoftMax))
                        {
                            if (UseWeightPartitioner)
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFParallel;
                            else
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSerial;
                        }
                        break;

                    case TrainingStrategy.MiniBatchSGDLevenbergMarquardtModA:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDLMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDLMSerial(batchSize); };

                        if ((Network.LossFunction == LossFunctions.CrossEntropy) && (ActivationFunctionId == ActivationFunctions.SoftMax))
                        {
                            if (UseWeightPartitioner)
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxParallel;
                            else
                                BackpropagateSecondDerivativesAction = BackpropagateSecondDerivativesCCFSoftMaxSerial;
                        }
                        break;

                    case TrainingStrategy.MiniBatchSGDM:
                        if (UseWeightPartitioner)
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDMParallel(batchSize); };
                        else
                            UpdateWeights = delegate(int batchSize) { UpdateWeightsMiniBatchSGDMSerial(batchSize); };
                        break;
                }
            }
            else
            {
                UpdateWeights = null;
                UpdateWeights = delegate(int batchSize) { NoUpdate(batchSize); };
            }
		}
			
		public override string ToString()
		{
			return Name;
		}

		private static double Abs(double value)
		{
			return Math.Abs(value);
		}

		private static double DAbs(double value)
		{
			return value > 0 ? 1.0D : -1D;
		}

		private static double AbsTanh(double value)
		{
			return Math.Abs(MathUtil.Tanh(value));
		}

		private static double DAbsTanh(double value)
		{
			double x = MathUtil.Tanh(value);
			return ((1D - MathUtil.Pow2(x)) * x) / Math.Abs(x);
		}

		private static double BReLU1(double value)
		{
			return value < 0D ? 0D : value > 1 ? 1 : value;
		}

		private static double DBReLU1(double value)
		{
			return value < 0D ? 0 : value > 1 ? 0 : 1;
		}

        private static double BReLU6(double value)
        {
            return value < 0D ? 0D : value > 6 ? 6 : value;
        }

        private static double DBReLU6(double value)
        {
            return value < 0D ? 0 : value > 6 ? 0 : 1;
        }

		private static double Ramp(double value)
		{
			return Math.Min(-1.0D, Math.Max(1.0D, value));
		}

		private static double DRamp(double value)
		{
			if ((-1D < value) && (value < 1D))
				return 1D;
			else
				return 0D;
		}

		private static double Ident(double value)
		{
			return value;
		}

		private static double DIdent(double value)
		{
			return 1D;
		}

		private static double Tanh(double value)
		{
			return MathUtil.Tanh(value);
		}

		private static double DTanh(double value)
		{
			return 1D - (value * value);
		}

		private static double STanh(double value)
		{
			return MathUtil.SymmetricTanh(value);
		}

		private static double DSTanh(double value)
		{
			return MathUtil.DSymmetricTanh(value);
		}

		private static double Logistic(double value)
		{
			//return 1D / (1D + MathUtil.Exp(-value));
			return 0.5D * (1D + (MathUtil.Tanh(value * 0.5D))); // non overflowing version of logistic
		}

		private static double DLogistic(double value)
		{
			return (value * (1D - value));
			//return value - (value * value);
		}

		private static double DSoftMax(double value)
		{
			return (value * (1D - value));
			//return value - (value * value);
		}

		private static double Gaussian(double value)
		{
			return MathUtil.Exp(-0.5D*value);
		}

		private static double DGaussian(double value)
		{
			return -(MathUtil.Exp(-0.5D*value) * 0.5D);
		}

		private static double Square(double value)
		{
			return value * value;
		}

		private static double DSquare(double value)
		{
			return 2D * value;
		}

		private static double ReLU(double value)
		{
			//return value < 0D ? 0D : value > 6 ? 6 : value;

			return value < 0D ? 0D : value;
		}

		private static double DReLU(double value)
		{
			//return (value > 0D) && (value <= 6) ? 1D : 0D;

			return value > 0D ? 1D : 0D;
		}

		private static double SoftReLU(double value)   // aka SoftPlus
		{
			if (value > 4D)
				return value;
			else
				return MathUtil.Log(1D + MathUtil.Exp(value));
		}

		private static double DSoftReLU(double value)
		{
			//return Logistic(value);
			if (value > 4D)
				return value;
			else
				return 0.5D * (1D + (MathUtil.Tanh(value * 0.5D))); // non overflowing version of logistic
		}

		private static double SoftSign(double value)
		{
			return value / (1D + Math.Abs(value));
		}

		private static double DSoftSign(double value)
		{
			double temp = Math.Abs(value);
			return temp / (MathUtil.Pow3(temp) + temp + (2D * MathUtil.Pow2(temp)));
		}

		public static double SquareRoot(double value)
		{
			return Math.Sqrt(value);
		}

		public static double DSquareRoot(double value)
		{
			if (value != 0D)
				return 1D / (2D * Math.Sqrt(value));

			return 0D;
		}

		public void InitializeWeights()
		{
			checked
			{
				switch (LayerType)
				{
					case LayerTypes.RBF:
						int index = 0;
						for (int i = 0; i < NeuronCount; i++)
						{
							byte[] weightImage = new byte[12];
							weightImage = Network.RbfWeights[index++].ToArray();
							double[] realWeights = new double[7 * 12];
							int row = 0;
							for (int y = 0; y < 12; y++)
							{
								row = (int)weightImage[y];
								realWeights[0 + (7 * y)] = (((128 & ~row) / 128) * 2) - 1;
								realWeights[1 + (7 * y)] = (((64 & ~row) / 64) * 2) - 1;
								realWeights[2 + (7 * y)] = (((32 & ~row) / 32) * 2) - 1;
								realWeights[3 + (7 * y)] = (((16 & ~row) / 16) * 2) - 1;
								realWeights[4 + (7 * y)] = (((8 & ~row) / 8) * 2) - 1;
								realWeights[5 + (7 * y)] = (((4 & ~row) / 4) * 2) - 1;
								realWeights[6 + (7 * y)] = (((2 & ~row) / 2) * 2) - 1;
							}

							foreach (Connection connection in Connections[i]) //84x
								Weights[connection.ToWeightIndex].Value = (realWeights[connection.ToNeuronIndex] == 1D) ? Network.TrainToValue : -Network.TrainToValue;
						}
						break;

					default:
						if (Weights != null && WeightCount > 0)
						{
							switch (ActivationFunctionId)
							{
								case ActivationFunctions.ReLU:
								case ActivationFunctions.BReLU1:
                                case ActivationFunctions.BReLU6:
								case ActivationFunctions.SoftReLU:
									{
										Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
										{
											double stdDev = 1D / (double)Math.Sqrt(Connections[i].Length);
											foreach (Connection connection in Connections[i])
												if (connection.ToNeuronIndex == int.MaxValue)
													Weights[connection.ToWeightIndex].Value = Network.RandomGenerator.NextDouble(stdDev, 0.00001D); //(LayerIndex == 1) ? 1.0D : Network.RandomGenerator.NextDouble(stdDev, 0D);
												else
													Weights[connection.ToWeightIndex].Value = Network.RandomGenerator.NextDouble(stdDev);
										});
									}
									break;

								default:
									{
										Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
										{
											double stdDev = 1D / (double)Math.Sqrt(Connections[i].Length);
											foreach (Connection connection in Connections[i])
												if (connection.ToNeuronIndex == int.MaxValue)
													Weights[connection.ToWeightIndex].Value = Network.RandomGenerator.NextDouble(stdDev);
												else
													Weights[connection.ToWeightIndex].Value = Network.RandomGenerator.NextDouble(stdDev);
										});
									}
									break;
							}
						}
						break;
				}
			}
		}
		
		public void CalculateSoftMax()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				//double bias = Weights[i * (PreviousLayer.NeuronCount+1)].Value;
				int idx = i * (PreviousLayer.NeuronCount + 1);
				double dSum = Weights[idx++].Value;
				for (int c = 0; c < PreviousLayer.NeuronCount; c++)
					dSum += Weights[c + idx].Value * PreviousLayer.Neurons[c].Output;
					//foreach (Connection connection in Connections[i])
					//{
					//    if (connection.ToNeuronIndex == int.MaxValue)
					//        dSum += Weights[connection.ToWeightIndex].Value;
					//    else
					//        dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
					//}
				Neurons[i].Output = dSum;
			});

			//Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			//{
			//    double dSum = 0D;
			//    foreach (Connection connection in Connections[i])
			//    {
			//        if (connection.ToNeuronIndex == int.MaxValue)
			//            dSum += Weights[connection.ToWeightIndex].Value;
			//        else
			//            dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
			//    }
			//    Neurons[i].Output = dSum;
			//});

			double hi = Neurons[0].Output;
			for (int i = 1; i < NeuronCount; i++)
				if (Neurons[i].Output > hi)
					hi = Neurons[i].Output;
			//double hi = Neurons.Max(neuron => neuron.Output);

			double total = 0D;
			for (int i = 0; i < NeuronCount; i++)
				total += MathUtil.Exp(Neurons[i].Output - hi);

			double logsumexp = hi + MathUtil.Log(total);
			for (int i = 0; i < NeuronCount; i++)
				Neurons[i].Output = MathUtil.Exp(Neurons[i].Output - logsumexp);        
		}

		public void CalculateSoftMaxRBF()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double dSum = 0D;
				foreach (Connection connection in Connections[i])
					dSum += MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
				Neurons[i].Output = dSum;
			});

			double hi = Neurons.Max(neuron => neuron.Output);

			//double total = 0D;
			//foreach (Neuron neuron in Neurons)
			//    total += MathUtil.Exp(-0.5D * (neuron.Output - hi));

			//double logsumexp = hi + MathUtil.Log(total);

			//for (int i = 0; i < NeuronCount; i++)
			//    Neurons[i].Output = MathUtil.Exp((-0.5*Neurons[i].Output) - logsumexp);        
		   
			double total = 0D;
			for (int i = 0; i < NeuronCount; i++)
			{
				Neurons[i].Output = MathUtil.Exp(-0.5D * (Neurons[i].Output - hi));
				total += Neurons[i].Output;
			}

			for (int i = 0; i < NeuronCount; i++)
				Neurons[i].Output /= total;
		}

		public void CalculateFullyConnected()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				int idx = i * (PreviousLayer.NeuronCount + 1);
				double dSum = Weights[idx++].Value;
				for (int c = 0; c < PreviousLayer.NeuronCount; c++)
					dSum += Weights[idx + c].Value * PreviousLayer.Neurons[c].Output;
				Neurons[i].Output = ActivationFunction(dSum);
			});
		}

		//public void CalculateLocalConnected()
		//{
		//    Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//    {
		//        int idx = i * (PreviousLayer.NeuronCount + 1);
		//        double dSum = Weights[idx++].Value;
		//        for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//            dSum += Weights[idx + c].Value * PreviousLayer.Neurons[c].Output;
		//        Neurons[i].Output = ActivationFunction(dSum);
		//    });
		//}

		//public void CalculateLocalConnectedDropOut()
		//{
		//    if (Network.OperationState == NetworkStates.Training)
		//    {
		//        Parallel.For(0, NeuronCount, Network.ParallelOption, i => NeuronActive[i] = (Network.RandomGenerator.NextPercentage() < DropOutPercentage) ? 1 : 0);

		//        Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//        {
		//            if (NeuronActive[i] == 1)
		//            {
		//                int idx = i * (PreviousLayer.NeuronCount + 1);
		//                double dSum = Weights[idx++].Value;
		//                for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//                    dSum += Weights[idx + c].Value * PreviousLayer.Neurons[c].Output;
		//                Neurons[i].Output = ActivationFunction(dSum);
		//            }
		//            else
		//                Neurons[i].Output = 0D;
		//        });
		//    }
		//    else
		//    {
		//        Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//        {
		//            int idx = i * (PreviousLayer.NeuronCount + 1);
		//            double dSum = Weights[idx++].Value;
		//            for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//                dSum += Weights[idx + c].Value * PreviousLayer.Neurons[c].Output;
		//            Neurons[i].Output = ActivationFunction(DropOutFactor * dSum);
		//        });
		//    }
		//}

		public void CalculateCCF()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double dSum = 0D;
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						dSum += Weights[connection.ToWeightIndex].Value;
					else
						dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
				}
				Neurons[i].Output = ActivationFunction(dSum);
			});
		}

		public void CalculateCCFDropOut()
		{
			if (Network.OperationState == NetworkStates.Training)
			{
				Parallel.For(0, NeuronCount, Network.ParallelOption, i => NeuronActive[i] = Network.RandomGenerator.NextPercentage() < DropOutPercentage ? 1 : 0);

				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					if (NeuronActive[i] == 1)
					{
						double dSum = 0D;
						foreach (Connection connection in Connections[i])
						{
							if (connection.ToNeuronIndex == int.MaxValue)
								dSum += Weights[connection.ToWeightIndex].Value;
							else
								dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
						}
						Neurons[i].Output = ActivationFunction(dSum);
					}
					else
						Neurons[i].Output = 0D;
				});
			}
			else
			{
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					double dSum = 0D;
					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							dSum += Weights[connection.ToWeightIndex].Value;
						else
							dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
					}
					Neurons[i].Output = ActivationFunction(DropOutFactor * dSum);
				});
			}
		}

		public void CalculateAveragePooling()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / (Connections[i].Length-1);
				double dSum = 0D;
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						dSum += Weights[connection.ToWeightIndex].Value;
					else
						dSum += Weights[connection.ToWeightIndex].Value * PreviousLayer.Neurons[connection.ToNeuronIndex].Output * sf;
				}
				Neurons[i].Output = ActivationFunction(dSum);
			});
		}

		public void CalculateMaxPooling()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double bias = 0D;
				double weight = 0D;
				double max = double.MinValue;
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						bias = Weights[connection.ToWeightIndex].Value;
					else
					{
						if (PreviousLayer.Neurons[connection.ToNeuronIndex].Output >= max)
						{
							weight = Weights[connection.ToWeightIndex].Value;
							max = PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
							NeuronActive[i] = connection.ToNeuronIndex;
						}
					}
				}
				Neurons[i].Output = ActivationFunction((max * weight) + bias);
			});
		}

		public void CalculateAveragePoolingWeightless()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double dSum = 0D;
				foreach (Connection connection in Connections[i])
					dSum += PreviousLayer.Neurons[connection.ToNeuronIndex].Output;

				Neurons[i].Output = ActivationFunction(dSum / Connections[i].Length);
			});

		}

		public void CalculateMaxPoolingWeightless()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double max = double.MinValue;
				foreach (Connection connection in Connections[i])
				{
					if (PreviousLayer.Neurons[connection.ToNeuronIndex].Output > max)
					{
						max = PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
						NeuronActive[i] = connection.ToNeuronIndex;
					}
				}
				Neurons[i].Output = ActivationFunction(max);
			});
		}

		public void CalculateStochasticPooling()
		{
			if ((Network.OperationState == NetworkStates.Training) || (Network.OperationState == NetworkStates.CalculatingHessian))
			{
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					double dSum = 0D;
					foreach (Connection connection in Connections[i])
						if (PreviousLayer.Neurons[connection.ToNeuronIndex].Output > 0D)
							dSum += PreviousLayer.Neurons[connection.ToNeuronIndex].Output;

					if (dSum > 0D)
					{
						int count = 0;
						List<Tuple<double, int>> elements = new List<Tuple<double, int>>();
						foreach (Connection connection in Connections[i])
							if (PreviousLayer.Neurons[connection.ToNeuronIndex].Output > 0D)
								elements.Add(new Tuple<double, int>(PreviousLayer.Neurons[connection.ToNeuronIndex].Output / dSum, count++));
							else
								elements.Add(new Tuple<double, int>(0D, count++));
						elements.Sort();
						
						double diceRoll = 1D - Network.RandomGenerator.NextDouble();
						double cumulative = 0D;
						for (int j = 0; j < count; j++)
						{
							cumulative += elements[j].Item1;
							if (diceRoll < cumulative)
							{
								NeuronActive[i] = Connections[i][elements[j].Item2].ToNeuronIndex;
								break;
							}
						}
					}
					else
						NeuronActive[i] = Connections[i][Network.RandomGenerator.Next(Connections[i].Length)].ToNeuronIndex;

					Neurons[i].Output = ActivationFunction(PreviousLayer.Neurons[NeuronActive[i]].Output);
				});
			}
			else
			{
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					double sum = 0D;
					foreach (Connection connection in Connections[i])
						sum += PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
				   
					double prob = 0D;
					if (sum > 0D)
					{
						foreach (Connection connection in Connections[i])
							prob += PreviousLayer.Neurons[connection.ToNeuronIndex].Output * (PreviousLayer.Neurons[connection.ToNeuronIndex].Output / sum);

						Neurons[i].Output = ActivationFunction(prob);
					}
					else
						Neurons[i].Output = ActivationFunction(0D);
				}); 
			}
		}

		public void CalculateL2Pooling()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double dSum = 0D;
				for (int c=0; c < Connections[i].Length; c++)
					dSum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				Neurons[i].Output = ActivationFunction(dSum / Connections[i].Length);
			});
		}

		public void CalculateLocalResponseNormalization()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0.0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				Neurons[i].Output = PreviousLayer.Neurons[i].Output / Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow);
			});
		}

		public void CalculateLocalResponseNormalizationCM()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				Neurons[i].Output = PreviousLayer.Neurons[i].Output / Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow);
			});
		}

		public void CalculateLocalContrastNormalization()
		{
            Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
                double mean = 0D;
                double[] a = new double[ReceptiveFieldSize];
                //double b = 0D;
                int x = 0;
                int y = 0;
                for (int c = 0; c < Connections[i].Length; c++)
                {
                    if (c % ReceptiveFieldWidth == 0)
                    {
                        x = 0;
                        y++;
                    }
                    mean += PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output * Gaussian2DKernel[x++][y];                   
                }
				mean /= Connections[i].Length;

				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output - mean);

				Neurons[i].Output = PreviousLayer.Neurons[i].Output / Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow);
			});
		}

		public void CalculateRBF()
		{
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double dSum = 0D;
				foreach (Connection connection in Connections[i])
					dSum += MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
				Neurons[i].Output = ActivationFunction(dSum);
			});
		}

		public void EraseGradientsWeightsSerial()
		{
			for (int i = 0; i < WeightCount; i++)
				Weights[i].D1Err = 0D;
		}

		public void EraseGradientsWeightsParallel()
		{
			Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].D1Err = 0D);
		}

		public void UpdateWeightsSGDLMSerial(int batchSize = 1)
		{
			for (int i = 0; i < WeightCount; i++)
				Weights[i].Value -= (Network.TrainingRate.Rate / (Weights[i].DiagonalHessian + Network.dMicron)) * Weights[i].D1Err;
		}

		public void UpdateWeightsSGDLMParallel(int batchSize = 1)
		{
			Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= (Network.TrainingRate.Rate / (Weights[i].DiagonalHessian + Network.dMicron)) * Weights[i].D1Err);
		}

		public void UpdateWeightsMiniBatchSGDLMSerial(int batchSize = 1)
		{
			for (int i = 0; i < WeightCount; i++)
				Weights[i].Value -= (Network.TrainingRate.Rate / (Weights[i].DiagonalHessian + Network.dMicron)) * (Weights[i].D1Err / batchSize);
		}

		public void UpdateWeightsMiniBatchSGDLMParallel(int batchSize = 1)
		{
			Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= ((Network.TrainingRate.Rate / (Weights[i].DiagonalHessian + Network.dMicron)) * (Weights[i].D1Err / batchSize)));
		}

		public void UpdateWeighsSGDSerial(int batchSize = 1)
		{
			for (int i = 0; i < WeightCount; i++)
                Weights[i].Value -= (Network.TrainingRate.Rate * Weights[i].D1Err) - (Network.TrainingRate.Rate * Network.TrainingRate.WeightDecayFactor * Weights[i].Value);
		}

		public void UpdateWeightsSGDParallel(int batchSize = 1)
		{
            Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= (Network.TrainingRate.Rate * Weights[i].D1Err) - (Network.TrainingRate.Rate * Network.TrainingRate.WeightDecayFactor * Weights[i].Value));
		}

        public void UpdateWeighsSGDMSerial(int batchSize = 1)
        {
            for (int i = 0; i < WeightCount; i++)
                Weights[i].Value -= (Network.TrainingRate.Momentum * Weights[i].Value) - (Network.TrainingRate.Rate * Weights[i].D1Err);
        }

        public void UpdateWeightsSGDMParallel(int batchSize = 1)
        {
            Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= (Network.TrainingRate.Momentum * Weights[i].Value) - (Network.TrainingRate.Rate * Weights[i].D1Err));
        }

		public void UpdateWeightsMiniBatchSGDSerial(int batchSize = 1)
		{
			for (int i = 0; i < WeightCount; i++)
                Weights[i].Value -= (Network.TrainingRate.Rate * (Weights[i].D1Err / batchSize)) - (Network.TrainingRate.Rate * Network.TrainingRate.WeightDecayFactor * Weights[i].Value);
		}

		public void UpdateWeightsMiniBatchSGDParallel(int batchSize = 1)
		{
            Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= (Network.TrainingRate.Rate * (Weights[i].D1Err / batchSize)) - (Network.TrainingRate.Rate * Network.TrainingRate.WeightDecayFactor * Weights[i].Value));
		}

        public void UpdateWeightsMiniBatchSGDMSerial(int batchSize = 1)
        {
            for (int i = 0; i < WeightCount; i++)
                Weights[i].Value -= (Network.TrainingRate.Momentum * Weights[i].Value) - (Network.TrainingRate.Rate * (Weights[i].D1Err / batchSize));
        }

        public void UpdateWeightsMiniBatchSGDMParallel(int batchSize = 1)
        {
            Parallel.For(0, WeightCount, Network.ParallelOption, i => Weights[i].Value -= (Network.TrainingRate.Momentum * Weights[i].Value) - (Network.TrainingRate.Rate * (Weights[i].D1Err / batchSize)));
        }

		public void NoUpdate(int batchSize = 1)
		{
			// do nothing
		}

		public void NoErase()
		{
			// do nothing
		}

		//public void BackpropagateFullyConnectedSerial()
		//{
		//    for (int i = 0; i < PreviousLayer.NeuronCount; i++)
		//        PreviousLayer.Neurons[i].D1ErrX = 0D;

		//    Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//    {
		//        double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
		//        int idx = i * (PreviousLayer.NeuronCount + 1);
		//        Weights[idx++].D1Err += neuronD1ErrY;
		//        for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//        {
		//            Weights[idx + c].D1Err += neuronD1ErrY * PreviousLayer.Neurons[c].Output;
		//            PreviousLayer.Neurons[c].D1ErrX += neuronD1ErrY * Weights[idx + c].Value;
		//        }
		//    });
		//}

		//public void BackpropagateFullyConnectedParallel()
		//{
		//    Parallel.For(0, PreviousLayer.NeuronCount, i => PreviousLayer.Neurons[i].D1ErrX = 0D);

		//    Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//    {
		//        
		//        double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
		//        int idx = i * (PreviousLayer.NeuronCount + 1);
		//        Weights[idx++].D1Err += neuronD1ErrY;
		//        for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//        {
		//            Weights[idx + c].D1Err += neuronD1ErrY * PreviousLayer.Neurons[c].Output;
		//            PreviousLayer.Neurons[c].D1ErrX += neuronD1ErrY * Weights[idx + c].Value;
		//        }
		//    });
		//}

		public void BackpropagateCCFSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
			});
		}

		public void BackpropagateCCFParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
			});
		}

		//public void BackpropagateFullyConnectedSerial()
		//{
		//    for (int i = 0; i < WeightCount; i++)
		//        Weights[i].D1Err = 0D;

		//    for (int i = 0; i < PreviousLayer.NeuronCount; i++)
		//        PreviousLayer.Neurons[i].D1ErrX = 0D;

		//    Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
		//    {
		//        Neurons[i].D1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
		//        //int idx = i * (PreviousLayer.NeuronCount + 1);
		//        //Weights[idx++].D1Err = Neurons[i].D1ErrY;
		//        //for (int c = 0; c < PreviousLayer.NeuronCount; c++)
		//        //{
		//        //    Weights[idx + c].D1Err = Neurons[i].D1ErrY * PreviousLayer.Neurons[c].Output;
		//        //    PreviousLayer.Neurons[c].D1ErrX += Neurons[i].D1ErrY * Weights[idx + c].Value;
		//        //}

		//        foreach (Connection connection in Connections[i])
		//        {
		//            if (connection.ToNeuronIndex == int.MaxValue)
		//                Weights[connection.ToWeightIndex].D1Err += Neurons[i].D1ErrY;
		//            else
		//            {
		//                Weights[connection.ToWeightIndex].D1Err += Neurons[i].D1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
		//                PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += Neurons[i].D1ErrY * Weights[connection.ToWeightIndex].Value;
		//            }
		//        }
		//    });

		//    for (int i = 0; i < WeightCount; i++)
		//        Weights[i].Value -= Network.TrainingRate.Rate / (Weights[i].DiagonalHessian + Network.dMicron) * Weights[i].D1Err;
		//}
			
		public void BackpropagateCCFDropOutSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				if (NeuronActive[i] == 1)
				{
					double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
						else
						{
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
							PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
						}
					}
				}
			});
		}

		public void BackpropagateCCFDropOutParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				if (NeuronActive[i] == 1)
				{
					double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
						else
						{
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
							PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
						}
					}
				}
			});
		}

		public void BackpropagateAveragePoolingSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / (Connections[i].Length - 1);
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output * sf;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
			});
		}

		public void BackpropagateAveragePoolingParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / (Connections[i].Length - 1);
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * PreviousLayer.Neurons[connection.ToNeuronIndex].Output * sf;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
			});
		}

		public void BackpropagateMaxPoolingSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
				int idx = 0;
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						idx = connection.ToWeightIndex;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
				Weights[idx].D1Err += neuronD1ErrY * PreviousLayer.Neurons[NeuronActive[i]].Output;
			});
		}

		public void BackpropagateMaxPoolingParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
				int idx = 0;
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
					else
					{
						idx = connection.ToWeightIndex;
						PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * Weights[connection.ToWeightIndex].Value;
					}
				}
				Weights[idx].D1Err += neuronD1ErrY * PreviousLayer.Neurons[NeuronActive[i]].Output;
			});
		}

		public void BackpropagateAveragePoolingWeightlessSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / Connections[i].Length;
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				for (int c = 0; c < Connections[i].Length; c++)
					PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].D1ErrX += neuronD1ErrY * sf;
			});
		}

		public void BackpropagateAveragePoolingWeightlessParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / Connections[i].Length;
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;

				for (int c = 0; c < Connections[i].Length; c++)
					PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].D1ErrX += neuronD1ErrY * sf;
			});
		}

		public void BackpropagateMaxPoolingWeightlessSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				PreviousLayer.Neurons[NeuronActive[i]].D1ErrX += DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateMaxPoolingWeightlessParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				PreviousLayer.Neurons[NeuronActive[i]].D1ErrX += DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
			});
		}
		
		public void BackpropagateStochasticPoolingSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				PreviousLayer.Neurons[NeuronActive[i]].D1ErrX += DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateStochasticPoolingParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				PreviousLayer.Neurons[NeuronActive[i]].D1ErrX += DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
			});
		}
	   
		public void BackpropagateL2PoolingSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / Connections[i].Length;
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
				for (int c = 0; c < Connections[i].Length; c++)
					PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].D1ErrX += neuronD1ErrY * sf;
			});
		}

		public void BackpropagateL2PoolingParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / Connections[i].Length;
				double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
				for (int c = 0; c < Connections[i].Length; c++)
					PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].D1ErrX += neuronD1ErrY * sf;
			});
		}

		public void BackpropagateLocalResponseNormalizationSerial()
		{
		   for (int i = 0; i < PreviousLayer.NeuronCount; i++)
			   PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0.0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateLocalResponseNormalizationParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0.0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagationLocalResponseNormalizationCMSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0f;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagationLocalResponseNormalizationCMParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateLocalContrastNormalizationSerial()
		{
			for (int i = 0; i < PreviousLayer.NeuronCount; i++)
				PreviousLayer.Neurons[i].D1ErrX = 0D;
            double[][] w = MathUtil.CreateGaussian2DKernel(ReceptiveFieldWidth, ReceptiveFieldHeight);

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
                double mean = 0D;
                int x = 0;
                int y = 0;
                for (int c = 0; c < Connections[i].Length; c++)
                {
                    if (c % ReceptiveFieldWidth == 0)
                    {
                        x = 0;
                        y++;
                    }
                    mean += PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output * w[x++][y];

                }
                mean /= Connections[i].Length;

				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output - mean);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateLocalContrastNormalizationParallel()
		{
			Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
            double[][] w = MathUtil.CreateGaussian2DKernel(ReceptiveFieldWidth, ReceptiveFieldHeight);

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
                double mean = 0D;
                int x = 0;
                int y = 0;
                for (int c = 0; c < Connections[i].Length; c++)
                {
                    if (c % ReceptiveFieldWidth == 0)
                    {
                        x = 0;
                        y++;
                    }
                    mean += PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output * w[x++][y];

                }
                mean /= Connections[i].Length;

				double sum = 0.0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output - mean);

				PreviousLayer.Neurons[i].D1ErrX += (1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * Neurons[i].D1ErrX;
			});
		}

		public void BackpropagateRBF()
		{
			if (PreviousLayer.UseNeuronPartitioner)
				Parallel.For(0, PreviousLayer.NeuronCount, Network.ParallelOption, i => PreviousLayer.Neurons[i].D1ErrX = 0D);
			else
				for (int i = 0; i < PreviousLayer.NeuronCount; i++)
					PreviousLayer.Neurons[i].D1ErrX = 0D;

			if (ActivationFunctionId == ActivationFunctions.SoftMax)
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					double neuronD1ErrY = DGaussian(Neurons[i].Output) * Neurons[i].D1ErrX;
					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
						else
						{
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * (PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
							PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * (PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
						}
					}
				});
			else
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					double neuronD1ErrY = DerivativeActivationFunction(Neurons[i].Output) * Neurons[i].D1ErrX;
					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY;
						else
						{
							Weights[connection.ToWeightIndex].D1Err += neuronD1ErrY * (PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
							PreviousLayer.Neurons[connection.ToNeuronIndex].D1ErrX += neuronD1ErrY * (PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
						}
					}
				});
		}

		public double[] BackpropagateSecondDerivativesCCFParallel(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];
		   
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output);
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
			});

			Parallel.For(0, WeightCount, Network.ParallelOption, i =>
			{
				Weights[i].DiagonalHessian += weightsD2Err[i];
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesCCFSerial(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output);
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
			});

			for (int i = 0; i < WeightCount; i++)
				Weights[i].DiagonalHessian += weightsD2Err[i];

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesCCFSoftMaxParallel(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];
		   
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				neuronsD2ErrY[i] = neuronsD2ErrX[i];

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output);
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
			});

			Parallel.For(0, WeightCount, Network.ParallelOption, i =>
			{
				Weights[i].DiagonalHessian += weightsD2Err[i];
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesCCFSoftMaxSerial(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];
		   
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				neuronsD2ErrY[i] = neuronsD2ErrX[i];
				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output);
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
			});

			for (int i = 0; i < WeightCount; i++)
				Weights[i].DiagonalHessian += weightsD2Err[i];

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesAveragePooling(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / (Connections[i].Length - 1);
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output * sf);
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
			});

			for (int i = 0; i < WeightCount; i++)
				Weights[i].DiagonalHessian += weightsD2Err[i];
			
			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesMaxPooling(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				int weightIndex = 1;
				double max = double.MinValue;
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];

				foreach (Connection connection in Connections[i])
				{
					if (connection.ToNeuronIndex == int.MaxValue)
						weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
					else
					{
						weightIndex = connection.ToWeightIndex;
						if (max < PreviousLayer.Neurons[connection.ToNeuronIndex].Output)
							max = PreviousLayer.Neurons[connection.ToNeuronIndex].Output;
						prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(Weights[connection.ToWeightIndex].Value);
					}
				}
				weightsD2Err[weightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(max);
			});

			for (int i = 0; i < WeightCount; i++)
				Weights[i].DiagonalHessian += weightsD2Err[i];

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesAveragePoolingWeightless(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sf = 1D / MathUtil.Pow2(Connections[i].Length);
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
				for (int c = 0; c < Connections[i].Length; c++)
					prevLayerNeuronsD2ErrX[Connections[i][c].ToNeuronIndex] += neuronsD2ErrY[i] * sf;
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesMaxPoolingWeightless(double[] neuronsD2ErrX)
		{
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];
		   
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				prevLayerNeuronsD2ErrX[NeuronActive[i]] += MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesStochasticPooling(double[] neuronsD2ErrX)
		{
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				prevLayerNeuronsD2ErrX[NeuronActive[i]] += MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativesL2Pooling(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
				double sf = 1D / MathUtil.Pow2(Connections[i].Length);
				for (int c = 0; c < Connections[i].Length; c++)
					prevLayerNeuronsD2ErrX[Connections[i][c].ToNeuronIndex] += neuronsD2ErrY[i] * sf;
			});

			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativeLocalResponseNormalization(double[] neuronsD2ErrX)
		{
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				prevLayerNeuronsD2ErrX[i] = MathUtil.Pow2(1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * neuronsD2ErrX[i];
			});
		   
			return prevLayerNeuronsD2ErrX;
		}
	   
		public double[] BackpropagateSecondDerivativeLocalResponseNormalizationCM(double[] neuronsD2ErrX)
		{
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output);

				prevLayerNeuronsD2ErrX[i] = MathUtil.Pow2(1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * neuronsD2ErrX[i];
			});
			
			return prevLayerNeuronsD2ErrX;
		}

		public double[] BackpropagateSecondDerivativeLocalContrastNormalization(double[] neuronsD2ErrX)
		{
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];
            double[][] w = MathUtil.CreateGaussian2DKernel(ReceptiveFieldWidth, ReceptiveFieldHeight);
			Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
			{
                double mean = 0D;
                int x = 0;
                int y = 0;
                for (int c = 0; c < Connections[i].Length; c++)
                {
                    if (c % ReceptiveFieldWidth == 0)
                    {
                        x = 0;
                        y++;
                    }
                    mean += PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output * w[x++][y];

                }
                mean /= Connections[i].Length;

				double sum = 0D;
				for (int c = 0; c < Connections[i].Length; c++)
					sum += MathUtil.Pow2(PreviousLayer.Neurons[Connections[i][c].ToNeuronIndex].Output - mean);

				prevLayerNeuronsD2ErrX[i] = MathUtil.Pow2(1D / (double)Math.Pow(((MathUtil.Scale * sum / Connections[i].Length) + 1D), MathUtil.Pow)) * neuronsD2ErrX[i];
			});
			
			return prevLayerNeuronsD2ErrX;
		}
	   
		public double[] BackpropagateSecondDerivativesRBF(double[] neuronsD2ErrX)
		{
			double[] neuronsD2ErrY = new double[NeuronCount];
			double[] weightsD2Err = new double[WeightCount];
			double[] prevLayerNeuronsD2ErrX = new double[PreviousLayer.NeuronCount];

			if (ActivationFunctionId == ActivationFunctions.SoftMax)
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					neuronsD2ErrY[i] = MathUtil.Pow2(DGaussian(Neurons[i].Output)) * neuronsD2ErrX[i];
					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
						else
						{
							weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
							prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
						}
					}
				});
			else
				Parallel.For(0, NeuronCount, Network.ParallelOption, i =>
				{
					neuronsD2ErrY[i] = MathUtil.Pow2(DerivativeActivationFunction(Neurons[i].Output)) * neuronsD2ErrX[i];
					foreach (Connection connection in Connections[i])
					{
						if (connection.ToNeuronIndex == int.MaxValue)
							weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i];
						else
						{
							weightsD2Err[connection.ToWeightIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
							prevLayerNeuronsD2ErrX[connection.ToNeuronIndex] += neuronsD2ErrY[i] * MathUtil.Pow2(PreviousLayer.Neurons[connection.ToNeuronIndex].Output - Weights[connection.ToWeightIndex].Value);
						}
					}
				});

			for (int i = 0; i < WeightCount; i++)
				Weights[i].DiagonalHessian += weightsD2Err[i];
			
			return prevLayerNeuronsD2ErrX;
		}
	}
}