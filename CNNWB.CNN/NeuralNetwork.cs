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
using CNNWB.Common;
using CNNWB.Data;

namespace CNNWB.CNN
{
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
                    _handler.Invoke(this, e);
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
                    _handler.Invoke(this, e);  
		}

		public NeuralNetwork(DataProvider dataProvider, string name = "Neural Network", int classCount = 10, double trainTo = 0.8D, LossFunctions lossFunction = LossFunctions.MeanSquareError, DataProviderSets dataProviderSet = DataProviderSets.MNIST, TrainingStrategy trainingStrategy = TrainingStrategy.SGDLevenbergMarquardt, double dmicron = 0.02D, int hessianSamples = 2000, bool substractMean = false, double min= -1.0, double max = 1.0)
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
                    throw new Exception("Failed load definition");

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

		public List<List<byte>> GetRBFLeNet5WeightSamples()
		{
            //ResourceManager rm = new ResourceManager("CNNWB.Properties.Resources", Assembly.GetExecutingAssembly());
            //CultureInfo ci = Thread.CurrentThread.CurrentCulture;

            //List<List<byte>> rbfImages = new List<List<byte>>(10);
            //int width = 7;
            //int height = 12;

            //for (int i = 0; i < 10; ++i)
            //{
            //    Image image = rm.GetObject("Image" + i.ToString(), ci) as System.Drawing.Image;
            //    BitmapSource bitmapSource = BitmapToSource(image);
            //    if (bitmapSource.CanFreeze)
            //        bitmapSource.Freeze();

            //    int rawStride = (width * bitmapSource.Format.BitsPerPixel + 7) / 8;
            //    byte[] rawImage = new byte[rawStride * height];
            //    bitmapSource.CopyPixels(rawImage, rawStride, 0);
				
            //    List<byte> inputPattern = new List<byte>(rawStride * 12);
            //    for (int j = 0; j < (rawStride * height); j++)
            //        inputPattern.Add(rawImage[j]);
				
            //    rbfImages.Add(inputPattern);
            //}

            //return rbfImages;

            //throw new NotImplementedException("It took weights from images in resources");
            return null;
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

			OnRaiseProgressEvent(EventArgs.Empty);			
		}

		public async Task StartTraining(CancellationToken cnToken)
		{            
			if (CurrentTaskState == TaskState.Stopped)
			{                   
				workerTimer = new System.Timers.Timer(1000);
				workerTimer.Elapsed += new ElapsedEventHandler(WorkerTimerElapsed);
				TaskDuration.Reset();
				SampleSpeedTimer.Reset();

				CurrentTaskState = TaskState.Running;
				EpochDuration = TimeSpan.Zero;

				TaskDuration.Start();
				SampleSpeedTimer.Start();
				workerTimer.Start();

                await Task.Factory.StartNew(TrainingTask, cnToken).ContinueWith(t => { StopTraining(); });
			}            
		}

		public void StopTraining()
		{
			if (CurrentTaskState != TaskState.Stopped)
			{
				CurrentTaskState = TaskState.Stopped;				
                
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

        public void TrainingTask(object obj)
		{
            CancellationToken token = (CancellationToken)obj;

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

                if (token.IsCancellationRequested)
                {
                    StopTraining();
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

        public async Task StartTesting(CancellationToken cnToken)
		{
			if (CurrentTaskState == TaskState.Stopped)
			{   
				workerTimer = new System.Timers.Timer(1000);
				workerTimer.Elapsed += new ElapsedEventHandler(WorkerTimerElapsed);
				TaskDuration.Reset();
				SampleSpeedTimer.Reset();

				CurrentTaskState = TaskState.Running;			

				TaskDuration.Start();
				SampleSpeedTimer.Start();
				workerTimer.Start();

                await Task.Factory.StartNew(TestingTask, cnToken).ContinueWith(t => { StopTesting(); });
			}
		}

		public void StopTesting()
		{
			if (CurrentTaskState != TaskState.Stopped)
			{
				CurrentTaskState = TaskState.Stopped;

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

        public void TestingTask(object obj)
		{
            CancellationToken token = (CancellationToken)obj;            

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
					OnRaiseTestErrorEvent(new AddUnrecognizedTestSampleEventArgs(SampleIndex, bestIndex, CurrentSample));
				}

				if (CurrentTaskState != TaskState.Running)
					if (!CheckStateChange())
						break;

                if (token.IsCancellationRequested)
                {
                    StopTesting();
                    break;
                }
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
 
	
}