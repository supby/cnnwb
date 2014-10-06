using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace CNNWB.Data
{
    public enum DataProviderSets
    {
        CIFAR10 = 0,
        MNIST = 1
    }

    public enum CIFAR10Classes
    {
        airplane = 0,
        automobile = 1,
        bird = 2,
        cat = 3,
        deer = 4,
        dog = 5,
        frog = 6,
        horse = 7,
        ship = 8,
        truck = 9
    }

	public class DatasetDownloadInformation
	{
		public DataProviderSets DataProviderSet { get; set; }
		public string Name { get; set; }
		public string Path { get; set; }
		public string Url { get; set; }
		public List<Tuple<string, string>> Files;

		public DatasetDownloadInformation(DataProviderSets dataproviderSet, string name, string path, string url, List<Tuple<string, string>> files)
		{
			DataProviderSet = dataproviderSet;
			Name = name;
			Path = path;
			Url = url;
			Files = files;
		}
	}

	public class DataProviderEventArgs : EventArgs
	{
		public int Result { get; private set; }
		public string Message { get; private set; }
		public TimeSpan Time { get; private set; }

		public DataProviderEventArgs(int result, string msg, TimeSpan time)
		{
			Result = result;
			Message = msg;
			Time = time;
		}
	}

	public class DataProvider : IDisposable
	{
		public event EventHandler<DataProviderEventArgs> RaiseDataProgressEvent = delegate { };
		public event EventHandler<DataProviderEventArgs> RaiseDataLoadedEvent = delegate { };

		public DataProviderSets CurrentDataSet { get; private set; }
		public int ClassCount { get; private set; }
		public int SampleWidth { get; private set; }
		public int SampleHeight { get; private set; }
		public int SampleSize { get; private set; }
		public int SampleChannels { get; private set; }
		public PixelFormat UsedPixelFormat { get; private set; }
		public int GaussianFieldSize { get; private set; }
		public int TrainingSamplesCount { get; private set; }
		public int TestingSamplesCount { get; private set; }
		public int[] RandomTrainingSample;
		public ImageData[] TrainingSamples;
		public ImageData[] TestingSamples;
		public double[] Mean;
		public ParallelOptions ParallelOption { get; private set; }
		public string StorageDirectory { get; private set; }
		public List<String> DataSetPath { get; private set; }
		public List<DatasetDownloadInformation> DataSetDownloadInfo { get; private set; }
		public ThreadSafeRandom RandomGenerator;
        //double SeverityFactor;
        //double MaxScaling;
        //double MaxRotation;
        //double ElasticSigma;
        //double ElasticScaling;
        //public double[] uniformH;
        //public double[] uniformV;
        //public double[] DispH;
        //public double[] DispV;
        //public double[][] GaussianKernel;
        //int Mid;
        //int iMid;
        //double ElasticScale;
        //double AngleFixed;
        //double ScaleFixed;
        //double TwoSigmaSquared;
        //double TwoPiSigma;

		private System.Diagnostics.Stopwatch stopwatch;
		
		private int _maxDegreeOfParallelism = Environment.ProcessorCount;
		public int MaxDegreeOfParallelism
		{
			get
			{
				return _maxDegreeOfParallelism;
			}

			set
			{
				if (value == _maxDegreeOfParallelism)
					return;

				if ((value == 0) || (value > Environment.ProcessorCount))
					_maxDegreeOfParallelism = Environment.ProcessorCount;
				else
					_maxDegreeOfParallelism = value;

				ParallelOption.MaxDegreeOfParallelism = _maxDegreeOfParallelism;
			}
		}

	   //private enum DataType : byte
		//{
		//    typeUnsignedByte = (byte)0x08,
		//    typeSignedByte = (byte)0x09,
		//    typeShort = (byte)0x0B,
		//    typeInteger = (byte)0x0C,
		//    typeFloat = (byte)0x0D,
		//    typeDouble = (byte)0x0E,
		//};

		protected virtual void OnRaiseDataLoadedEvent(DataProviderEventArgs e)
		{
			// Make a temporary copy of the event to avoid possibility of
			// a race condition if the last subscriber unsubscribes
			// immediately after the null check and before the event is raised.
			EventHandler<DataProviderEventArgs> handler = null;

			lock (this)
			{
				handler = RaiseDataLoadedEvent;
			}
			// Event will be null if there are no subscribers

			if (handler != null)
			{
				foreach (EventHandler<DataProviderEventArgs> _handler in handler.GetInvocationList())
				{
					try
					{
						//System.Windows.Application.Current.Dispatcher.Invoke(_handler, new object[] { this, e });
					}
					catch (Exception ex)
					{
						Debug.WriteLine("Error in the handler {0}: {1}", handler.Method.Name, ex.Message);
					}
				}
			}
		}

		protected virtual void OnRaiseDataProgressEvent(DataProviderEventArgs e)
		{
			// Make a temporary copy of the event to avoid possibility of
			// a race condition if the last subscriber unsubscribes
			// immediately after the null check and before the event is raised.
			EventHandler<DataProviderEventArgs> handler = null;

			lock (this)
			{
				handler = RaiseDataProgressEvent;
			}
			// Event will be null if there are no subscribers

			if (handler != null)
			{
				foreach (EventHandler<DataProviderEventArgs> _handler in handler.GetInvocationList())
				{
					try
					{
						//System.Windows.Application.Current.Dispatcher.Invoke(_handler, new object[] { this, e });
					}
					catch (Exception ex)
					{
						Debug.WriteLine("Error in the handler {0}: {1}", handler.Method.Name, ex.Message);
					}
				}
			}
		}

		public DataProvider(string storageDirectory)
		{
			stopwatch = new Stopwatch();
			ParallelOption = new ParallelOptions();
			ParallelOption.TaskScheduler = null;
			ParallelOption.MaxDegreeOfParallelism = Environment.ProcessorCount;
			RandomGenerator = new ThreadSafeRandom();
			
			StorageDirectory = storageDirectory;
			
			DataSetPath = new List<string>();
			DataSetPath.Add(StorageDirectory + @"\CIFAR-10 Dataset");
			DataSetPath.Add(StorageDirectory + @"\MNIST Dataset");
			
			DataSetDownloadInfo = new List<DatasetDownloadInformation>();
			List<Tuple<string, string>> files = new List<Tuple<string, string>>();
			files.Add(new Tuple<string, string>("cifar-10-binary.tar.gz", "All images and labels"));
			DataSetDownloadInfo.Add(new DatasetDownloadInformation(DataProviderSets.CIFAR10, "CIFAR-10 Dataset", DataSetPath[(int)DataProviderSets.CIFAR10], @"http://www.cs.toronto.edu/~kriz/", files));
			files = new List<Tuple<string, string>>();
			files.Add(new Tuple<string, string>("train-images-idx3-ubyte.gz", "Training images"));
			files.Add(new Tuple<string, string>("t10k-images-idx3-ubyte.gz", "Testing images"));
			files.Add(new Tuple<string, string>("train-labels-idx1-ubyte.gz", "Training labels"));
			files.Add(new Tuple<string, string>("t10k-labels-idx1-ubyte.gz", "Testing labels"));
			DataSetDownloadInfo.Add(new DatasetDownloadInformation(DataProviderSets.MNIST, "MNIST Dataset", DataSetPath[(int)DataProviderSets.MNIST], @"http://yann.lecun.com/exdb/mnist/", files));
			
			// Create the directories if necessary
			foreach(string path in DataSetPath)
				if (!Directory.Exists(path))
					Directory.CreateDirectory(path);
		}

		public bool DataSetFilesAvailable(DataProviderSets set)
		{
			bool available = true;
			string path = DataSetPath[(int)set];

			switch (set)
			{
				case DataProviderSets.CIFAR10:
					if (!File.Exists(path + @"\data_batch_1.bin") || !File.Exists(path + @"\data_batch_2.bin") || !File.Exists(path + @"\data_batch_3.bin") || !File.Exists(path + @"\data_batch_4.bin") || !File.Exists(path + @"\data_batch_5.bin") || !File.Exists(path + @"\test_batch.bin") || !File.Exists(path + @"\batches.meta.txt"))
						available = false;
					break;

				case DataProviderSets.MNIST:
					if (!File.Exists(path + @"\t10k-labels-idx1-ubyte") || !File.Exists(path + @"\train-labels-idx1-ubyte") || !File.Exists(path + @"\t10k-images-idx3-ubyte") || !File.Exists(path + @"\train-images-idx3-ubyte"))
						available = false;
					break;
			}

			return available;
		}

		public bool AllDataSetFilesAvailable()
		{
			return DataSetFilesAvailable(DataProviderSets.CIFAR10) && DataSetFilesAvailable(DataProviderSets.MNIST);
		}

		~DataProvider()
		{
			// In case the client forgets to call
			// Dispose , destructor will be invoked for
			Dispose(false);
		}

		protected virtual void Dispose(bool disposing)
		{
			if (disposing)
			{
				// dispose managed resources
				if (RandomGenerator != null)
					RandomGenerator.Dispose();
				RandomGenerator = null;
				TrainingSamples = null;
				TrainingSamples = null;
				RandomTrainingSample = null;
				ParallelOption = null;
				stopwatch = null;
				RaiseDataLoadedEvent = null;
				RaiseDataProgressEvent = null;
			}
			// free native resources
		}

		public void Dispose()
		{
			Dispose(true);
			GC.SuppressFinalize(this);
		}

		public void LoadDataset(DataProviderSets set)
		{
			switch (set)
			{
				case DataProviderSets.CIFAR10:
					Task.Factory.StartNew(() => LoadCIFAR10DataSet());
					break;

				case DataProviderSets.MNIST:
					Task.Factory.StartNew(() => LoadMNISTDataSet());
					break;
			}
		}

		public void ScrambleTrainingSamples()
		{
			int k, l;
			for (int i = 0; i < TrainingSamplesCount; i++)
			{
				l = RandomGenerator.Next(TrainingSamplesCount);
				k = RandomTrainingSample[i];
				RandomTrainingSample[i] = RandomTrainingSample[l];
				RandomTrainingSample[l] = k;
			}
		}

		public void LoadMNISTDataSet()
		{
			stopwatch.Restart();
			ClassCount = 10;
			SampleWidth = 32;
			SampleHeight = 32;
			SampleSize = SampleWidth * SampleHeight;
			SampleChannels = 1;
			UsedPixelFormat = PixelFormats.Gray8;
			GaussianFieldSize = 25;

			int counter = 0;
			List<Task> tasks = new List<Task>();
			tasks.Add(new Task(() => LoadMNISTTrainingSamples(ref counter)));
			tasks.Add(new Task(() => LoadMNISTTestingSamples(ref counter)));

			foreach (Task task in tasks)
				task.Start();

			while (!tasks.TrueForAll(task => task.IsCompleted))
			{
				OnRaiseDataProgressEvent(new DataProviderEventArgs((int)(counter), "MNIST Dataset loading...", stopwatch.Elapsed));
				Thread.Sleep(10);
			}

			foreach (Task task in tasks)
				task.Dispose();
			tasks = null;

			OnRaiseDataLoadedEvent(new DataProviderEventArgs(100, String.Format("Dataset loaded at {0}", DateTime.Now.ToString()), stopwatch.Elapsed));
			CurrentDataSet = DataProviderSets.MNIST;
			stopwatch.Reset();

			GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
			GC.WaitForPendingFinalizers();
			GC.Collect();
		}

		private void LoadMNISTTrainingSamples(ref int counter)
		{
			string path = DataSetPath[(int)DataProviderSets.MNIST];
			int MNistWidth;
			int MNistHeight;
			int MNistSize;
			try
			{

				byte[] TrainPatterns;
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\train-images-idx3-ubyte", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					byte[] magicNumber = new byte[16];
					fileStreamTrainPatterns.Read(magicNumber, 0, 16);
					//DataType dataType = (DataType)magicNumber[2];
					//int dimensions = (int)magicNumber[3];

					TrainingSamplesCount = magicNumber[7] + (magicNumber[6] * 256) + (magicNumber[5] * 65536) + (magicNumber[4] * 16777216);
					MNistHeight = magicNumber[11] + (magicNumber[10] * 256) + (magicNumber[9] * 65536) + (magicNumber[8] * 16777216);
					MNistWidth = magicNumber[15] + (magicNumber[14] * 256) + (magicNumber[13] * 65536) + (magicNumber[12] * 16777216);
					MNistSize = MNistWidth * MNistHeight;

					TrainPatterns = new byte[TrainingSamplesCount * MNistSize];
					fileStreamTrainPatterns.Seek(16, SeekOrigin.Begin);
					fileStreamTrainPatterns.Read(TrainPatterns, 0, TrainingSamplesCount * MNistSize);
				}

				byte[] TrainLabels;
				using (FileStream fileStreamTrainLabels = System.IO.File.Open(path + @"\train-labels-idx1-ubyte", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					TrainLabels = new byte[TrainingSamplesCount];
					fileStreamTrainLabels.Seek(8, SeekOrigin.Begin);
					fileStreamTrainLabels.Read(TrainLabels, 0, TrainingSamplesCount);
				}

				int xOffset = (SampleWidth - MNistWidth) / 2;
				int yOffset = (SampleHeight - MNistHeight) / 2;
				TrainingSamples = new ImageData[TrainingSamplesCount];
				RandomTrainingSample = new int[TrainingSamplesCount];
				counter += 25;
				Parallel.For(0, TrainingSamplesCount, ParallelOption, j =>
				{
					RandomTrainingSample[j] = j;

					ImageData sample = new ImageData(TrainLabels[j], new byte[SampleSize]);
					for (int i = 0; i < SampleSize; i++)
						sample.Image[i] = 0;

					for (int y = 0; y < MNistHeight; y++)
						for (int x = 0; x < MNistWidth; x++)
							sample.Image[(x + xOffset) + ((y + yOffset) * SampleWidth)] = TrainPatterns[x + (y * MNistWidth) + (j * MNistSize)];

					TrainingSamples[j] = sample; //.InvertColors(PatternWidth, PatternHeight, PatternChannels);
				});
				counter += 25;
				
				// Fix wrongly assigned numbers in the MNIST training dataset
				TrainingSamples[7080].Label = 5;
				TrainingSamples[10994].Label = 9;
				TrainingSamples[30751].Label = 0;
				TrainingSamples[35310].Label = 6;
				TrainingSamples[40144].Label = 3;
				TrainingSamples[43454].Label = 3;
				TrainingSamples[59915].Label = 7;
			}
			catch (Exception)
			{
				throw;
			}
		}

		private void LoadMNISTTestingSamples(ref int counter)
		{
			string path = DataSetPath[(int)DataProviderSets.MNIST];
			int MNistWidth;
			int MNistHeight;
			int MNistSize;
			try
			{
				byte[] TestPatterns;
				using (FileStream fileStreamTestPatterns = System.IO.File.Open(path + @"\t10k-images-idx3-ubyte", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					byte[] magicNumber = new byte[16];
					fileStreamTestPatterns.Read(magicNumber, 0, 16);
					//DataType dataType = (DataType)magicNumber[2];
					//int dimensions = (int)magicNumber[3];

					TestingSamplesCount = magicNumber[7] + (magicNumber[6] * 256) + (magicNumber[5] * 65536) + (magicNumber[4] * 16777216);
					MNistHeight = magicNumber[11] + (magicNumber[10] * 256) + (magicNumber[9] * 65536) + (magicNumber[8] * 16777216);
					MNistWidth = magicNumber[15] + (magicNumber[14] * 256) + (magicNumber[13] * 65536) + (magicNumber[12] * 16777216);
					MNistSize = MNistWidth * MNistHeight;

					TestPatterns = new byte[TestingSamplesCount * MNistSize];
					fileStreamTestPatterns.Seek(16, SeekOrigin.Begin);
					fileStreamTestPatterns.Read(TestPatterns, 0, TestingSamplesCount * MNistSize);
				}

				byte[] TestLabels;
				using (FileStream fileStreamTestLabels = System.IO.File.Open(path + @"\t10k-labels-idx1-ubyte", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					TestLabels = new byte[TestingSamplesCount];
					fileStreamTestLabels.Seek(8, SeekOrigin.Begin);
					fileStreamTestLabels.Read(TestLabels, 0, TestingSamplesCount);
				}

				int xOffset = (SampleWidth - MNistWidth) / 2;
				int yOffset = (SampleHeight - MNistHeight) / 2;
				TestingSamples = new ImageData[TestingSamplesCount];
				counter += 25;
				Parallel.For(0, TestingSamplesCount, ParallelOption, j =>
				{
					ImageData sample = new ImageData(TestLabels[j], new byte[SampleSize]);
					for (int i = 0; i < SampleSize; i++)
						sample.Image[i] = 0;

					for (int y = 0; y < MNistHeight; y++)
						for (int x = 0; x < MNistWidth; x++)
							sample.Image[(x + xOffset) + ((y + yOffset) * SampleWidth)] = TestPatterns[x + (y * MNistWidth) + (j * MNistSize)];

					TestingSamples[j] = sample; //.InvertColors(PatternWidth, PatternHeight, PatternChannels);
				});
				counter += 25;
			}
			catch (Exception)
			{
				throw;
			}
		}

		public void LoadCIFAR10DataSet(bool normalizeContrast = false)
		{
			string path = DataSetPath[(int)DataProviderSets.CIFAR10];

			stopwatch.Restart();
			SampleWidth = 32;
			SampleHeight = 32;
			SampleSize = SampleWidth * SampleHeight;
			SampleChannels = 3;
			UsedPixelFormat = PixelFormats.Rgb24;
			ClassCount = 10;
			GaussianFieldSize = 25;

			int counter = 0;
			ImageData[] CIFAR10Training = new ImageData[50000];
			ImageData[] CIFAR10Testing = new ImageData[10000];

			List<Task> tasks = new List<Task>(6);

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\data_batch_1.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 0; index < 10000; index++)
				{
					CIFAR10Training[index].Label = samples[3073 * index];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * index) + 1) + i];
					}
					CIFAR10Training[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\data_batch_2.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 10000; index < 20000; index++)
				{
					CIFAR10Training[index].Label = samples[3073 * (index - 10000)];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * (index - 10000)) + 1) + i];
					}
					CIFAR10Training[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\data_batch_3.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 20000; index < 30000; index++)
				{
					CIFAR10Training[index].Label = samples[3073 * (index - 20000)];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * (index - 20000)) + 1) + i];
					}
					CIFAR10Training[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\data_batch_4.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 30000; index < 40000; index++)
				{
					CIFAR10Training[index].Label = samples[3073 * (index - 30000)];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * (index - 30000)) + 1) + i];
					}
					CIFAR10Training[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\data_batch_5.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 40000; index < 50000; index++)
				{
					CIFAR10Training[index].Label = samples[3073 * (index - 40000)];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * (index - 40000)) + 1) + i];
					}
					CIFAR10Training[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			tasks.Add(new Task(() =>
			{
				byte[] samples = new byte[30730000];
				using (FileStream fileStreamTrainPatterns = System.IO.File.Open(path + @"\test_batch.bin", FileMode.Open, FileAccess.Read, FileShare.Read))
				{
					fileStreamTrainPatterns.Read(samples, 0, 30730000);
				}
				for (int index = 0; index < 10000; index++)
				{
					CIFAR10Testing[index].Label = samples[3073 * index];
					byte[] image = new byte[3072];
					for (int i = 0; i < 3072; i++)
					{
						image[i] = samples[((3073 * index) + 1) + i];
					}
					CIFAR10Testing[index].Image = image;
					System.Threading.Interlocked.Increment(ref counter);
				}
			}));

			foreach (Task task in tasks)
				task.Start();

			while (!tasks.TrueForAll(task => task.IsCompleted))
			{
				OnRaiseDataProgressEvent(new DataProviderEventArgs((int)(counter / 1200), "CIFAR-10 Dataset loading...", stopwatch.Elapsed));
				Thread.Sleep(10);
			}

			foreach (Task task in tasks)
				task.Dispose();
			tasks = null;

			OnRaiseDataProgressEvent(new DataProviderEventArgs(50, "Creating Training Samples...", stopwatch.Elapsed));

			// substractive mean
			Mean = new double[SampleSize*SampleChannels];

			TrainingSamplesCount = 100000;
			RandomTrainingSample = new int[TrainingSamplesCount];
			TrainingSamples = new ImageData[TrainingSamplesCount];
			int part = TrainingSamplesCount / 2;

            for (int i = 0; i < SampleSize; i++)
                for (int c = 0; c < SampleChannels; c++)
                    Mean[i + (c * SampleSize)] = 0D;

            OnRaiseDataProgressEvent(new DataProviderEventArgs(70, String.Format("Creating Training Samples..."), stopwatch.Elapsed));

			Parallel.For(0, part, ParallelOption, j =>
			{
				RandomTrainingSample[j] = j;
				TrainingSamples[j] = CIFAR10Training[j];
			});

			Parallel.For(part, 2 * part, ParallelOption, j =>
			{
				RandomTrainingSample[j] = j;
				TrainingSamples[j] = TrainingSamples[j - part].Mirror(32, 32, 3);
			});

            Parallel.For(0, part, ParallelOption, j =>
            {
                for (int i = 0; i < SampleSize; i++)
                    for (int c = 0; c < SampleChannels; c++)
                        Mean[i + (c * SampleSize)] += (double)TrainingSamples[j].Image[i + (c * SampleSize)] / 255D;
            });

            Parallel.For(part, 2 * part, ParallelOption, j =>
            {
                for (int i = 0; i < SampleSize; i++)
                    for (int c = 0; c < SampleChannels; c++)
                        Mean[i + (c * SampleSize)] += (double)TrainingSamples[j].Image[i + (c * SampleSize)] / 255D;
            });

            for (int i = 0; i < SampleSize; i++)
                for (int c = 0; c < SampleChannels; c++)
                    Mean[i + (c * SampleSize)] /= (double)TrainingSamplesCount;

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(60, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(2 * part, 3 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (2 * part)].GetReducedSample(RGBImageData.Position.TopRight, 24, 24).Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(65, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(3 * part, 4 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (3 * part)].GetReducedSample(RGBImageData.Position.TopLeft, 24, 24).GetMirroredImage().Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(70, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(4 * part, 5 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (4 * part)].GetReducedSample(RGBImageData.Position.BottomRight, 24, 24).Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(75, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(5 * part, 6 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (5 * part)].GetReducedSample(RGBImageData.Position.BottomLeft, 24, 24).GetMirroredImage().Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(80, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(6 * part, 7 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (6 * part)].GetReducedSample(RGBImageData.Position.Center, 24, 24).GetMirroredImage().Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			//OnRaiseDataProgressEvent(new DataProviderEventArgs(85, String.Format("Creating Training Patterns..."), stopwatch.Elapsed));
			//Parallel.For(7 * part, 8 * part, ParallelOption, j =>
			//{
			//    RandomTrainingPattern[j] = j;
			//    CIFAR10TrainingPatterns[j] = CIFAR10TrainingPatterns[j - (7 * part)].GetReducedSample(RGBImageData.Position.Center, 24, 24).Resize(32, 32, RGBImageData.Interpolation.Bilinear);
			//});

			OnRaiseDataProgressEvent(new DataProviderEventArgs(90, "Creating Testing Samples...", stopwatch.Elapsed));
			TestingSamplesCount = 10000;
			TestingSamples = new ImageData[TestingSamplesCount];
			Parallel.For(0, TestingSamplesCount, ParallelOption, j =>
			{
				TestingSamples[j] = CIFAR10Testing[j];
			});

			CIFAR10Training = null;
			CIFAR10Testing = null;
			CurrentDataSet = DataProviderSets.CIFAR10;

			OnRaiseDataLoadedEvent(new DataProviderEventArgs(100, String.Format("Dataset loaded at {0}", DateTime.Now.ToString()), stopwatch.Elapsed));

			stopwatch.Stop();
			GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced);
			GC.WaitForPendingFinalizers();
			GC.Collect();
			stopwatch.Reset();
		}
	}
}