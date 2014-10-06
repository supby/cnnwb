cnnwb
=====

Example:

	NeuralNetwork cnn = new NeuralNetwork(DataProvider, "LeNet-5", 10, 0.8D, LossFunctions.MeanSquareError, 
											DataProviderSets.MNIST, TrainingStrategy.SGDLevenbergMarquardt, 0.02D);
	cnn.AddLayer(LayerTypes.Input, 1, 32, 32);
	cnn.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 6, 28, 28, 5, 5);
	cnn.AddLayer(LayerTypes.AveragePooling, ActivationFunctions.Tanh, 6, 14, 14, 2, 2);

	bool[] maps = new bool[6 * 16] 
	{
	 true, false,false,false,true, true, true, false,false,true, true, true, true, false,true, true,
	 true, true, false,false,false,true, true, true, false,false,true, true, true, true, false,true,
	 true, true, true, false,false,false,true, true, true, false,false,true, false,true, true, true,
	 false,true, true, true, false,false,true, true, true, true, false,false,true, false,true, true,
	 false,false,true, true, true, false,false,true, true, true, true, false,true, true, false,true,
	 false,false,false,true, true, true, false,false,true, true, true, true, false,true, true, true
	};

	cnn.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 16, 10, 10, 5, 5, new Mappings(maps));
	cnn.AddLayer(LayerTypes.AveragePooling, ActivationFunctions.Tanh, 16, 5, 5, 2, 2);
	cnn.AddLayer(LayerTypes.Convolutional, ActivationFunctions.Tanh, 120, 1, 1, 5, 5);
	cnn.AddLayer(LayerTypes.FullyConnected, ActivationFunctions.Tanh, 10);
	cnn.InitializeWeights(); 
