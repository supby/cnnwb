using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CNNWB.Common;

namespace CNNWB.CNN
{
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
        public Func<double, double> DerivativeActivationFunction;
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
                    for (int row = 0; row < ReceptiveFieldHeight; row++)
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
                    for (int row = 0; row < ReceptiveFieldHeight; row++)
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
                    for (int row = 0; row < ReceptiveFieldHeight; row++)
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
                        for (int y = 0; y < NeuronCount; y++)
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
            return MathUtil.Exp(-0.5D * value);
        }

        private static double DGaussian(double value)
        {
            return -(MathUtil.Exp(-0.5D * value) * 0.5D);
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
                double sf = 1D / (Connections[i].Length - 1);
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
                for (int c = 0; c < Connections[i].Length; c++)
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
