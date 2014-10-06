using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CNNWB.Common
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
}
