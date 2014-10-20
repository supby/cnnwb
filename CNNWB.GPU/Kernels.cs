using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Cudafy;
using Cudafy.Host;
using Cudafy.Translator;

namespace CNNWB.GPU
{
    class Kernels
    {
        [Cudafy]
        public static void CalculateFullyConnected(GThread thread, double[] weights, double[] previousLayerNeurons, double[] layerNeurons)
        {
            int bid = thread.blockIdx.x;
            int tid = thread.threadIdx.x;

            //int i = thread.blockDim.x * bid + tid;

            //c[tid] += a[i] + b[tid];
        }
    }
}
