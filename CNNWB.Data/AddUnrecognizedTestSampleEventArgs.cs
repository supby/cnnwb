using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNNWB.Data
{
    public sealed class AddUnrecognizedTestSampleEventArgs : EventArgs
    {
        public int SampleIndex { get; private set; }
        public int WrongValue { get; private set; }
        public ImageData Sample { get; private set; }

        public AddUnrecognizedTestSampleEventArgs(int sampleIndex, int wrongValue, ImageData sample)
        {
            SampleIndex = sampleIndex;
            WrongValue = wrongValue;
            Sample = sample;
        }
    }
}
