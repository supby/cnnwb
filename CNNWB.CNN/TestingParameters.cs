using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNNWB.CNN
{
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
}
