using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNNWB.CNN
{
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
        public bool Distorted { get; set; }
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
}
