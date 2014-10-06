using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using CNNWB.Common;

namespace CNNWB.Data
{
    public struct ImageData
    {
        public int Label;
        public byte[] Image;

        public ImageData(int label, byte[] image)
            : this()
        {
            Label = label;
            Image = image;
        }

        public enum Position
        {
            TopLeft = 0,
            TopRight = 1,
            BottomLeft = 2,
            BottomRight = 3,
            Center = 4
        }

        public enum Interpolation
        {
            NearestNeighbor = 0,
            Bilinear = 1
        }

        public ImageData Resize(int sourceWidth, int sourceHeight, int destWidth, int destHeight, int channels, Interpolation interpolation)
        {
            int sourceSize = sourceWidth * sourceHeight;
            int destSize = destWidth * destHeight;
            byte[] destImage = new byte[destSize * channels];

            double xs = (double)sourceWidth / destWidth;
            double ys = (double)sourceHeight / destHeight;
            double sx, sy;
            int x0, y0;

            if (interpolation == Interpolation.NearestNeighbor)
            {
                for (int y = 0; y < destHeight; y++)
                    for (int x = 0; x < destWidth; x++)
                    {
                        sx = x * xs;
                        sy = y * ys;
                        x0 = (int)sx;
                        y0 = (int)sy;

                        for (int c = 0; c < channels; c++)
                            destImage[x + (y * destWidth) + (c * destSize)] = Image[x0 + (y0 * sourceWidth) + (c * sourceSize)];
                    }
            }
            else if (interpolation == Interpolation.Bilinear)
            {
                double fracx, fracy, ifracx, ifracy;
                double l0, l1;
                double[] cf = new double[channels];
                int x1, y1;
                byte[] c1 = new byte[channels];
                byte[] c2 = new byte[channels];
                byte[] c3 = new byte[channels];
                byte[] c4 = new byte[channels];

                for (int y = 0; y < destHeight; y++)
                    for (int x = 0; x < destWidth; x++)
                    {
                        sx = x * xs;
                        sy = y * ys;
                        x0 = (int)sx;
                        y0 = (int)sy;

                        // Calculate coordinates of the 4 interpolation points
                        fracx = sx - x0;
                        fracy = sy - y0;
                        ifracx = 1f - fracx;
                        ifracy = 1f - fracy;

                        x1 = x0 + 1;
                        if (x1 >= sourceWidth)
                            x1 = x0;

                        y1 = y0 + 1;
                        if (y1 >= sourceHeight)
                            y1 = y0;

                        for (int c = 0; c < channels; c++)
                        {
                            c1[c] = Image[(y0 * sourceWidth + x0) + (c * sourceSize)];
                            c2[c] = Image[(y0 * sourceWidth + x1) + (c * sourceSize)];
                            c3[c] = Image[(y1 * sourceWidth + x0) + (c * sourceSize)];
                            c4[c] = Image[(y1 * sourceWidth + x1) + (c * sourceSize)];

                            l0 = ifracx * c1[c] * 255 + fracx * c2[c] * 255;
                            l1 = ifracx * c3[c] * 255 + fracx * c4[c] * 255;
                            cf[c] = ifracy * l0 + fracy * l1;

                            // Divide by alpha
                            cf[c] /= 255;

                            destImage[x + (y * destWidth) + (c * destSize)] = (byte)cf[c];
                        }
                    }
            }

            return (new ImageData(Label, destImage));
        }

        public ImageData Resize(DataProvider dataProvider, int destWidth, int destHeight, Interpolation interpolation)
        {
            return Resize(dataProvider.SampleWidth, dataProvider.SampleHeight, destWidth, destHeight, dataProvider.SampleChannels, interpolation);
        }

        public ImageData InvertColors(int width, int height, int channels)
        {
            int size = width * height;
            byte[] destImage = new byte[size * channels];

            for (int c = 0; c < channels; c++)
                for (int i = 0; i < size; i++)
                    destImage[i + (c * size)] = (byte)(255 - Image[i + (c * size)]);

            return (new ImageData(Label, destImage));
        }

        public ImageData InvertColors(DataProvider dataProvider)
        {
            return InvertColors(dataProvider.SampleWidth, dataProvider.SampleHeight, dataProvider.SampleChannels);
        }

        public double[] SubstractMean(DataProvider dataProvider)
        {
            double[] meanSubstractedSample = new double[dataProvider.SampleSize * dataProvider.SampleChannels];

            for (int i = 0; i < dataProvider.SampleSize; i++)
                for (int c = 0; c < dataProvider.SampleChannels; c++)
                    meanSubstractedSample[i + (c * dataProvider.SampleSize)] = ((double)Image[i + (c * dataProvider.SampleSize)] / 255D) - dataProvider.Mean[i + (c * dataProvider.SampleSize)];

            return meanSubstractedSample;
        }

        public ImageData Mirror(int width, int height, int channels)
        {
            int size = width * height;
            byte[] destImage = new byte[size * channels];

            for (int c = 0; c < channels; c++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        destImage[x + (y * width) + (c * size)] = Image[((y * width) + width - 1) - x + (c * size)];

            return new ImageData(Label, destImage);
        }

        public ImageData Mirror(DataProvider dataProvider)
        {
            return Mirror(dataProvider.SampleWidth, dataProvider.SampleHeight, dataProvider.SampleChannels);
        }

        public ImageData Crop(Position position, int sourceWidth, int sourceHeight, int destWidth, int destHeight, int channels)
        {
            int sourceSize = sourceWidth * sourceHeight;
            int destSize = destWidth * destHeight;
            int dX = sourceWidth - destWidth;
            int dY = sourceHeight - destHeight;
            byte[] destImage = new byte[destSize * channels];

            switch (position)
            {
                case Position.TopLeft:
                    for (int c = 0; c < channels; c++)
                        for (int y = 0; y < destHeight; y++)
                            for (int x = 0; x < destWidth; x++)
                                destImage[x + (y * destWidth) + (c * destSize)] = Image[x + (y * sourceWidth) + (c * sourceSize)];
                    break;

                case Position.TopRight:
                    for (int c = 0; c < channels; c++)
                        for (int y = 0; y < destHeight; y++)
                            for (int x = 0; x < destWidth; x++)
                                destImage[x + (y * destWidth) + (c * destSize)] = Image[x + (dX - 1 + (y * sourceWidth)) + (c * sourceSize)];
                    break;

                case Position.BottomLeft:
                    for (int c = 0; c < channels; c++)
                        for (int y = 0; y < destHeight; y++)
                            for (int x = 0; x < destWidth; x++)
                                destImage[x + (y * destWidth) + (c * destSize)] = Image[x + ((y + dY - 1) * sourceWidth) + (c * sourceSize)];
                    break;

                case Position.BottomRight:
                    for (int c = 0; c < channels; c++)
                        for (int y = 0; y < destHeight; y++)
                            for (int x = 0; x < destWidth; x++)
                                destImage[x + (y * destWidth) + (c * destSize)] = Image[x + (dX - 1 + ((y + dY - 1) * sourceWidth)) + (c * sourceSize)];
                    break;

                case Position.Center:
                    for (int c = 0; c < channels; c++)
                        for (int y = 0; y < destHeight; y++)
                            for (int x = 0; x < destWidth; x++)
                                destImage[x + (y * destWidth) + (c * destSize)] = Image[x + ((dX / 2) - 1 + ((y + (dY / 2) - 1) * sourceWidth)) + (c * sourceSize)];
                    break;
            }

            return new ImageData(Label, destImage);
        }

        public ImageData Crop(DataProvider dataProvider, Position position, int destWidth, int destHeight)
        {
            return Crop(position, dataProvider.SampleWidth, dataProvider.SampleHeight, destWidth, destHeight, dataProvider.SampleChannels);
        }

        public ImageData Distorted(DataProvider dataProvider, double severityFactor, double maxScaling, double maxRotation, double elasticSigma, double elasticScaling)
        {
            int Mid = dataProvider.GaussianFieldSize / 2;
            int iMid = dataProvider.SampleHeight / 2;
            double elasticScale = severityFactor * elasticScaling;
            double angleFixed = severityFactor * maxRotation * Math.PI / 180.0D;
            double scaleFixed = severityFactor * maxScaling * 0.01D;
            double twoSigmaSquared = 1D / (2D * elasticSigma * elasticSigma);
            double twoPiSigma = 1D / (2D * Math.PI * elasticSigma * elasticSigma);
            //double twoPiSigma = 1D / elasticSigma * Math.Sqrt(2D * Math.PI);
            double[] uniformH = new double[dataProvider.SampleSize];
            double[] uniformV = new double[dataProvider.SampleSize];
            double[] DispH = new double[dataProvider.SampleSize];
            double[] DispV = new double[dataProvider.SampleSize];
            double[] GaussianKernel = new double[dataProvider.GaussianFieldSize * dataProvider.GaussianFieldSize];
            byte[] destImage = new byte[dataProvider.SampleSize * dataProvider.SampleChannels];
            byte[] sourceImage = Image;

            Parallel.For(0, dataProvider.GaussianFieldSize, row =>
            {
                for (int col = 0; col < dataProvider.GaussianFieldSize; col++)
                    GaussianKernel[(row * dataProvider.GaussianFieldSize) + col] = twoPiSigma * (MathUtil.Exp(-((((row - Mid) * (row - Mid)) + ((col - Mid) * (col - Mid))) * twoSigmaSquared)));
            });

            for (int i = 0; i < dataProvider.SampleSize; i++)
            {
                uniformH[i] = dataProvider.RandomGenerator.NextDouble(1);
                uniformV[i] = dataProvider.RandomGenerator.NextDouble(1);
            }

            Parallel.For(0, dataProvider.SampleWidth, col =>
            {
                double fSampleH;
                double fSampleV;
                int xxxDisp;
                int yyyDisp;
                for (int row = 0; row < dataProvider.SampleHeight; row++)
                {
                    double fConvolvedH = 0.0f;
                    double fConvolvedV = 0.0f;

                    for (int xxx = 0; xxx < dataProvider.GaussianFieldSize; xxx++)
                        for (int yyy = 0; yyy < dataProvider.GaussianFieldSize; yyy++)
                        {
                            xxxDisp = col - Mid + xxx;
                            yyyDisp = row - Mid + yyy;

                            if (xxxDisp < 0 || xxxDisp >= dataProvider.SampleWidth || yyyDisp < 0 || yyyDisp >= dataProvider.SampleHeight)
                            {
                                fSampleH = 0.0D;
                                fSampleV = 0.0D;
                            }
                            else
                            {
                                fSampleH = uniformH[xxxDisp + (yyyDisp * dataProvider.SampleWidth)];
                                fSampleV = uniformV[xxxDisp + (yyyDisp * dataProvider.SampleWidth)];
                            }

                            fConvolvedH += fSampleH * GaussianKernel[yyy * (dataProvider.GaussianFieldSize) + xxx];
                            fConvolvedV += fSampleV * GaussianKernel[yyy * (dataProvider.GaussianFieldSize) + xxx];
                        }

                    DispH[col + (row * dataProvider.SampleWidth)] = elasticScale * fConvolvedH;
                    DispV[col + (row * dataProvider.SampleWidth)] = elasticScale * fConvolvedV;
                }
            });

            // next, the scaling of the image by a random scale factor
            // Horizontal and vertical directions are scaled independently
            double dSFHoriz = scaleFixed * dataProvider.RandomGenerator.NextDouble(1);  // MaxScaling is a percentage
            double dSFVert = scaleFixed * dataProvider.RandomGenerator.NextDouble(1);   // MaxScaling is a percentage
            double angle = angleFixed * dataProvider.RandomGenerator.NextDouble(1);
            double cosAngle = Math.Cos(angle);
            double sinAngle = Math.Sin(angle);
            Parallel.For(0, dataProvider.SampleHeight, row =>
            {
                for (int col = 0; col < dataProvider.SampleWidth; col++)
                {
                    DispH[col + (row * dataProvider.SampleWidth)] += dSFHoriz * (col - iMid);
                    DispV[col + (row * dataProvider.SampleWidth)] -= dSFVert * (iMid - row);  // negative because of top-down bitmap
                    // Apply a rotation
                    DispH[col + (row * dataProvider.SampleWidth)] += (col - iMid) * (cosAngle - 1) - (iMid - row) * sinAngle;
                    DispV[col + (row * dataProvider.SampleWidth)] -= (iMid - row) * (cosAngle - 1) + (col - iMid) * sinAngle;  // negative because of top-down bitmap
                }
            });

            Parallel.For(0, dataProvider.SampleWidth, col =>
            {
                for (int row = 0; row < dataProvider.SampleHeight; row++)
                    for (int c = 0; c < dataProvider.SampleChannels; c++)
                    {
                        double sourceRow = (double)row - DispV[col + (row * dataProvider.SampleWidth)];
                        double sourceCol = (double)col - DispH[col + (row * dataProvider.SampleWidth)];

                        double fracRow = sourceRow - (int)sourceRow;
                        double fracCol = sourceCol - (int)sourceCol;

                        double w1 = (1.0D - fracRow) * (1.0D - fracCol);
                        double w2 = (1.0D - fracRow) * fracCol;
                        double w3 = fracRow * (1.0D - fracCol);
                        double w4 = fracRow * fracCol;

                        bool skipOutOfBounds = false;

                        if ((sourceRow + 1.0D) >= dataProvider.SampleHeight) skipOutOfBounds = true;
                        if (sourceRow < 0D) skipOutOfBounds = true;

                        if ((sourceCol + 1.0D) >= dataProvider.SampleWidth) skipOutOfBounds = true;
                        if (sourceCol < 0D) skipOutOfBounds = true;

                        double sourceValue = -1.0D;   // Background
                        if (!skipOutOfBounds)
                        {
                            // the supporting pixels for the "phantom" source pixel are all within the 
                            // bounds of the character grid.
                            // Manufacture its value by bi-linear interpolation of surrounding pixels

                            int sRow = (int)sourceRow;
                            int sCol = (int)sourceCol;

                            int sRowp1 = sRow + 1;
                            int sColp1 = sCol + 1;

                            while (sRowp1 >= dataProvider.SampleHeight) sRowp1 -= dataProvider.SampleHeight;
                            while (sRowp1 < 0) sRowp1 += dataProvider.SampleHeight;

                            while (sColp1 >= dataProvider.SampleWidth) sColp1 -= dataProvider.SampleWidth;
                            while (sColp1 < 0) sColp1 += dataProvider.SampleWidth;

                            // perform bi-linear interpolation
                            sourceValue = (w1 * (double)sourceImage[sCol + (sRow * dataProvider.SampleWidth) + (c * dataProvider.SampleSize)]) + (w2 * (double)sourceImage[(sRow * dataProvider.SampleWidth) + sColp1 + (c * dataProvider.SampleSize)]) + (w3 * (double)sourceImage[(sRowp1 * dataProvider.SampleWidth) + sCol + (c * dataProvider.SampleSize)]) + (w4 * (double)sourceImage[(sRowp1 * dataProvider.SampleWidth) + sColp1 + (c * dataProvider.SampleSize)]);
                        }
                        else
                            sourceValue = sourceImage[col + (row * dataProvider.SampleWidth) + (c * dataProvider.SampleSize)];

                        destImage[col + (row * dataProvider.SampleWidth) + (c * dataProvider.SampleSize)] = (byte)sourceValue;
                    }
            });

            return new ImageData(Label, destImage);
        }

        public BitmapSource GetBitmapSource(int width, int height, int channels, PixelFormat pixelFormat)
        {
            int size = width * height;
            int rawStride = (width * pixelFormat.BitsPerPixel + 7) / 8;  // should be the same value as channels * width
            byte[] rawImage = new byte[rawStride * height];

            for (int c = 0; c < channels; c++)
                for (int i = 0; i < size; i++)
                    rawImage[(i * channels) + c] = Image[i + (c * size)];

            BitmapSource output = BitmapSource.Create(width, height, 96, 96, pixelFormat, null, rawImage, rawStride);
            if (output.CanFreeze)
                output.Freeze();

            return (output);
        }

        public BitmapSource GetBitmapSource(DataProvider dataProvider)
        {
            return GetBitmapSource(dataProvider.SampleWidth, dataProvider.SampleHeight, dataProvider.SampleChannels, dataProvider.UsedPixelFormat);
        }
    }
}
