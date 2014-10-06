//--------------------------------------------------------------------------
// 
//  Copyright (c) Microsoft Corporation.  All rights reserved. 
// 
//  File: ThreadSafeRandom.cs
//
//--------------------------------------------------------------------------

using System.Security.Cryptography;

namespace System.Threading
{
	/// <summary>
	/// Represents a thread-safe, pseudo-random number generator.
	/// </summary>
	public class ThreadSafeRandom : Random, IDisposable
	{
		/// <summary>Seed provider.</summary>
		private static readonly RNGCryptoServiceProvider _global = new RNGCryptoServiceProvider();
	   
		/// <summary>The underlyin provider of randomness, one instance per thread, initialized with _global.</summary>
		private ThreadLocal<Random> _local = new ThreadLocal<Random>(() =>
		{
			byte[] buffer = new byte[4];
			
			_global.GetBytes(buffer); // RNGCryptoServiceProvider is thread-safe for use in this manner
			return new Random(BitConverter.ToInt32(buffer, 0));
		});

		
		/// <summary>Returns a nonnegative random number.</summary>
		/// <returns>A 32-bit signed integer greater than or equal to zero and less than MaxValue.</returns>
		public override int Next()
		{
			return _local.Value.Next();
		}

		/// <summary>Returns a nonnegative random number less than the specified maximum.</summary>
		/// <param name="maxValue">
		/// The exclusive upper bound of the random number to be generated. maxValue must be greater than or equal to zero. 
		/// </param>
		/// <returns>
		/// A 32-bit signed integer greater than or equal to zero, and less than maxValue; 
		/// that is, the range of return values ordinarily includes zero but not maxValue. However, 
		/// if maxValue equals zero, maxValue is returned.
		/// </returns>
		public override int Next(int maxValue)
		{
			return _local.Value.Next(maxValue);
		}

		/// <summary>Returns a random number within a specified range.</summary>
		/// <param name="minValue">The inclusive lower bound of the random number returned.</param>
		/// <param name="maxValue">The exclusive upper bound of the random number returned. maxValue must be greater than or equal to minValue.</param>
		/// <returns>
		/// A 32-bit signed integer greater than or equal to minValue and less than maxValue; 
		/// that is, the range of return values includes minValue but not maxValue. 
		/// If minValue equals maxValue, minValue is returned.
		/// </returns>
		public override int Next(int minValue, int maxValue)
		{
			return _local.Value.Next(minValue, maxValue);
		}

		/// <summary>Returns a random integer percentage number within the range (0-100).</summary>
		/// <returns>
		/// A 32-bit signed integer between 0 and 100
		/// </returns>
		public int NextPercentage()
		{
			return _local.Value.Next(101);
		}

		/// <summary>Returns a random number between 0.0 and 1.0.</summary>
		/// <returns>A double-precision floating point number greater than or equal to 0.0, and less than 1.0.</returns>
		public override double NextDouble()
		{
			double t = _local.Value.NextDouble();

			while ((t < 0D) && (t >= 1D))
				t =  _local.Value.NextDouble();
			
			return t;
		}

		public double NextDouble(double stdDev, double mean)
		{
			double t = _local.Value.NextDouble();

			while ((t < 0D) && (t >= 1D))
				t = _local.Value.NextDouble();

			t *= stdDev;
			t += mean;

			return t;
		}

		public double NextDouble(double stdDev)
		{
			double t = _local.Value.NextDouble();

			while ((t < 0D) && (t >= 1D))
				t = _local.Value.NextDouble();
	   
			t *= 2D;
			t -= 1D;
			t *= stdDev;
									
			return t;
		}

		/// <summary>Fills the elements of a specified array of bytes with random numbers.</summary>
		/// <param name="buffer">An array of bytes to contain random numbers.</param>
		public override void NextBytes(byte[] buffer)
		{
			_local.Value.NextBytes(buffer);
		}

		#region IDisposable Members
		~ThreadSafeRandom()
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
				if (_local != null)
					_local.Dispose();
			}
			// free native resources
		}

		public void Dispose()
		{
			Dispose(true);
			GC.SuppressFinalize(this);
		}
	   
		#endregion
	}
}