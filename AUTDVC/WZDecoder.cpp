#include "opencv2/core/core.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "Helpers.h"
using namespace cv;

namespace AUTDVC {
namespace WZDecoder {
	void SI_SimpleAverage(Mat PrevKeyQuad,Mat NextKeyQuad,Mat& SI)
	{
		static StereoBM SBM;
		static bool FirstRun = true;
		if(FirstRun)
		{
			FirstRun = false;
			SBM.init(0);
		}
		Mat disp;
		SBM(PrevKeyQuad,NextKeyQuad,SI);
				
		// Estimate the SI through a simple average
		//SI = (PrevKeyQuad + NextKeyQuad)/2.0;
		
		SI.convertTo(SI,CV_8U);
	}


	void CorrelationNoiseModeling(Mat PrevKeyQuad,Mat NextKeyQuad,vector<vector<double>>& Alphas)
	{
		//Correlation Noise Modeling
		// R = (X_F - X_B)/2
		Mat ResidualQuad = (NextKeyQuad - PrevKeyQuad)/2;

		// abs( Tn(u,v) ) = abs( DCT(Rn) ) for n-th block
		vector<Mat> ResiBlocks = AUTDVC::MapData::Quad2Blocks(ResidualQuad,Consts::DCTSize);
		for(int iBlock=0 ; iBlock<ResiBlocks.size() ; iBlock++)
		{
			Mat m;
			m = ResiBlocks[iBlock];
			m.convertTo(m,CV_64F);
			dct(m,m);
			ResiBlocks[iBlock] = abs( m );
		}

		// calculate mean and variance at each band
		vector<vector<double>> ResiCoeffs = AUTDVC::MapData::Blocks2Coeffs<double>(ResiBlocks,Consts::DCTSize);

		vector<double> BandsMean;
		vector<double> BandsVariance;
		BandsMean.resize( ResiCoeffs.size() );
		BandsVariance.resize( ResiCoeffs.size() );
		for(int iBand=0 ; iBand<ResiCoeffs.size() ; iBand++)
		{
			double e_t2 = 0.0;
			double e_t = 0.0;
			
			for(int iBlock=0 ; iBlock<ResiCoeffs[iBand].size() ; iBlock++)
			{
				e_t2 += ResiCoeffs[iBand][iBlock]*ResiCoeffs[iBand][iBlock];
				e_t += ResiCoeffs[iBand][iBlock];
			}
			e_t = e_t / ResiCoeffs[iBand].size();
			e_t2 = e_t2 / ResiCoeffs[iBand].size();
			BandsMean[iBand] = e_t;
			BandsVariance[iBand]  = e_t2 - (e_t*e_t);
		}

		// calculate distance of each coefficient from the average band value
		vector<vector<int>> DistCoeffs;
		DistCoeffs.resize( ResiCoeffs.size() );
		for(int iBand=0 ; iBand<ResiCoeffs.size() ; iBand++)
		{
			DistCoeffs[iBand].resize( ResiCoeffs[iBand].size() );
			for(int iCoeff=0 ; iCoeff<ResiCoeffs[iBand].size() ; iCoeff++)
			{
				DistCoeffs[iBand][iCoeff] = (int)pow(ResiCoeffs[iBand][iCoeff] - BandsMean[iBand],2.0);
			}
		}

		// calculate Laplacian parameter (alpha) for each coefficient
		Alphas.resize( ResiCoeffs.size() );
		for(int iBand=0 ; iBand<ResiCoeffs.size() ; iBand++)
		{
			Alphas[iBand].resize( ResiCoeffs[iBand].size() );
			for(int iCoeff=0 ; iCoeff<ResiCoeffs[iBand].size() ; iCoeff++)
			{
				if( DistCoeffs[iBand][iCoeff] < BandsVariance[iBand] )
					Alphas[iBand][iCoeff] = sqrt( 2.0/BandsVariance[iBand] );
				else
					Alphas[iBand][iCoeff] = sqrt( 2.0/DistCoeffs[iBand][iCoeff] );
			}
		}
	}


	void Reconstruction(vector<vector<int>> DecodedCoeffs,int iQuant,vector<vector<Point2d>> QuantRanges,const vector<vector<double>>& Alphas,
		const vector<vector<double>>& SideBands,const vector<vector<int>>& SideBandsQuant,vector<Mat>& RecBlocks)
	{		
		vector<vector<double>> RecCoeffs;
		RecCoeffs.resize( SideBands.size() );

		for(int iBand=0 ; iBand<RecCoeffs.size() ; ++iBand)
		{
			RecCoeffs[iBand].resize( SideBands[iBand].size() );
			for(int iCof=0 ; iCof<DecodedCoeffs[iBand].size() ; ++iCof)
			{
				double u = QuantRanges[iBand][ DecodedCoeffs[iBand][iCof] ].y;
				double l = QuantRanges[iBand][ DecodedCoeffs[iBand][iCof] ].x;
				double y = SideBands[iBand][iCof];
				double delta = u - l;
				double a = Alphas[iBand][iCof];
				
				double b = 1.0/a + delta/(1.0-exp(a*delta));
				double gamma = y - l;
				double sigma = u - y;

				if(y<l)
					RecCoeffs[iBand][iCof] = l + b;
				else if(y>=u)
					RecCoeffs[iBand][iCof] = u - b;
				else
					RecCoeffs[iBand][iCof] = y + ( (gamma+1/a)*exp(-a*gamma) - (sigma+1/a)*exp(-a*delta) ) / (2.0-exp(-a*gamma)-exp(-a*sigma));
			}
			// If this band is not WZ encoded use coefficients
			// from corresponding Side Information band
			if(DecodedCoeffs[iBand].size()==0) 
				RecCoeffs[iBand] = SideBands[iBand];
		}

		RecBlocks = MapData::Coeffs2Blocks(RecCoeffs,Consts::DCTSize);

		for(int i=0 ; i<RecBlocks.size() ; i++)
		{
			Mat m;
			RecBlocks[i].convertTo(m,CV_32F);
			cv::idct(m,RecBlocks[i]);
		}
	}
};
};