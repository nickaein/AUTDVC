#include "opencv2/core/core.hpp"
#include "opencv2/video/tracking.hpp"
#include "Helpers.h"
#include "opencv2/highgui/highgui.hpp"

#include "OpticalFlow/OpticalFlow.h"
#include <fstream>
using namespace std;

using namespace cv;

namespace AUTDVC {
namespace WZDecoder {
	void SI_OpticalFlow1(Mat PrevKeyQuad,Mat NextKeyQuad,Mat& SI,Mat& Residual)
	{
		Mat PK = PrevKeyQuad.clone();
		Mat NK = NextKeyQuad.clone();

		PK.convertTo(PK,CV_64F);
		NK.convertTo(NK,CV_64F);
		PK = PK / 255.0;
		NK = NK / 255.0;


		DImage Im1(PK.size[1],PK.size[0]);
		DImage Im2(NK.size[1],NK.size[0]);
		memcpy(Im1.pData,PK.data,sizeof(double)*Im1.npixels());
		memcpy(Im2.pData,NK.data,sizeof(double)*Im2.npixels());

		double alpha=0.012;
		double ratio=0.8;
		int minWidth= 5;
		int nOuterFPIterations = 7;
		int nInnerFPIterations = 3;
		int nSORIterations= 30;
		OpticalFlow::IsDisplay = false;
		DImage vx1,vy1,warpI2;
		DImage vx2,vy2;
		OpticalFlow::Coarse2FineFlow(vx1,vy1,warpI2,Im1,Im2,alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
		OpticalFlow::Coarse2FineFlow(vx2,vy2,warpI2,Im2,Im1,alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
		
		//DImage tmp;
		//// map frame (t-1) to fram (t+1) 
		//// using calculated optical flow
		//ReconstructByOF(tmp,Im1,vx1,vy1);
		//PrevToNext = Mat(PrevKeyQuad.size(),CV_64F);
		//memcpy(PrevToNext.data,tmp.pData,sizeof(double)*tmp.npixels());
		//PrevToNext.convertTo(PrevToNext,CV_8U);

		//// map frame (t+1) to fram (t-1)
		//// using calculated optical flow
		//ReconstructByOF(tmp,Im2,vx2,vy2);
		//NextToPrev = Mat(NextKeyQuad.size(),CV_64F);
		//memcpy(NextToPrev.data,tmp.pData,sizeof(double)*tmp.npixels());
		//NextToPrev.convertTo(NextToPrev,CV_8U);

		vx1.Multiplywith(0.5);
		vy1.Multiplywith(0.5);
		vx2.Multiplywith(0.5);
		vy2.Multiplywith(0.5);
		//OpticalFlow::showFlow(vx,"fl1.bmp");
		//OpticalFlow::showFlow(vy,"fl2.bmp");

		DImage Out1,Out2;
		OpticalFlow::warpFL(Out1,Im1,Im2,vx1,vy1);
		//vx1.Multiplywith(-1.0);
		//vy1.Multiplywith(-1.0);
		//OpticalFlow::warpFL(Out2,Im2,Im1,vx1,vy1);
		OpticalFlow::warpFL(Out2,Im2,Im1,vx2,vy2);

		memcpy(PK.data,Out1.pData,sizeof(double)*Out1.npixels());
		memcpy(NK.data,Out2.pData,sizeof(double)*Out2.npixels());

		//Generate side information by averaging f_1 and f_2
		SI = ((PK + NK)/2.0) * 255.0;
		SI.convertTo(SI,CV_8U);

		//imwrite("F20.bmp",SI);


		Residual = 255.0 * (NK - PK)/1.0;
	}


	void CorrelationNoiseModeling(Mat Residual,vector<vector<double>>& Alphas)
	{
		// abs( Tn(u,v) ) = abs( DCT(Rn) ) for n-th block
		vector<Mat> ResiBlocks = AUTDVC::MapData::Quad2Blocks(Residual,Consts::DCTSize);
		for(int iBlock=0 ; iBlock<ResiBlocks.size() ; iBlock++)
		{
			Mat m;
			m = ResiBlocks[iBlock];
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

	void SI_SimpleAverage(Mat PrevKeyQuad,Mat NextKeyQuad,Mat& SI)
	{
		//static StereoBM SBM;
		//static bool FirstRun = true;
		//if(FirstRun)
		//{
		//	FirstRun = false;
		//	SBM.init(0);
		//}
		//Mat disp;
		//SBM(PrevKeyQuad,NextKeyQuad,SI);		

		// Estimate the SI through a simple average
		SI = (PrevKeyQuad + NextKeyQuad)/2.0;
		
		SI.convertTo(SI,CV_8U);
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