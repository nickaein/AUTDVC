#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "Helpers.h"
#include <string>
#include <vector>

using std::string;
using std::vector;
using namespace cv;

namespace AUTDVC {
namespace Misc {
int LoadVideo(string pVideoFileName, cv::Size2i Resolution, const int nStartFrame, const int nEndFrame, vector<Mat>& VidY, vector<Mat>& VidU, vector<Mat>& VidV)
{
	bool bSuccess = true;

	// Determine dimensions
	int nWidth = Resolution.width;
	int nHeight = Resolution.height;
	// Allocate space
	int nFrames = nEndFrame - nStartFrame + 1;
	if (nFrames <= 0) {
		fprintf(stderr, "Invalid frame bounds: %d %d\n", nStartFrame, nEndFrame);
		return -1;
	}
	
	VidY.clear();
	VidU.clear();
	VidV.clear();

	FILE* pFile = fopen(pVideoFileName.c_str(), "rb");
	if (pFile == NULL) {
		fprintf(stderr, "Video file cannot be opened: %s\n", pVideoFileName);
		return -1;
	}

	int nPixelsPerChannel = nWidth*nHeight;
	uchar* pBuffer = new uchar[nPixelsPerChannel];
	
	// Skip frames before start of subsequence
	for (int nFrame = 1; nFrame < nStartFrame; nFrame++) {
		// Read Y frame
		size_t nRead = fread(pBuffer, 1, nPixelsPerChannel, pFile);
		// Read U frame
		nRead = fread(pBuffer, 1, nPixelsPerChannel/4, pFile);
		// Read V frame
		nRead = fread(pBuffer, 1, nPixelsPerChannel/4, pFile);
	} // end nFrame

	// Process frames of subsequence
	for (int nFrame = nStartFrame; nFrame <= nEndFrame; nFrame++) {
		size_t nRead;
		Mat m;
		
		// Read Y frame
		nRead = fread(pBuffer, 1, nPixelsPerChannel, pFile);
		m = Mat(nHeight,nWidth,CV_8U,pBuffer);
		VidY.push_back( m.clone() );

		// Read U frame
		nRead = fread(pBuffer, 1, nPixelsPerChannel/2/2, pFile);		
		m = Mat(nHeight/2,nWidth/2,CV_8U,pBuffer);
		VidU.push_back( m.clone() );

		// Read V frame
		nRead = fread(pBuffer, 1, nPixelsPerChannel/2/2, pFile);
		m = Mat(nHeight/2,nWidth/2,CV_8U,pBuffer);
		VidV.push_back( m.clone() );
	}
	fclose(pFile);

	delete[] pBuffer;

	return bSuccess ? 0 : -1;
}

double calcPSNR(const Mat& I1, const Mat& I2)
{
 Mat s1;
 absdiff(I1, I2, s1);       // |I1 - I2|
 s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
 s1 = s1.mul(s1);           // |I1 - I2|^2

 Scalar s = sum(s1);         // sum elements per channel

 double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels

 if( sse <= 1e-20) // for small values return zero
     return 0;
 else
 {
     double  mse = sse /(double)(I1.channels() * I1.total());
     double psnr = 10.0*log10((255*255)/mse);
     return psnr;
 }
}

Mat EncodeDecodeJpeg(Mat im,int Quality)
{
	//Set Quality
	vector<int> params;
	params.push_back(CV_IMWRITE_JPEG_QUALITY);
	params.push_back(Quality);
	
	//Encode the image into JPEG
	vector<uchar> buf_encodedjpeg;
	cv::imencode(".jpg",im,buf_encodedjpeg,params);
	Mat im2 = cv::imdecode(buf_encodedjpeg,0);
	return im2;
}

void calcJpegBitrate(Mat im,int Quality,double& PSNR,int& nTransmittedBits)
{
	//Set Quality
	vector<int> params;
	params.push_back(CV_IMWRITE_JPEG_QUALITY);
	params.push_back(Quality);
	
	//Encode the image into JPEG
	vector<uchar> buf_encodedjpeg;
	cv::imencode(".jpg",im,buf_encodedjpeg,params);
	Mat im2 = cv::imdecode(buf_encodedjpeg,0);
	
	nTransmittedBits = (int) (8 * (buf_encodedjpeg.size()));
	PSNR = calcPSNR(im,im2);
}


vector<Mat> convYUVtoQuads(Mat Y,Mat U,Mat V, bool ChromaQuads)
{
	int w = U.cols;
	int h = U.rows;
	
	Rect R1(0,0,w,h);
	Rect R2(w,0,w,h);
	Rect R3(0,h,w,h);
	Rect R4(w,h,w,h);

	vector<Mat> Quads;

	Quads.push_back( Y(R1).clone() );
	Quads.push_back( Y(R2).clone() );
	Quads.push_back( Y(R3).clone() );
	Quads.push_back( Y(R4).clone() );

	if(ChromaQuads)
	{
		Quads.push_back( U.clone() );
		Quads.push_back( V.clone() );
	}

	return Quads;
}
};

namespace MapData {
	Mat Blocks2Quad(const vector<Mat>& Blocks,int QuadWidth,int QuadHeight)
	{
		int blockWidth = Blocks[0].cols;
		int blockHeight = Blocks[0].rows;
		int nBlocksInRow = QuadWidth / blockWidth;
		int nBlocksInCol = QuadHeight / blockHeight;
		Mat Quad(QuadHeight,QuadWidth,Blocks[0].type());

		for(int i=0 ; i<nBlocksInRow ; i++)
			for(int j=0 ; j<nBlocksInCol ; j++)
			{
				cv::Rect rct(i*blockWidth,j*blockHeight,blockWidth,blockHeight);
				Blocks[i*nBlocksInCol + j].copyTo( Quad(rct) );
			}
		return Quad;
	}

	vector<Mat> Quad2Blocks(Mat Quad,int TargetBlockSize)
	{
		vector<Mat> Splitted;
		int w = Quad.cols;
		int h = Quad.rows;
		int nBlocks = (w*h) / (TargetBlockSize*TargetBlockSize);
		for(int i=0;i<w;i+=TargetBlockSize)
		{
			for(int j=0;j<h;j+=TargetBlockSize)
			{
				Rect BlockRect(i,j,TargetBlockSize,TargetBlockSize);
				Splitted.push_back( Quad(BlockRect) );
			}
		}

		return Splitted;
	}

	template<class T> vector<vector<T>> Blocks2Coeffs(const vector<Mat>& Blocks,int DCTSize)
	{
		vector<vector<T>> Coeffs;

		int nBands = DCTSize*DCTSize;
		Coeffs.resize( nBands );

		for(int iBand=0 ; iBand<nBands ; ++iBand)
		{
			int bandx = iBand / DCTSize;
			int bandy = iBand % DCTSize;
			for(int iblock=0 ; iblock<Blocks.size() ; iblock++)
			{
				Coeffs[iBand].push_back( Blocks[iblock].at<T>(bandx,bandy) );
			}
		}
		return Coeffs;
	}


	vector<Mat> Coeffs2Blocks(const vector<vector<double>>& Coeffs,int DCTSize)
	{
		vector<Mat> Blocks;

		for(int iblock=0 ; iblock<Coeffs[0].size() ; iblock++)
		{
			Mat blk(DCTSize,DCTSize,CV_64F);
			
			for(int iband=0 ; iband<Coeffs.size(); iband++)
			{
				int bandx = iband / DCTSize;
				int bandy = iband % DCTSize;
				blk.at<double>(bandx,bandy) = Coeffs[iband][iblock];
			}
			Blocks.push_back( blk );
		}

		return Blocks;
	}

	void ExtractBitplanes(const vector<int>& BandCoeffs,int nBitplanes,vector<double*>& Bitplanes)
	{
		Bitplanes.clear();
		Bitplanes.resize( nBitplanes );
		for(int iBitplane=0 ; iBitplane<nBitplanes ; iBitplane++)
		{
			Bitplanes[iBitplane] = new double[BandCoeffs.size()];

			int bitMask = 1<<iBitplane;

			for(int i=0;i<BandCoeffs.size();i++)
			{
				Bitplanes[iBitplane][i] = (double)((BandCoeffs[i] & bitMask)!=0);
			}
		}
	}

};
};