#include "opencv2/core/core.hpp"
#include <string>
#include <vector>
using std::string;
using std::vector;
using namespace cv;

namespace AUTDVC {

/* Quantization table is from: Esmaili, "Wyner–Ziv Video Coding With Classified Correlation Noise Estimation 
		and Key Frame Coding Mode Selection", IEEE Trans. on Image Processing, 2011 */

namespace Consts {
	const int ThreadCount = 10;
	const int LDPCALength = 396;
	const int FramesPerSecond = 15;

	const short QLevels[8][16] = {
		{	16, 8, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{	32, 8, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{	32, 8, 4, 0, 8, 4, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0 },
		{	32,16, 8, 4,16, 8, 4, 0, 8, 4, 0, 0, 4, 0, 0, 0 },
		{	32,16, 8, 4,16, 8, 4, 4, 8, 4, 4, 0, 4, 4, 0, 0 },
		{	64,16, 8, 8,16, 8, 8, 4, 8, 8, 4, 4, 8, 4, 4, 0 },
		{	64,32,16, 8,32,16, 8, 4,16, 8, 4, 4, 8, 4, 4, 0 },
		{  128,64,32,16,16,32,16, 8,32,16, 8, 4,16, 8, 4, 0 }
	};
	const int JPEGQualityIndex[] = { 7, 20, 40, 60, 70, 80, 90, 95 };
	const int DCTSize = 4;
};
namespace Misc {
	int LoadVideo(string pVideoFileName, cv::Size2i Resolution, const int nStartFrame, const int nEndFrame, vector<Mat>& VidY, vector<Mat>& VidU, vector<Mat>& VidV);
	double calcPSNR(const Mat& I1, const Mat& I2);
	Mat EncodeDecodeJpeg(Mat im,int Quality);
	void calcJpegBitrate(Mat im,int Quality,double& PSNR,int& nTransmittedBits);
	vector<Mat> convYUVtoQuads(Mat Y,Mat U,Mat V, bool ChromaQuads);

	struct SSWDECODEBANDINFO {
		vector<double>* SideInfoDCT;
		vector<int>* SideInfoQuantizedBands;
		vector<int>* DecodedCoeffs;
		int nBitplanes;
		int iBand;
		vector<double*>* AccumSyndromes;
		vector<double*>* OriginalBitplanes;
		vector<double>* Alphas;
		vector<Point2d>* QuantRanges;
		int nTransmittedBits;
	};
};
namespace MapData {
	Mat Blocks2Quad(const vector<Mat>& Blocks,int QuadWidth,int QuadHeight);
	vector<Mat> Quad2Blocks(Mat Quad,int TargetBlockSize);
	
	template<class T> vector<vector<T>> Blocks2Coeffs(const vector<Mat>& Blocks,int DCTSize);
	template vector<vector<int>> Blocks2Coeffs<int>(const vector<Mat>&,int);
	template vector<vector<double>> Blocks2Coeffs<double>(const vector<Mat>&,int);

	vector<Mat> Coeffs2Blocks(const vector<vector<double>>& Coeffs,int DCTSize);

	// extract bitplanes from coefficients from LSB to MSB
	void ExtractBitplanes(const vector<int>& BandCoeffs,int nBitplanes,vector<double*>& Bitplanes);
};
namespace ImageProc {
	Mat calcDCTandQuantize(Mat Block,int iQuant);
};
};