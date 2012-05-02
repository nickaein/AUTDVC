#include "opencv2/core/core.hpp"
#include <iostream>
#include "Helpers.h"
#include "SWCoding.h"
#include "WZDecoder.h"
#include "LDPCA.h"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;
using namespace std;

#ifdef _DEBUG
#pragma comment(lib,"opencv_calib3d231d")
#pragma comment(lib,"opencv_imgproc231d")
#else
#pragma comment(lib,"opencv_calib3d231")
#pragma comment(lib,"opencv_imgproc231")
#endif

enum CodingMode { codeJpegIntra, codeWZ, codeH264};
struct CodecSettings {
	string resourcespath;
	string outputpath;
	string videoname;
	string videofile;
	string RDlogfile;
	int iQuant;
	CodingMode codingMode;
};

void CodingWZ(int iFrame,int iQuant,const vector<vector<Mat>>& AllQs,double& framePSNRY,double& framePSNRU,double& framePSNRV,double& frameBitrate,double& frameLDPCArate)
{
	// encode and decode each quad independently
	for(int iquad=0 ; iquad<AllQs[iFrame].size() ; ++iquad)
	{
		vector<int> nBits;
		vector<vector<double*>> AccumSynd;
		vector<vector<double*>> AllBitplanes;
		vector<vector<Point2d>> QuantRanges;
		AUTDVC::EncodeSWQuad(AllQs[iFrame][iquad],iQuant,QuantRanges,AllBitplanes,AccumSynd);

		//Generate SideInformation
		Mat SideInformation;
		//Mat PrevKey = AUTDVC::Misc::EncodeDecodeJpeg(AllQs[iFrame-1][iquad],AUTDVC::Consts::JPEGQualityIndex[iQuant]);
		//Mat NextKey = AUTDVC::Misc::EncodeDecodeJpeg(AllQs[iFrame+1][iquad],AUTDVC::Consts::JPEGQualityIndex[iQuant]);
		//AUTDVC::WZDecoder::SI_SimpleAverage(PrevKey,NextKey,SideInformation);

		//////zzzzzzzzzzzzzzzzzz
		SideInformation = AUTDVC::Misc::EncodeDecodeJpeg(AllQs[iFrame][iquad], AUTDVC::Consts::JPEGQualityIndex[iQuant]);
		vector<Mat> SideBlocks;
		vector<vector<double>> SideBands;
		vector<vector<int>> SideBandsQuant;
		SideBlocks = AUTDVC::MapData::Quad2Blocks(SideInformation,AUTDVC::Consts::DCTSize);
		for(int i=0 ; i<SideBlocks.size() ; i++)
		{
			Mat m;
			m = SideBlocks[i];
			m.convertTo(m,CV_64F);
			dct(m,m);
			SideBlocks[i] = m;
		}
		SideBands = AUTDVC::MapData::Blocks2Coeffs<double>(SideBlocks,AUTDVC::Consts::DCTSize);
		SideBandsQuant = AUTDVC::Quantizer(SideBands,QuantRanges);

		vector<vector<double>> Alphas;
		AUTDVC::WZDecoder::CorrelationNoiseModeling(AllQs[iFrame-1][iquad],AllQs[iFrame+1][iquad],Alphas);

		// Decode the quad using SideInformation, accumulated syndromes in 
		// the encoder buffer, and Laplacian coefficients for CNM		
		vector<vector<int>> DecodedCoeffs;
		int nTransmittedBits;
		AUTDVC::DecodeSWQuad(SideInformation,Alphas,AllBitplanes,iQuant,QuantRanges,AccumSynd,DecodedCoeffs,SideBands,SideBandsQuant,nTransmittedBits);
		
		double nbits = 0.0;
		for(int i=0;i<16;i++)
		{
			if(AUTDVC::Consts::QLevels[iQuant][i]!=0)
				nbits += log((double)AUTDVC::Consts::QLevels[iQuant][i]) / log(2.0);
		}

		double LDPCAbitrate = nTransmittedBits / nbits / (double)AUTDVC::Consts::LDPCALength;

		// reconstruct the quad using decoded bitplanes,
		// side information and 
		vector<Mat> RecBlocks;
		AUTDVC::WZDecoder::Reconstruction(DecodedCoeffs,iQuant,QuantRanges,Alphas,SideBands,SideBandsQuant,RecBlocks);
		Mat ReconstructedQuad = AUTDVC::MapData::Blocks2Quad(RecBlocks,SideInformation.cols,SideInformation.rows);
		ReconstructedQuad.convertTo(ReconstructedQuad,CV_8U);
			
		double PSNR = AUTDVC::Misc::calcPSNR(AllQs[iFrame][iquad],ReconstructedQuad);

		// free up the memory allocated by WZ encoder
		// for the original bitplanes and syndromes
		for(int i=0;i<AllBitplanes.size();i++)
			for(int j=0;j<AllBitplanes[i].size();j++)
				delete[] AllBitplanes[i][j];
		for(int i=0;i<AccumSynd.size();i++)
			for(int j=0;j<AccumSynd[i].size();j++)
				delete[] AccumSynd[i][j];
				
		frameBitrate += nTransmittedBits;
		frameLDPCArate += LDPCAbitrate;
		if(iquad < 4) //luma
			framePSNRY += PSNR;
		else if(iquad == 4) //chroma 1
			framePSNRU += PSNR;
		else if(iquad == 5)
			framePSNRV += PSNR;
	}
	framePSNRY = framePSNRY / 3;
	frameLDPCArate = frameLDPCArate / AllQs[iFrame].size();
}

void CodingJpegIntra(Mat Y,Mat U,Mat V,int iQuant,double& framePSNRY,double& framePSNRU,double& framePSNRV,double& frameBitrate)
{
	int transbits;

	AUTDVC::Misc::calcJpegBitrate(Y,AUTDVC::Consts::JPEGQualityIndex[iQuant],framePSNRY,transbits);
	frameBitrate += transbits;

	AUTDVC::Misc::calcJpegBitrate(U,AUTDVC::Consts::JPEGQualityIndex[iQuant],framePSNRU,transbits);
	frameBitrate += transbits;

	AUTDVC::Misc::calcJpegBitrate(V,AUTDVC::Consts::JPEGQualityIndex[iQuant],framePSNRV,transbits);
	frameBitrate += transbits;
}

void Initialize(CodecSettings& cs)
{
	srand( (unsigned int)getTickCount() );
	char buf[1000];

	//generate fullpath for the video file
	cs.videofile = cs.resourcespath + cs.videoname;

	//genereate logfile name
	string logfile = cs.outputpath;
	
	//video name
	logfile += cs.videoname;
	//coding mode
	if(cs.codingMode==codeWZ)
		logfile += "_WZ";
	else if(cs.codingMode==codeJpegIntra)
		logfile += "_JPEG";
	else if(cs.codingMode==codeH264)
		logfile += "_H264Intra";
	//quantization level
	_itoa_s(cs.iQuant,buf,10);
	logfile += "_Q";
	logfile += buf;

	logfile += ".log";
	cs.RDlogfile = logfile;

	// initialize log file
	FILE* flog;
	if (fopen_s(&flog,cs.RDlogfile.c_str(),"wt")!=0)
	{
		printf("cannot open log file...");
		exit(-1);
	}
	fclose(flog);

	// initialize LDPCA encoder and decoder by loading 
	// the ladder file (to eliminate race-condition in multi-thread calls)
	string ladderfile = cs.resourcespath;
	ladderfile = ladderfile + "396_regDeg3.lad";
	LDPCAencodeBits(ladderfile.c_str(),0,0);
	LDPCAdecodeBits(ladderfile.c_str(),0,0,0,0,0,0);
}

void AppendToFile(string filename,string text)
{
	FILE* fout;
	fopen_s(&fout,filename.c_str(),"at");
	fprintf(fout,"%s",text.c_str());
	fclose(fout);
}

int main(int argc, char ** argv)
{
	CodecSettings cs;
	if(argc < 2)
	{
		cs.codingMode = codeWZ;
		printf("Enter Q parameter: ");
		cs.iQuant = 0;
		scanf("%d", &cs.iQuant); ////
		cs.resourcespath = "E:\\Thesis\\Resources\\";
		cs.outputpath = "E:\\Thesis\\Output\\";
		cs.videoname = "foreman_qcif.yuv";
	}
	
	printf("\nInitializing the Codec...\n");
	Initialize(cs);

	vector<Mat> VidY,VidU,VidV;
	// read the frames
	AUTDVC::Misc::LoadVideo(cs.videofile, Size(176,144), 1, 145, VidY, VidU, VidV);
	
	// divide the frames into the quads
	vector<vector<Mat>> AllQs;
	AllQs.resize( VidY.size() );
	for(int iFrame=0 ; iFrame<VidY.size() ; ++iFrame)
		AllQs[iFrame] = AUTDVC::Misc::convYUVtoQuads(VidY[iFrame],VidU[iFrame],VidV[iFrame],true);

	double running_time = (double)getTickCount();

	char buf[1000];

	// encode WZ Frame
	for(int iFrame=0 ; iFrame<VidY.size(); iFrame++)
	{
		double framePSNRY = 0.0;
		double framePSNRU = 0.0;
		double framePSNRV = 0.0;
		double frameBitrate = 0.0;
		double frameLDPCArate = 0.0;

		
		if(cs.codingMode==codeWZ)
		{
			if(iFrame % 2 == 1) //WZ frame
			{
				CodingWZ(iFrame,cs.iQuant,AllQs,framePSNRY,framePSNRU,framePSNRV,frameBitrate,frameLDPCArate);
			}
			else //Key-frame
			{
				CodingJpegIntra(VidY[iFrame],VidU[iFrame],VidV[iFrame],cs.iQuant,framePSNRY,framePSNRU,framePSNRV,frameBitrate);
			}
		}else if(cs.codingMode==codeJpegIntra)
		{
			CodingJpegIntra(VidY[iFrame],VidU[iFrame],VidV[iFrame],cs.iQuant,framePSNRY,framePSNRU,framePSNRV,frameBitrate);
		}
		frameBitrate = (frameBitrate / 1000.0) * AUTDVC::Consts::FramesPerSecond;

		sprintf(buf,"%05d  %05.2f %05.2f %05.2f  %07.3f %7.5f\n",iFrame,framePSNRY,framePSNRU,framePSNRV,frameBitrate,frameLDPCArate);
		printf(buf);
		AppendToFile(cs.RDlogfile,buf);
	}
	

	running_time = ((double)getTickCount() - running_time)/(double)getTickFrequency();

	printf("\r\n\r\n Running time (sec): %.2f",running_time);
	return 0;

	//cvNamedWindow("Image:", CV_WINDOW_AUTOSIZE);
	//for(int i=0;i<10;++i)
	//{
	//	Mat M;
	//	VidY[i].convertTo(M,CV_8U);

	//	stringstream s;
	//	s << "D:\\Docs\\Visual Studio 2010\\Projects\\alaki\\alaki\\out\\" << i << ".jpg";
	//	double bpp = 0,psnr = 0;
	//	AUTDVC::Misc::calcJpegBitrate(M,75,psnr,bpp);
	//}


	//cvWaitKey(0);
	//return 0;
	//CvCapture* capt = cvCreateFileCapture(fname.c_str());

	//Mat input = Mat(cvQueryFrame(capt));
 //   cv::SiftFeatureDetector detector;
 //   std::vector<cv::KeyPoint> keypoints;
 //   detector.detect(input, keypoints);

 //   // Add results to image and save.
 //   cv::Mat output;
 //   cv::drawKeypoints(input, keypoints, output);

 //   // Display the image.
 //   cvNamedWindow("Image:", CV_WINDOW_AUTOSIZE);

	//cv::imshow("Image:",output);

 //   // Wait for the user to press a key in the GUI window.
 //   cvWaitKey(0);

 //   // Free the resources.
 //   cvDestroyWindow("Image:");
 //       
 //   return 0;
}
