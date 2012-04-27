#include "opencv2/core/core.hpp"
using cv::Mat;

namespace AUTDVC {
namespace WZDecoder {
	void SI_SimpleAverage(Mat PrevKeyQuad,Mat NextKeyQuad,Mat& SI);
	void CorrelationNoiseModeling(Mat PrevKeyQuad,Mat NextKeyQuad,vector<vector<double>>& Alphas);
	void Reconstruction(vector<vector<int>> DecodedCoeffs,int iQuant,vector<vector<Point2d>> QuantRanges,const vector<vector<double>>& Alphas,
		const vector<vector<double>>& SideBands,const vector<vector<int>>& SideBandsQuant,vector<Mat>& RecBlocks);

};
};