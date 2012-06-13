#include <vector>
#include "opencv2/core/core.hpp"

using std::vector;
using cv::Mat;
namespace AUTDVC {
	vector<vector<int>> Quantizer(const vector<vector<double>>& DCTBands,const vector<vector<Point2d>>& Ranges);
	void EncodeSWQuad(Mat Quad,int iQuant,vector<vector<Point2d>>& QuantRanges,vector<vector<double*>>& AllBitplanes,vector<vector<double*>>& AccumSyndromes);

	void DecodeSWQuad(Mat SideInfo,vector<vector<double>>& Alphas,vector<vector<double*>>& OriginalBitplanes,int iQuant,
		vector<vector<Point2d>>& QuantRanges,vector<vector<double*>>& AccumSyndromes,vector<vector<int>>& DecodedCoeffs,
		vector<vector<double>>& SideBands,vector<vector<int>>& SideBandsQuant,int& nTransmittedbits);
};