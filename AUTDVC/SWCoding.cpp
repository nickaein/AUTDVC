#include <Windows.h>
#include <vector>
#include "opencv2/core/core.hpp"
#include "Helpers.h"
#include "LDPCA.h"

using namespace std;
using namespace cv;

namespace AUTDVC {
DWORD WINAPI ThreadSWDecodeBand(void* par);

// encodes the quads and returns # of bitplanes for each band and accumulated syndromes
//           Quad: quad to encode
//         iQuant: index for quantization level
//   AllBitplanes: will be filled by bitplanes for each band (Band,Bitplane,Coefficient)
//     nBitplanes: Number for allocated bitplanes for each band (calculated based on coefficient values of each band)
// AccumSyndromes: generated LDPCA syndromes for bitplanes of each band
vector<vector<int>> Quantizer(const vector<vector<double>>& DCTBands,const vector<vector<Point2d>>& Ranges) 
{
	vector<vector<int>> DCTIndicis;
	DCTIndicis.resize( DCTBands.size() );
	for(int iBand=0 ; iBand<DCTBands.size() ; iBand++)
	{
		DCTIndicis[iBand].resize( DCTBands[iBand].size() );
		for(int iCoeff=0 ; iCoeff<DCTBands[iBand].size() ; iCoeff++)
		{
			double v = DCTBands[iBand][iCoeff];
			int index = 0;
			for(; index<Ranges[iBand].size() ; index++)
			{
				if( v >= Ranges[iBand][index].x && v < Ranges[iBand][index].y )
					break;
			}
			DCTIndicis[iBand][iCoeff] = index;
		}
	}

	return DCTIndicis;
}

void EncodeSWQuad(Mat Quad,int iQuant,vector<vector<Point2d>>& QuantRanges,vector<vector<double*>>& AllBitplanes,vector<vector<double*>>& AccumSyndromes)
{
	// Divide quad into NxN blocks (N is DCT tranformation size)
	vector<Mat> Blocks = AUTDVC::MapData::Quad2Blocks(Quad,AUTDVC::Consts::DCTSize);

	// Calculate DCT of the blocks
	vector<Mat> DCTBlocks;
	DCTBlocks.resize(Blocks.size());
	for(int iblock=0 ; iblock!=Blocks.size() ; ++iblock)
	{
		Mat m = Blocks[iblock].clone();
		m.convertTo(m,CV_64F);
		dct(m,m);
		DCTBlocks[iblock] = m;
	}

	// Group the coefficients to create DCT bands
	int nBands = AUTDVC::Consts::DCTSize * AUTDVC::Consts::DCTSize;
	vector<vector<double>> DCTBands = AUTDVC::MapData::Blocks2Coeffs<double>(DCTBlocks,AUTDVC::Consts::DCTSize);

	// Calculate the quantization step size for each band
	vector<double> QStepSize;
	QStepSize.resize( nBands );
	//maxv is assumed 1024 for the DC band
	QStepSize[0] = 1024.0 / AUTDVC::Consts::QLevels[iQuant][0]; 

	for(int iBand=1 ; iBand<nBands ; iBand++)
	{
		// Find maximum amplitude for the band
		double maxv = 0;
		for(vector<double>::const_iterator v=DCTBands[iBand].begin() ; v!=DCTBands[iBand].end() ; ++v)
		{
			if( abs(*v) > maxv )
				maxv = abs(*v);
		}
		
		QStepSize[iBand] = 2 * maxv / (AUTDVC::Consts::QLevels[iQuant][iBand]);
	}

	QuantRanges.clear();
	QuantRanges.resize( nBands );

	//Uniform scalar quantizer for DC band
	for(int i=0 ; i<Consts::QLevels[iQuant][0] ; i++)
		QuantRanges[0].push_back( Point2d( i*QStepSize[0], (i+1)*QStepSize[0]) );

	// Uniform quantizer symmetric around zero for AC bands
	for(int iBand=1 ; iBand<nBands ; iBand++)
	{		
		// negative bins
		for(int i=-(int)floor(Consts::QLevels[iQuant][iBand]/2.0) ; i<-1 ; i++)
		{
			QuantRanges[iBand].push_back( Point2d(i*QStepSize[iBand], (i+1)*QStepSize[iBand]) );
		}

		// "zero" bin which has a doubled width
		QuantRanges[iBand].push_back( Point2d(-QStepSize[iBand], QStepSize[iBand]) );

		// positive bins
		for(int i=1 ; i<=(int)floor(Consts::QLevels[iQuant][iBand]/2.0) ; i++)
		{
			QuantRanges[iBand].push_back( Point2d(i*QStepSize[iBand], (i+1)*QStepSize[iBand]) );
		}
	}

	// Encode each coefficient by its index
	vector<vector<int>> DCTIndicis = Quantizer(DCTBands,QuantRanges);

	// Extract bitplanes for DCT indicis
	size_t bitplaneLength = DCTBlocks.size();
	AllBitplanes.clear();
	AllBitplanes.resize( nBands );
	for(int iBand=0 ; iBand<nBands ; iBand++)
	{
		// check if the band should be encoded 
		// and transmitted or not (be ignored)
		if(AUTDVC::Consts::QLevels[iQuant][iBand]>0)
		{
			int nbits = (int)ceil(log((double)AUTDVC::Consts::QLevels[iQuant][iBand]) / log(2.0));
			AUTDVC::MapData::ExtractBitplanes(DCTIndicis[iBand],nbits,AllBitplanes[iBand]);
		}
	}

	// LDPCA Encode each bitplane of each band
	AccumSyndromes.resize( nBands );
	for(int iBand=0 ; iBand<nBands ; iBand++)
	{
		AccumSyndromes[iBand].resize( AllBitplanes[iBand].size() );
		for(int iBitplane=0 ; iBitplane<AllBitplanes[iBand].size() ; iBitplane++)
		{
			AccumSyndromes[iBand][iBitplane] = new double[bitplaneLength];
			LDPCAencodeBits(0,AllBitplanes[iBand][iBitplane],AccumSyndromes[iBand][iBitplane]);
		}
	}
	
	return;
}

void DecodeSWQuad(Mat SideInfo,vector<vector<double>>& Alphas,vector<vector<double*>>& OriginalBitplanes,int iQuant,
	vector<vector<Point2d>>& QuantRanges,vector<vector<double*>>& AccumSyndromes,vector<vector<int>>& DecodedCoeffs,
	vector<vector<double>>& SideBands,vector<vector<int>>& SideBandsQuant,int& nTransmittedbits)
{	
	// LDPCA DECODING
	nTransmittedbits = 0;

	// Decode each band independently in a thread
	DecodedCoeffs.clear();
	DecodedCoeffs.resize( SideBandsQuant.size() );
	
	for(int iBand=0 ; iBand<SideBandsQuant.size() ;)
	{
		HANDLE* hThreads = new HANDLE[AUTDVC::Consts::ThreadCount];
		Misc::SSWDECODEBANDINFO* info = new Misc::SSWDECODEBANDINFO[AUTDVC::Consts::ThreadCount];

		int nCreatedThreads = 0;
		for(;nCreatedThreads<AUTDVC::Consts::ThreadCount && iBand<SideBandsQuant.size();)
		{
			//if the band should not be skiped
			if(AccumSyndromes[iBand].size()>0)
			{
				info[nCreatedThreads].AccumSyndromes = &AccumSyndromes[iBand];
				info[nCreatedThreads].DecodedCoeffs = &DecodedCoeffs[iBand];
				info[nCreatedThreads].iBand = iBand;
				info[nCreatedThreads].nBitplanes = (int)ceil( log((double)Consts::QLevels[iQuant][iBand])/log(2.0) );
				info[nCreatedThreads].OriginalBitplanes = &OriginalBitplanes[iBand];
				info[nCreatedThreads].Alphas = &Alphas[iBand];
				info[nCreatedThreads].SideInfoDCT = &SideBands[iBand];
				info[nCreatedThreads].SideInfoQuantizedBands = &SideBandsQuant[iBand];
				info[nCreatedThreads].QuantRanges = &QuantRanges[iBand];
				hThreads[nCreatedThreads] = CreateThread(0,0,ThreadSWDecodeBand,&info[nCreatedThreads],0,0);
				nCreatedThreads++;
			}
			iBand++;
		}

		WaitForMultipleObjects(nCreatedThreads,hThreads,TRUE,INFINITE);
		
		for(int i=0 ; i<nCreatedThreads ; i++)
		{
			nTransmittedbits += info[i].nTransmittedBits;
		}
		delete[] hThreads;
		delete[] info;
	}
}

DWORD WINAPI ThreadSWDecodeBand(void* par)
{
	Misc::SSWDECODEBANDINFO* info = (Misc::SSWDECODEBANDINFO*)par;
	vector<double>* SideInfoDCT = info->SideInfoDCT;
	vector<int>* SideInfoQuantizedBands = info->SideInfoQuantizedBands;
	vector<int>* DecodedCoeffs = info->DecodedCoeffs;
	int nBitplanes = info->nBitplanes;
	int iBand = info->iBand;
	vector<double*>* AccumSyndromes = info->AccumSyndromes;
	vector<double*>* OriginalBitplanes = info->OriginalBitplanes;
	vector<double>* Alphas = info->Alphas;
	vector<Point2d>* QuantRanges = info->QuantRanges;
	vector<double> QuantCenters;

	for(int i=0 ; i<QuantRanges->size() ; i++)
		QuantCenters.push_back( (((*QuantRanges)[i]).x + ((*QuantRanges)[i]).y)/2.0 );

	int nTransmittedbits = 0;

	size_t NCoeffs = SideInfoQuantizedBands->size();
	DecodedCoeffs->resize( NCoeffs );
	
	int maxLevel = (int)(QuantRanges->size() - 1);

	// decode starting from MSB to LSB
	for(int iBitplane=nBitplanes-1 ; iBitplane>=0 ; iBitplane--)
	{
		int currentbitMask = 1 << iBitplane;
		int decodedbitMask = (maxLevel) - ((currentbitMask<<1) - 1);
		double* p1 = new double[NCoeffs];
		double* p0 = new double[NCoeffs];

		for(unsigned int iCoeff=0 ; iCoeff<NCoeffs ; iCoeff++)
		{
			p0[iCoeff] = p1[iCoeff] = FLT_MIN;
			double sidevalue = (*SideInfoDCT)[iCoeff];

			int decodedcoeff = (*DecodedCoeffs)[iCoeff];
				
			double alpha = (*Alphas)[iCoeff];		

			for(int level=0 ; level<QuantRanges->size() ; level++)
			{
				double levelvalue = QuantCenters[level];
				// check if this level is compatible with the
				// already decoded bitplanes
				if( ((level ^ decodedcoeff) & decodedbitMask) == 0)
				{
					if( level & currentbitMask )
						p1[iCoeff] += (alpha/2) * exp( -alpha*abs(levelvalue - sidevalue)/(double)1 ) ;
					else
						p0[iCoeff] += (alpha/2) * exp( -alpha*abs(levelvalue - sidevalue)/(double)1 ) ;
				}
			}
		}

		double* pLLR = new double[NCoeffs];
		for(int iCoeff=0 ; iCoeff<NCoeffs ; iCoeff++)
			pLLR[iCoeff] = log( (p0[iCoeff])/p1[iCoeff] );
			
		delete[] p0;
		delete[] p1;

		double* decodedbits = new double[NCoeffs];
		double rate = 0;
		double numerrors = 0;
		LDPCAdecodeBits(0,pLLR,(*AccumSyndromes)[iBitplane],(*OriginalBitplanes)[iBitplane],decodedbits,&rate,&numerrors);

		nTransmittedbits += (int) (rate * AUTDVC::Consts::LDPCALength);

		delete[] pLLR;

		// TODO: maybe you need to replace "DecodedCoeffs" variable with 
		// "SideInfoQuantizedBands" or something else to update side 
		// information at each bitplane decoding.
		for(unsigned int iCoeff=0 ; iCoeff<NCoeffs ; iCoeff++)
		{
			if(decodedbits[iCoeff])
				(*DecodedCoeffs)[iCoeff] |= currentbitMask;
			else
				(*DecodedCoeffs)[iCoeff] &= (~currentbitMask);
		}
			
		delete[] decodedbits;
	}

	for(unsigned int iCoeff=0 ; iCoeff<NCoeffs ; iCoeff++)
		if( (*DecodedCoeffs)[iCoeff] >= maxLevel )
			(*DecodedCoeffs)[iCoeff] = maxLevel - 1;

	info->nTransmittedBits = nTransmittedbits;
	return 0;
}
};