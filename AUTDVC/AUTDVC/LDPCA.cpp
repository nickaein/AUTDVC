//LDPCA ENCODER
//Author:  David Varodayan (varodayan@stanford.edu)
//Date:    May 8, 2006

#include <stdio.h>
#include <math.h>
#include <string.h>

//returns entire accumulated syndrome of source wrt code of source rate 1
void LDPCAencodeBits(const char *InitByLadderFile, double *source, double *accumulatedSyndrome)
{
    FILE* fp = NULL;
	static int n = 0, m = 0;
    static int  *ir = NULL, *jc = NULL, *txSeq = NULL;
    int k, l;
    
	bool static FirstRun = true;

	// Read LDPCA ladder file at the 
	// first run
	if(InitByLadderFile)
	{
		int numInc = 0;
		int totalNumInc = 0;
		int numCodes = 0;
		int nzmax = 0;

		if( fopen_s(&fp,InitByLadderFile,"r") != 0)
			return;
 
		fscanf_s(fp, "%d", &numCodes);
		fscanf_s(fp, "%d", &n);
		fscanf_s(fp, "%d", &nzmax);
		fscanf_s(fp, "%d", &totalNumInc);
        
		ir = new int[nzmax];
		jc = new int[n+1];
		txSeq = new int[totalNumInc];
    
		for(k=0; k<n+1; k++)
			fscanf_s(fp, "%d", jc+k);    
		for(k=0; k<numCodes; k++)
		{
			fscanf_s(fp, "%d", &numInc);
			for(l=0; l<numInc; l++)
				fscanf_s(fp, "%d", txSeq+l);
			for(l=0; l<nzmax; l++)
				fscanf_s(fp, "%d", ir+l);
		}
	    m = (n/totalNumInc)*numInc;
		fclose(fp);

		return;
	}
    
    
    for(k=0; k<m; k++)
        accumulatedSyndrome[k] = (double) 0;
    
    //source*H'
    for(k=0; k<n; k++)
        for(l=jc[k]; l<jc[k+1]; l++)
            accumulatedSyndrome[ir[l]] += source[k];
    
    //accumulate
    for(k=1; k<m; k++)
        accumulatedSyndrome[k] += accumulatedSyndrome[k-1];
    
    //mod 2
    for(k=0; k<m; k++)
        accumulatedSyndrome[k] = (double) ((int) accumulatedSyndrome[k] % 2);

    //delete[] ir;
    //delete[] jc;
    //delete[] txSeq;
}


//LDPCA DECODER

//Author:  David Varodayan (varodayan@stanford.edu)
//Date:    May 8, 2006

#define max(a,b) (((a)>(b))?(a):(b))

int beliefPropagation(int *ir, int *jc, int m, int n, int nzmax, 
                       double *LLR_intrinsic, double *syndrome,
                       double *decoded);

//decodeBits() finds the minimum rate for which the decoded bitstream matches
//the transmitted portion of the accumulated syndrome.
//The number of residual bit errors is also calculated.
void LDPCAdecodeBits(const char* InitByLadderFile, double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors)
{
	static int numCodes = 0, n = 0, nzmax = 0, totalNumInc = 0;
	static int*  numIncAll = NULL;
	static int** txSeqAll = NULL;
	static int** irAll = NULL;
	static int *jc = NULL;

	int *ir = NULL, *txSeq = NULL;

	int m;
	int numInc;
    int code, k, currIndex, prevIndex;

	if(InitByLadderFile)
	{
		FILE *fp;
		if( fopen_s(&fp,InitByLadderFile,"r") != 0)
			return;
 
		fscanf_s(fp, "%d", &numCodes);
		fscanf_s(fp, "%d", &n);
		fscanf_s(fp, "%d", &nzmax);
		fscanf_s(fp, "%d", &totalNumInc);
        
		jc = new int[n+1];
		for(k=0; k<n+1; k++)
			fscanf_s(fp, "%d", &jc[k]);

		numIncAll = new int[numCodes];
		txSeqAll = new int*[numCodes];
		irAll = new int*[numCodes];
		for(code=0; code<numCodes; code++)
		{
			txSeqAll[code] = new int[totalNumInc]; //actual length: numInc
			irAll[code] = new int[nzmax];
			
			fscanf_s(fp, "%d", &numIncAll[code]);
			for(k=0; k<numIncAll[code]; k++)
				fscanf_s(fp, "%d", &txSeqAll[code][k]);
			for(k=0; k<nzmax; k++)
				fscanf_s(fp, "%d", &irAll[code][k]);
		}
		fclose(fp);
		return;
	}

	double* syndrome = new double[n]; //actual length: m

    //iterate through codes of increasing rate
    for(code=0; code<numCodes; code++)
    {
		numInc = numIncAll[code];
		txSeq = txSeqAll[code];
		ir = irAll[code];

		//txSeq = new int[totalNumInc];
		//ir = new int[nzmax];
		//memcpy_s(txSeq,totalNumInc*sizeof(int),txSeqAll[code],totalNumInc*sizeof(int));
		//memcpy_s(ir,nzmax*sizeof(int),irAll[code],nzmax*sizeof(int));

		m = (n/totalNumInc)*numInc;
        
        rate[0] = ((double) m)/((double) n);

        if(rate[0]==1)
        {
            for(k=0; k<n; k++)
                decoded[k]=source[k]; //result of Gaussian elimination
            numErrors[0] = 0;
			delete[] syndrome;
            return;
        }
        
        currIndex = txSeq[0];
        syndrome[0] = accumulatedSyndrome[currIndex];
        for(k=1; k<m; k++)
        {
            prevIndex = currIndex;
            currIndex = txSeq[k%numInc] + (k/numInc)*totalNumInc;
            syndrome[k] = (double) (((int) (accumulatedSyndrome[currIndex] + accumulatedSyndrome[prevIndex])) % 2);
        }

        if(beliefPropagation(ir, jc, m, n, nzmax, LLR_intrinsic, syndrome, decoded))
        {
            numErrors[0] = 0;
            for(k=0; k<n; k++)
                numErrors[0] += (double) (decoded[k]!=source[k]);
			// TODO: Remove this line and check 
			// correctness of decoded bits using a 
			// hash-based method. This below solution to
			// check the correctness is impractical since the 
			// source bits are not presen at the decoder
			//if(numErrors[0] == 0)
			{
				delete[] syndrome;
				return;
			}
        }   
    }
	delete[] syndrome;
	return;
}


//For implementation outline of beliefPropagation(), refer to 
//W. E. Ryan, "An Introduction to LDPC Codes," in CRC Handbook for Coding 
//and Signal Processing for Recording Systems (B. Vasic, ed.) CRC Press, 2004.
//available online (as of May 8, 2006) at 
//http://www.ece.arizona.edu/~ryan/New%20Folder/ryan-crc-ldpc-chap.pdf

//beliefPropagation() runs several iterations belief propagation until
//either the decoded bitstream agrees with the transmitted portion of 
//accumulated syndrome or convergence or the max number of iterations.
//Returns 1 if decoded bitstream agrees with 
//transmitted portion of accumulated syndrome.
int beliefPropagation(int *ir, int *jc, int m, int n, int nzmax, 
                       double *LLR_intrinsic, double *syndrome,
                       double *decoded)
{
    int iteration, k, l, sameCount;
    double *LLR_extrinsic, *check_LLR, *check_LLR_mag, *rowTotal, *LLR_overall;
    
    LLR_extrinsic = new double[nzmax];
    check_LLR = new double[nzmax];
    check_LLR_mag = new double[nzmax];
    rowTotal = new double[m];
    LLR_overall = new double[n];
    
    sameCount = 0;
    for(k=0; k<n; k++)
        decoded[k] = 0;
    
    //initialize variable-to-check messages
    for(k=0; k<n; k++)
        for(l=jc[k]; l<jc[k+1]; l++)
            LLR_extrinsic[l] = LLR_intrinsic[k];
    
    for(iteration=0; iteration<100; iteration++)
    {
        //Step 1: compute check-to-variable messages
        
        for(k=0; k<nzmax; k++)
        {
            check_LLR[k] = (double) ((LLR_extrinsic[k]<0) ? -1 : 1);
            check_LLR_mag[k] = ((LLR_extrinsic[k]<0) ? -LLR_extrinsic[k] : LLR_extrinsic[k]);
        }
        
        for(k=0; k<m; k++)
            rowTotal[k] = (double) ((syndrome[k]==1) ? -1 : 1);
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] *= check_LLR[k];        
        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * rowTotal[ir[k]];
            //sign of check-to-variable messages
        
        for(k=0; k<nzmax; k++)
            check_LLR_mag[k] = -log( tanh( max(check_LLR_mag[k], 0.000000001)/2 ) );
        for(k=0; k<m; k++)
            rowTotal[k] = (double) 0;
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] += check_LLR_mag[k];        
        for(k=0; k<nzmax; k++)
            check_LLR_mag[k] = -log( tanh( max(rowTotal[ir[k]] - check_LLR_mag[k], 0.000000001)/2 ) );
            //magnitude of check-to-variable messages
            
        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * check_LLR_mag[k];
            //check-to-variable messages
            
        //Step 2: compute variable-to-check messages
        
        for(k=0; k<n; k++)
        {
            LLR_overall[k] = LLR_intrinsic[k];
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_overall[k] += check_LLR[l];
        }
            
        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_extrinsic[l] = LLR_overall[k] - check_LLR[l];
                //variable-to-check messages
            
        //Step 3: test convergence and syndrome condition
        
        l = 0;
        for(k=0; k<n; k++)
            if(decoded[k] == ((LLR_overall[k]<0) ? 1 : 0))
                l++;
            else
                decoded[k] = ((LLR_overall[k]<0) ? 1 : 0);
        
        sameCount = ((l==n) ? sameCount+1 : 0); 
        
        if(sameCount==5)
		{
			delete[] LLR_extrinsic;
			delete[] check_LLR;
			delete[] check_LLR_mag;
			delete[] rowTotal;
			delete[] LLR_overall;
            return 0; //convergence (to wrong answer)
		}
        
        for(k=0; k<m; k++)
            rowTotal[k] = syndrome[k];
        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                rowTotal[ir[l]] += decoded[k];
                
        for(k=0; k<m; k++)
            if(((int) rowTotal[k] % 2) != 0)
                break;
            else if(k==m-1)
			{
				delete[] LLR_extrinsic;
				delete[] check_LLR;
				delete[] check_LLR_mag;
				delete[] rowTotal;
				delete[] LLR_overall;
                return 1; //all syndrome checks satisfied
			}
           
    }
 
    delete[] LLR_extrinsic;
    delete[] check_LLR;
    delete[] check_LLR_mag;
    delete[] rowTotal;
    delete[] LLR_overall;

    return 0;
}