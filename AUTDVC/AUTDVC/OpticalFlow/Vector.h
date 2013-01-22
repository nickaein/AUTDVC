#pragma once

#include "stdio.h"
#include "project.h"
#include <vector>

using namespace std;

template <class T>
class VectorEx
{
protected:
	int nDim;
	T* pData;
public:
	VectorEx(void);
	VectorEx(int ndim,const T *data=NULL);
	VectorEx(const VectorEx<T>& vect);
	~VectorEx(void);
	void releaseData();
	void allocate(int ndim);
	void allocate(const VectorEx<T>& vect){allocate(vect.nDim);};
	void copyData(const VectorEx<T>& vect);
	void dimcheck(const VectorEx<T>& vect) const;
	void reset();
	double norm2() const;

	T sum() const;

	void printVector();

	// access the members
	const T* data() const{return (const T*)pData;};
	T* data() {return pData;};
	int dim() const {return nDim;};
	inline bool matchDimension(int _ndim) const {if(nDim==_ndim) return true;else return false;};
	inline bool matchDimension(const VectorEx<T>& vect) const {return matchDimension(vect.nDim);};

	// operators
	inline T operator[](int index) const {return pData[index];};
	inline T& operator[](int index){return *(pData+index);};
	VectorEx<T>& operator=(const VectorEx<T>& vect);

	//const VectorEx<T>& operator/(double val) const
	//{
	//	VectorEx<T> result(nDim);
	//	for(int i =0;i<nDim;i++)
	//		result.pData[i] = pData[i]/val;
	//	return result;
	//}
	
	VectorEx<T>& operator+=(const VectorEx<T>& vect);
	VectorEx<T>& operator*=(const VectorEx<T>& vect);
	VectorEx<T>& operator-=(const VectorEx<T>& vect);
	VectorEx<T>& operator/=(const VectorEx<T>& vect);

	VectorEx<T>& operator+=(double val);
	VectorEx<T>& operator*=(double val);
	VectorEx<T>& operator-=(double val);
	VectorEx<T>& operator/=(double val);
	
	//friend const VectorEx<T> operator+(const VectorEx<T>& vect1,const VectorEx<T>& vect2);
	//friend const VectorEx<T> operator*(const VectorEx<T>& vect1,const VectorEx<T>& vect2);
	//friend const VectorEx<T> operator-(const VectorEx<T>& vect1,const VectorEx<T>& vect2);
	//friend const VectorEx<T> operator/(const VectorEx<T>& vect1,const VectorEx<T>& vect2);

	//friend const VectorEx<T> operator+(const VectorEx<T>& vect1,double val);
	//friend const VectorEx<T> operator*(const VectorEx<T>& vect1,double val);
	//friend const VectorEx<T> operator-(const VectorEx<T>& vect1,double val);
	//friend VectorEx<T> operator/(const VectorEx<T>& vect,double val);
	
	friend double innerproduct(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
	{
		double result = 0;
		for(int i = 0;i<vect1.dim();i++)
			result += vect1[i]*vect2[i];
		return result;
	}
	
	void concatenate(const vector< VectorEx<T> >& vect);

	//friend const VectorEx<T> concatenate(const vector<VectorEx<T>>& vect){VectorEx<T> result; result.concatenate(vect); return result;};
	bool write(ofstream& myfile)
	{
		myfile.write((char *)&nDim,sizeof(int));
		myfile.write((char *)pData,sizeof(T)*nDim);
		return true;
	}
	bool read(ifstream& myfile)
	{
		myfile.read((char *)&nDim,sizeof(int));
		allocate(nDim);
		myfile.read((char *)pData,sizeof(T)*nDim);
		return true;
	}
	T mean(int N=-1) const
	{
		if(N==-1)
			N = nDim;
		T result = 0;
		for(int i = 0;i<N;i++)
			result += pData[i];
		return result/N;
	}
#ifdef _MATLAB
	void readVector(const mxArray* prhs);
	void writeVector(mxArray*& prhs) const;
#endif
};

//template <class T>
//double innerproduct(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
//{
//	double result = 0;
//	for(int i = 0;i<vect1.dim();i++)
//		result += vect1[i]*vect2[i];
//	return result;
//}

template <class T>
VectorEx<T>::VectorEx(void)
{
	nDim=0;
	pData=NULL;
}

template <class T>
VectorEx<T>::VectorEx(int ndim, const T *data)
{
	nDim=ndim;
	pData=new T[nDim];
	if(data!=NULL)
		memcpy(pData,data,sizeof(T)*nDim);
	else
		memset(pData,0,sizeof(T)*nDim);
}

template <class T>
VectorEx<T>::VectorEx(const VectorEx& vect)
{
	nDim=0;
	pData=NULL;
	copyData(vect);
}

template <class T>
VectorEx<T>::~VectorEx(void)
{
	releaseData();
}

template <class T>
void VectorEx<T>::releaseData()
{
	if(pData!=NULL)
		delete pData;
	pData=NULL;
	nDim=0;
}

template <class T>
void VectorEx<T>::allocate(int ndim)
{
	releaseData();
	nDim=ndim;
	if(nDim>0)
	{
		pData=new T[nDim];
		reset();
	}
}


template <class T>
void VectorEx<T>::copyData(const VectorEx &vect)
{
	if(nDim!=vect.nDim)
	{
		releaseData();
		nDim=vect.nDim;
		pData=new T[nDim];
	}
	memcpy(pData,vect.pData,sizeof(T)*nDim);
}

template <class T>
void VectorEx<T>::dimcheck(const VectorEx &vect) const
{
	if(nDim!=vect.nDim)
		cout<<"The dimensions of the vectors don't match!"<<endl;
}

template <class T>
void VectorEx<T>::reset()
{
	if(pData!=NULL)
		memset(pData,0,sizeof(T)*nDim);
}


template <class T>
T VectorEx<T>::sum() const
{
	T total = 0;
	for(int i=0;i<nDim;i++)
		total += pData[i];
	return total;
}

template <class T>
double VectorEx<T>::norm2() const
{
	double temp=0;
	for(int i=0;i<nDim;i++)
		temp+=pData[i]*pData[i];
	return temp;
}

template <class T>
void VectorEx<T>::printVector()
{
	for(int i=0;i<nDim;i++)
		cout<<pData[i]<<' ';
	cout<<endl;
}


//----------------------------------------------------------------------------------
// operators
//----------------------------------------------------------------------------------
template <class T>
VectorEx<T>& VectorEx<T>::operator =(const VectorEx<T> &vect)
{
	copyData(vect);
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator +=(const VectorEx<T> &vect)
{
	dimcheck(vect);
	for(int i=0;i<nDim;i++)
		pData[i]+=vect.data()[i];
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator *=(const VectorEx<T> &vect)
{
	dimcheck(vect);
	for(int i=0;i<nDim;i++)
		pData[i]*=vect.data()[i];
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator -=(const VectorEx<T> &vect)
{
	dimcheck(vect);
	for(int i=0;i<nDim;i++)
		pData[i]-=vect.data()[i];
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator /=(const VectorEx<T> &vect)
{
	dimcheck(vect);
	for(int i=0;i<nDim;i++)
		pData[i]/=vect.data()[i];
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator +=(double val)
{
	for(int i=0;i<nDim;i++)
		pData[i]+=val;
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator *=(double val)
{
	for(int i=0;i<nDim;i++)
		pData[i]*=val;
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator -=(double val)
{
	for(int i=0;i<nDim;i++)
		pData[i]-=val;
	return *this;
}

template <class T>
VectorEx<T>& VectorEx<T>::operator /=(double val)
{
	for(int i=0;i<nDim;i++)
		pData[i]/=val;
	return *this;
}


template<class T>
const VectorEx<T> operator+(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
{
	vect1.dimcheck(vect2);
	VectorEx<T> result(vect1);
	result+=vect2;
	return result;
}

template<class T>
const VectorEx<T> operator-(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
{
	vect1.dimcheck(vect2);
	VectorEx<T> result(vect1);
	result-=vect2;
	return result;
}

template<class T>
const VectorEx<T> operator*(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
{
	vect1.dimcheck(vect2);
	VectorEx<T> result(vect1);
	result*=vect2;
	return result;
}

template<class T>
const VectorEx<T> operator/(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
{
	vect1.dimcheck(vect2);
	VectorEx<T> result(vect1);
	result/=vect2;
	return result;
}

template <class T>
VectorEx<T> operator+(const VectorEx<T>& vect,double val)
{
	VectorEx<T> result(vect);
	result+=val;
	return result;
}

template <class T>
VectorEx<T> operator-(const VectorEx<T>& vect,double val)
{
	VectorEx<T> result(vect);
	result-=val;
	return result;
}

template <class T>
VectorEx<T> operator*(const VectorEx<T>& vect,double val)
{
	VectorEx<T> result(vect);
	result*=val;
	return result;
}

template <class T>
VectorEx<T> operator/(const VectorEx<T>& vect,double val)
{
	VectorEx<T> result(vect);
	result/=val;
	return result;
}


template <class T>
double innerproduct(const VectorEx<T>& vect1,const VectorEx<T>& vect2)
{
	vect1.dimcheck(vect2);
	double result=0;
	for(int i=0;i<vect1.nDim;i++)
		result+=vect1.pData[i]*vect2.pData[i];
	return result;
}

template <class T>
void VectorEx<T>::concatenate(const vector< VectorEx<T> >& vect)
{
	releaseData();
	nDim = 0;
	for(int i = 0;i<vect.size();i++)
		nDim += vect[i].dim();
	allocate(nDim);
	int dim = 0;
	for(int i = 0;i<vect.size(); i++)
	{
		for(int j = 0;j<vect[i].dim();j++)
			pData[dim+j] = vect[i][j];
		dim += vect[i].dim();
	}
}

#ifdef _QT

bool VectorEx::writeVector(QFile& file) const
{
	file.write((char *)&nDim,sizeof(int));
	if(file.write((char *)pData,sizeof(double)*nDim)!=sizeof(double)*nDim)
		return false;
	return true;
}

bool Vector::readVector(QFile &file)
{
	releaseData();
	file.read((char *)&nDim,sizeof(int));
	if(nDim<0)
		return false;
	if(nDim>0)
	{
		allocate(nDim);
		if(file.read((char *)pData,sizeof(double)*nDim)!=sizeof(double)*nDim)
			return false;
	}
	return true;
}

#endif


#ifdef _MATLAB

template <class T>
void VectorEx<T>::readVector(const mxArray* prhs)
{
	if(pData!=NULL)
		delete pData;
	int nElements = mxGetNumberOfDimensions(prhs);
	if(nElements>2)
		mexErrMsgTxt("A vector is expected to be loaded!");
	const int* dims = mxGetDimensions(prhs);
	nDim = dims[0]*dims[1];
	pData = new T[nDim];
	double* ptr = (double*)mxGetData(prhs);
	for(int i =0;i<nDim;i++)
		pData[i] = ptr[i];
}

template <class T>
void VectorEx<T>::writeVector(mxArray*& plhs) const
{
	int dims[2];
	dims[0]=nDim;dims[1]=1;
	plhs=mxCreateNumericArray(2, dims,mxDOUBLE_CLASS, mxREAL);
	double *ptr = (double*)mxGetData(plhs);
	for(int i =0;i<nDim;i++)
		ptr[i] = pData[i];
}
#endif

//*/