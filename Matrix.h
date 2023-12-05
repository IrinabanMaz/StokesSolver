#ifndef MATRIX
#define MATRIX

#include<array>
#include<iostream>

template<size_t N>
class vec : public std::array<double, N>
{
public:
	
	vec<N>(double x = 0.0)
	{
		for (int i = 0; i < N; i++)
			(*this)[i] = x;
	}


	vec<N> operator =(double c);
	//element wise operations.
	vec<N> operator +(double c);
	vec<N> operator +(vec<N> v);
	friend vec<N> operator +(double c, vec<N> v);

	vec<N> operator -(double c);
	vec<N> operator -(vec<N> v);
	friend vec<N> operator -(double c, vec<N> v);

	vec<N> operator *(double c);
	vec<N> operator *(vec<N> v);
	friend vec<N> operator *(double c, vec<N> v);

	vec<N> operator /(double c);
	vec<N> operator /(vec<N> v);

};

template<size_t N>
vec<N> vec<N>::operator =(double c)
{
	for (int i = 0; i < N; i++)
		(*this)[i] = c;

	return *this;
}

template<size_t N>
vec<N> vec<N>::operator +(double c)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] + c;
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator +(vec<N> v)
{

	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] + v[i];
	return temp;
}

template<size_t N>
vec<N> operator +(double c, vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = v[i] + c;
	return temp;
}

template<size_t N>
vec<N> operator -(double c, vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = c - v[i];
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator -(double c)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] - c;
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator -(vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] - v[i];
	return temp;
}

template<size_t N>
vec<N> operator *(double c, vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = v[i] * c;
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator *(double c)
{
	vec<N> temp;
	for (int i = 0; i <  N; i++)
		temp[i] = (*this)[i] * c;
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator *(vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] * v[i];
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator /(vec<N> v)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] / v[i];
	return temp;
}

template<size_t N>
vec<N> vec<N>::operator /(double c)
{
	vec<N> temp;
	for (int i = 0; i < N; i++)
		temp[i] = (*this)[i] / c;
	return temp;
}

template<size_t N>
double sum(vec<N> v)
{
	double temp = 0;
	for (int i = 0; i < N; i++)
		temp += v[i];
	return temp;
}

template<size_t N>
double dot(vec<N> u, vec<N> v)
{
	double temp = 0;
	for (int i = 0; i < N; i++)
		temp += u[i] * v[i];
	return temp;
}

template<size_t N>
double maximum(vec<N> v)
{
	double temp = v[0];
	for (int i = 1; i < N; i++)
		if (v[i] > temp)
			temp = v[i];
	return temp;
}

template<size_t N>
size_t argmax(vec<N> v)
{
	size_t temp = 0;
	double max = v[0];
	for(int i = 1; i < N; i++)
		if (v[i] > max)
		{
			max = v[i];
			temp = i;
		}
	return temp;
}

inline double toDouble(vec<1> v)
{
	return v[0];
}

template<size_t N>
std::ostream& operator <<(std::ostream& o, vec<N> v)
{
	o << "[ ";
	for (int i = 0; i < N - 1; i++)
		o << v[i] << ", ";
	o << v[N - 1] << " ]";
	return o;
}

template<size_t N2, size_t N1 >
vec<N2> subvector(vec<N1> A, size_t pos)
{

	if ( pos + N2 > N1)
		std::cout << "Warning, submatrix out of bounds." << std::endl;

	vec<N2> sub;

		for (int i = pos; i < pos + N2; i++)
			sub[i - pos] = A[i];

	return sub;
}


template<size_t M,size_t N>
class Matrix : public std::array<vec<N>, M>
{
public:

	Matrix<M, N>(vec<N> v = 0.0)
	{
		for (int i = 0; i < M; i++)
			(*this)[i] = v;
	}


	Matrix<M,N> operator +(double c);
	friend Matrix<M, N>  operator +(double c, Matrix<M, N> A);
	Matrix<M, N> operator +(vec<N> v);
	friend Matrix<M, N>  operator +(vec<N> v, Matrix<M, N> A);
	Matrix<M, N>  operator +(Matrix<M, N> A);

	Matrix<M, N> operator -(double c);
	friend Matrix<M, N> operator -(double c, Matrix<M, N> A);
	Matrix<M, N> operator -(vec<N> v);
    friend Matrix<M, N>  operator -(vec<N> v, Matrix<M, N> A);
	Matrix<M, N> operator -(Matrix<M, N> A);

	Matrix<M, N> operator *(double c);
	friend Matrix<M, N> operator *(double c, Matrix<M, N> A);
	Matrix<M, N> operator *(vec<N> v);
	friend Matrix<M, N>  operator *(vec<N> v, Matrix<M, N> A);
	Matrix<M, N> operator *(Matrix<M, N> A);

	Matrix<M, N> operator /(double c);
	Matrix<M, N> operator /(vec<N> v);
	Matrix<M, N> operator /(Matrix<M, N> A);

	Matrix<M, N> operator =(double c);

	Matrix<N, M> T();

};

template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator +(double c)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] + c;

	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator +(double c, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = A[i] + c;
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator +(vec<N> v)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] + v;
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator +(vec<N> v, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = v + A[i];
}
template<size_t M, size_t N>
Matrix<M, N> Matrix < M, N>::operator +(Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] + A[i];
	return temp;
}

template<size_t M, size_t N>
Matrix<M,N> Matrix<M,N>:: operator -(double c)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] - c;

	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator -(double c, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = c - A[i];
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator -(vec<N> v)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] - v;
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator -(vec<N> v, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = v - A[i];
}
template<size_t M, size_t N>
Matrix<M, N> Matrix < M, N>::operator -(Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] - A[i];
	return temp;
}

template<size_t M, size_t N>
Matrix<M,N> Matrix<M,N>::operator *(double c)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] * c;

	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator *(double c, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = A[i] * c;
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator *(vec<N> v)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] * v;
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> operator *(vec<N> v, Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = v * A[i];
}
template<size_t M, size_t N>
Matrix<M, N> Matrix < M, N>::operator *(Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] * A[i];
	return temp;
}

template<size_t M, size_t N>
Matrix<M,N> Matrix<M,N>::operator /(double c)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] / c;

	return temp;
}

template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator /(vec<N> v)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] / v;
	return temp;
}

template<size_t M, size_t N>
Matrix<M, N> Matrix < M, N>::operator /(Matrix<M, N> A)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		temp[i] = (*this)[i] / A[i];
	return temp;
}
template<size_t M, size_t N>
Matrix<M, N> Matrix<M, N>::operator =(double c)
{
	for (int i = 0; i < M; i++)
		(*this)[i] = c;

	return *this;
}


template<size_t M, size_t N>
Matrix<N, M> Matrix<M, N>::T()
{
	Matrix<N, M> temp;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			temp[j][i] = (*this)[i][j];

	return temp;
}

template<size_t M, size_t N, size_t P>
Matrix<M, P> dot(Matrix<M, N> A, Matrix<N, P> B)
{
	Matrix<P,N> BT = B.T();
	Matrix<M, P> temp;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < P; j++)
			temp[i][j] = dot(A[i], BT[j]);

	return temp;
}



template<size_t M, size_t N, size_t P>
Matrix<M, P> outer(Matrix<M, N> A, Matrix<P,N> BT)
{
	Matrix<M, P> temp;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < P; j++)
			temp[i][j] = dot(A[i], BT[j]);

	return temp;
}
template<size_t M, size_t N, size_t P>
Matrix<M, N> inner(Matrix<P, M> A, Matrix<P, N> B)
{
	return dot(A.T(), B);
}
template<size_t M, size_t N>
vec<M> dot(Matrix<M, N> A, vec<N> v)
{
	vec<M> temp;

	for (int i = 0; i < M; i++)
		temp[i] = dot(A[i], v);

	return temp;
}

template<size_t M, size_t N>
vec<M> colsum(Matrix<M, N> A)
{
	vec<M> temp;
	for (int i = 0; i < M; i++)
		temp[i] = sum(A[i]);
	return temp;
}

template<size_t M, size_t N>
vec<N> rowsum(Matrix<M, N> A)
{
	vec<N> temp;
	temp = 0.0;
	for (int j = 0; j < N; j++)
		for (int i = 0; i < M; i++)
			temp[j] += A[i][j];

	return temp;

}
template<size_t M, size_t N>
double sum(Matrix<M, N> A)
{
	return sum(colsum(A));
}
template<size_t M,size_t N>
vec<M> colsqsum(Matrix<M, N> A)
{
	vec<M> temp;
	for (int i = 0; i < M; i++)
		temp[i] = dot(A[i], A[i]);
	return temp;
}

template<size_t M, size_t N>
vec<N> rowsqsum(Matrix<M, N> A)
{
	vec<N> temp = 0;
	for (int j = 0; j < N; j++)
		for (int i = 0; i < M; i++)
			temp[j] += A[i][j] * A[i][j];

	return temp;
}

template<size_t M,size_t N>
double sqsum(Matrix<M, N> A)
{
	return sum(colsqsum(A));
}

template<size_t M,size_t N>
vec<M> colmax(Matrix<M, N> A)
{
	vec<M> temp;
	for (int i = 0; i < M; i++)
		temp[i] = maximum(A[i]);
	return temp;
}

template<size_t M, size_t N>
vec<N> rowmax(Matrix<M, N> A)
{
	vec<N> temp;
	double max = 0;
	for (int j = 0; j < N; j++)
	{
		max = A[0][j];
		for (int i = 1; i < M; i++)
			if (A[i][j] > max)
				max = A[i][j];
					
		temp[j] = max;
	}
}

template<size_t M,size_t N>
double maximum(Matrix<M, N> A)
{
	return maximum(colmax(A));
}

template<size_t M,size_t N>
std::ostream& operator <<(std::ostream& o, Matrix<M, N> A)
{
	for (int i = 0; i < M; i++)
		o << A[i] << "\n";

	return o;
}

template<size_t M, size_t N>
Matrix<M, N> outer(vec<M> u, vec<N> v)
{
	Matrix<M, N> temp;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			temp[i][j] = u[i] * v[j];

	return temp;
}

template<size_t M2 , size_t N2 , size_t M1 , size_t N1 > 
Matrix<M2, N2> submatrix(Matrix<M1, N1> A, size_t xpos, size_t ypos)
{

	if ((xpos + M2 > M1) || (ypos + N2 > N1))
	{
		std::cout << "Warning, submatrix out of bounds. " << std::endl;
		std::cout << "Matrix dims: " << M1 << " x " << N1 << std::endl;
		std::cout << "subMatrix dims: " << M2 << " x " << N2 << std::endl;
		std::cout << "at coordinate: (" << xpos << ", " << ypos << ")" << std::endl;
	}

	Matrix<M2, N2> sub;

	for (int i = xpos; i < xpos + M2; i++)
		for (int j = ypos; j < ypos + N2; j++)
			sub[i - xpos][j - ypos] = A[i][j];

	return sub;
}

#endif MATRIX