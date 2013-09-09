#ifndef FEM_MATVEC_H
#define FEM_MATVEC_H

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/lu.hpp> 
#include <iostream>

//using boost::numeric::ublas::prod;
//using boost::numeric::ublas::trans;
//using boost::numeric::ublas::permutation_matrix;

#include <Eigen/Dense>
//#define EIGEN_DONT_VECTORIZE 
//#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_DONT_ALIGN
// both above to avoid vectorization

//#define EIGEN_UNUSED False

#include "../macros.h"

using Eigen::Dynamic;
typedef Eigen::Matrix<double, Dynamic, Dynamic> mat;
typedef Eigen::Matrix<double, Dynamic, 1 >      vec;
typedef Eigen::Matrix<int   , Dynamic, 1 >     ivec;

//typedef Eigen::Matrix<double, 6, 6 > mat6;
//typedef Eigen::Matrix<double, 6, 1 > vec6;

inline double jacobian_det(mat const & J)
{
	if (J.rows() == J.cols()) return J.determinant();
	if (J.rows() == 2 && J.cols() == 3) 
	{
		double j1 = J(0,0)*J(1,1) - J(0,1)*J(1,0);
		double j2 = J(0,1)*J(1,2) - J(0,2)*J(1,1);
		double j3 = J(0,2)*J(1,0) - J(0,0)*J(1,2);
		return sqrt(j1*j1 + j2*j2 + j3*j3); // jacobian determinant
	}
	if (J.rows() == 1)
	{
		double sum = 0;
		for (int i=0; i<J.cols(); i++) sum += J(0,i)*J(0,i);
		return sqrt(sum);
	}
	std::cout << "jacobian_det: No rule to calculate non-square matrix jacobian" << std::endl;
	throw;
}

inline mat general_inv(mat const & M)
{
	if (M.rows() > M.cols()) return (M.transpose()*M).inverse()*M.transpose();
	if (M.rows() < M.cols()) return M.transpose()*(M*M.transpose()).inverse();
	return M.inverse();
}

inline vec make_vec(double v0, double v1)
{
	vec V(2);
	V << v0, v1;
	return V;
}

inline vec make_vec(double v0, double v1, double v2)
{
	vec V(3);
	V << v0, v1, v2;
	return V;
}

#ifdef USE_BOOST_PYTHON

namespace bp = boost::python;

vec list2vec(bp::list const & L)
{
	int size = len(L);
	vec V(size);

	for (int i=0; i<size; i++)
	{
		bp::extract<double> get_value(L[i]);
		if (get_value.check() == false) 
		{
			std::cerr << "Error: in list2vec: Invalid list of values" << std::endl;
			throw;
		}
		V(i) = get_value();
	}
	return V;
}

class py_mat: public mat
{
public:
	// Constructors
	py_mat(): mat() {};
	py_mat(mat M): mat(M) {};

	double item(int i, int j)
	{
		return this->operator()(i,j);
	}

	// Function for output
	boost::python::str __repr__()
	{
		std::ostringstream os;
		os << *this;
		return boost::python::str(os.str().c_str());
	}
};

#endif


void solve_gauss(mat & M, vec & X)
{
	int rows = M.rows();
	int cols = M.cols();
	ASSERT(rows == cols);

	// Forward elimination
	for (int i=0; i<rows; i++)
	{
		double p = M(i,i);
		ASSERT(std::abs(p) > 1.0E-8 && "solve_gauss needs pivoting");
		for (int j=i+1; j<rows; j++)
		{
			if (std::abs(M(j,i))< 1.0E-8) continue;
			double coef = M(j,i)/p;
			for (int k=i; k<cols; k++)
				M(j,k) -= coef*M(i,k);
			X(j) -= coef*X(i);
		}
	}

	// Back subtitutuion
	for (int i=rows-1; i>=0; i--)
	{
		for (int j=i+1; j<cols; j++)
		{
			X(i) -= M(i,j)*X(j);
			M(i,j) = 0.0;
		}
		X(i) /= M(i,i);
		M(i,i) = 1.0;
	}
	
}


/////////////////////////////////////////////////////////////////////////////////////////// Type definition

/*
//typedef boost::numeric::ublas::matrix<double> mat;
//typedef boost::numeric::ublas::vector<double> vec;
//typedef boost::numeric::ublas::matrix<int>    imat;
//typedef boost::numeric::ublas::vector<int>    ivec;


/////////////////////////////////////////////////////////////////////////////////////////// Matrix initiallizer

template <typename type>
class mat_setter
{
public:
	typedef typename boost::numeric::ublas::matrix<type> mat_t;

	// Constructor
	mat_setter(mat_t & M): _M(M) { }

	// Values accumulator class
	class comma_accum
	{
	public:
		comma_accum(mat_t & M, size_t rows, size_t cols, type const & first_value):
			_M(M),_rows(rows),_cols(cols),_row_idx(0), _col_idx(0) 
		{ 
			_M(0,0) = first_value;
		}
		comma_accum & operator, (type const & num) 
		{
			_col_idx++;
			if(_col_idx==_cols)
			{
				_col_idx=0;
				_row_idx++;
			}
			assert(_row_idx < _rows);
			_M(_row_idx, _col_idx) = num;
			return *this;
		}
	private:
		mat_t & _M;
		int   _rows;
		int   _cols;
		int   _row_idx;
		int   _col_idx;
	};
	
	comma_accum operator= (type const & num) 
	{
		return comma_accum(_M, _M.size1(), _M.size2(), num);
	} 

private:
	mat_t & _M;
	
};

// Helping setter function
// Usage: 
//       matrix<double> M(3,2);
//       set(M) = 1, 2, 3, 4, 5 ,6;
template <typename type>
mat_setter<type> set(boost::numeric::ublas::matrix<type> & M)
{
	return mat_setter<type>(M);
}


/////////////////////////////////////////////////////////////////////////////////////////// Vector initiallizer


template <typename type>
class vec_setter
{
public:
	typedef typename boost::numeric::ublas::vector<type> vec_t;

	// Constructor
	vec_setter(vec_t & V): _V(V) { }

	// Values accumulator class
	class comma_accum
	{
	public:
		comma_accum(vec_t & V, size_t size, type const & first_value):
			_V(V),_size(size), _idx(0)  
		{ 
			_V(0) = first_value;
		}
		comma_accum & operator, (type const & num) 
		{
			_idx++;
			assert(_idx < _size);
			_V(_idx) = num;
			return *this;
		}
	private:
		vec_t & _V;
		int   _size;
		int   _idx;
	};
	
	comma_accum operator= (type const & num) 
	{
		return comma_accum(_V, _V.size(), num);
	} 

private:
	vec_t & _V;
	
};

// Helping setter function
// Usage: 
//       vector<double> V(3);
//       set(V) = 1, 2, 3;
template <typename type>
vec_setter<type> set(boost::numeric::ublas::vector<type> & V)
{
	return vec_setter<type>(V);
}
*/

/////////////////////////////////////////////////////////////////////////////////////////// Inverse function


///** General matrix inversion routine. 
//* It uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
//* From: Thomas Lemaire (2005) */ 
//
//void inv(mat const & m, mat & inv) 
//{ 
//	if (m.size1() != m.size2()) throw "ublasExtra::lu_inv(): input matrix must be squared"; 
//	inv.resize(m.size1(), m.size1());
//	using namespace boost::numeric::ublas; 
//
//	// create a working copy of the input 
//	mat mLu(m); 
//	// perform LU-factorization 
//	lu_factorize(mLu); 
//	// create identity matrix of "inverse" 
//	inv.assign(identity_matrix<double>(m.size1())); 
//	// backsubstitute to get the inverse 
//	lu_substitute<mat const, mat >(mLu, inv); 
//} 
//
//
///////////////////////////////////////////////////////////////////////////////////////////// Determinant function
//
//
///** General matrix determinant. 
//* It uses lu_factorize in uBLAS.  */ 
//double det(mat const & m) 
//{ 
//	if (m.size1() != m.size2()) throw "ublasExtra::lu_det: matrix must be square"; 
//
//	// create a working copy of the input 
//	mat mLu(m); 
//	permutation_matrix<std::size_t> pivots(m.size1()); 
//	lu_factorize(mLu, pivots); 
//	double det = 1.0; 
//	for (size_t i=0; i < pivots.size(); ++i) 
//	{ 
//		if (pivots(i) != i) det *= -1.0; 
//		det *= mLu(i,i); 
//	} 
//	return det; 
//} 


#endif

