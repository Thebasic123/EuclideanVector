#ifndef EUCLIDEANVECTOR_HPP
#define EUCLIDEANVECTOR_HPP


#include <iostream>
#include <sstream>
#include <cmath>
#include <initializer_list>
#include <vector>
#include <list>



namespace evec{
	class EuclideanVector{
	private:
		unsigned int dim_;
		double *mag_;

	public:
		/*  Constructors
		 *  
		 * */
		EuclideanVector(const unsigned int d);
		EuclideanVector(const unsigned int d, const double m);
		EuclideanVector(std::vector<double>::iterator begin, std::vector<double>::iterator end);
		EuclideanVector(std::list<double>::iterator begin, std::list<double>::iterator end);
		EuclideanVector(std::initializer_list<double> arr);
		EuclideanVector(const EuclideanVector& ev);
		EuclideanVector(EuclideanVector&& ev);
		~EuclideanVector();

		/*  Operations
		 *  
		 * */
		EuclideanVector& operator=(const EuclideanVector& ev);
		EuclideanVector& operator=(EuclideanVector&& ev);
		unsigned int getNumDimensions() const;
		double get(unsigned int d) const;
		double getEuclideanNorm() const;
		EuclideanVector& createUnitVector() const;

		double& operator[] (int index);
		double operator[] (int index) const;

		EuclideanVector& operator+=(const EuclideanVector &ev);
		EuclideanVector& operator-=(const EuclideanVector &ev);
		EuclideanVector& operator*=(const double s);
		EuclideanVector& operator/=(const double s);
	
		operator std::vector<double>() const;
		operator std::list<double>() const;

		/*  Friends
		 *  
		 * */
		friend bool operator==(const EuclideanVector &lhs, const EuclideanVector &rhs);
		friend bool operator!=(const EuclideanVector &lhs, const EuclideanVector &rhs);
		friend EuclideanVector& operator+(const EuclideanVector &lhs, const EuclideanVector &rhs);
		friend EuclideanVector& operator-(const EuclideanVector &lhs, const EuclideanVector &rhs);
		friend double operator*(const EuclideanVector &lhs, const EuclideanVector &rhs);
		friend EuclideanVector& operator*(const EuclideanVector &lhs, double rhs);
		friend EuclideanVector& operator/(const EuclideanVector &lhs, double rhs);
		friend std::ostream& operator<<(std::ostream &os, const EuclideanVector &v);

	};
}


#endif