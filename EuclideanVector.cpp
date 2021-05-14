#include "EuclideanVector.h"


namespace evec{
		/*  Constructors
		 *  
		 * */

		/** Construct euclidean vector from the number of dimensions, set all magnitude to 0
	    * 
	    @param d the number of dimensions
	    @return
		*/
		EuclideanVector::EuclideanVector(const unsigned int d)
			: EuclideanVector(d, 0.0) 
		{}

		/** Construct euclidean vector from the number of dimensions and initialize the magnitude
			for all dimension as the same value
	    * 
	    @param d the number of dimensions
	    	   m the magnitude for each dimension
	    @return
		*/
		EuclideanVector::EuclideanVector(const unsigned int d, const double m)
			: dim_(d)
			, mag_(new double[dim_])
		{
			for(unsigned int i=0; i<dim_; ++i)
				mag_[i] = m;
		}

		/** Construct euclidean vector from a vector
	    * 
	    @param begin the iterator pointing to the beginning of the vector
	    	   end   the iterator pointing to the end of the vectir
	    @return
		*/
		EuclideanVector::EuclideanVector(std::vector<double>::iterator begin, std::vector<double>::iterator end)
			: dim_(std::distance(begin, end)) 
			, mag_(new double[dim_])
		{
			for(std::vector<double>::iterator it= begin; it!=end;++it)
				mag_[std::distance(begin, it)] =  *it;
		}

		/** Construct euclidean vector from a list
	    * 
	    @param begin the iterator pointing to the beginning of the list
	    	   end   the iterator pointing to the end of the list
	    @return
		*/
		EuclideanVector::EuclideanVector(std::list<double>::iterator begin, std::list<double>::iterator end)
			: dim_(std::distance(begin, end)) 
			, mag_(new double[dim_])
		{
			for(std::list<double>::iterator it= begin; it!=end;++it)
				mag_[std::distance(begin, it)] =  *it;

		}

		/** Construct euclidean vector from a list of double values
	    * 
	    @param arr the list of double values
	    @return
		*/
		EuclideanVector::EuclideanVector(std::initializer_list<double> arr)
			: dim_(arr.end()-arr.begin())
			, mag_(new double[dim_])
		{
			for(auto it= arr.begin(); it!=arr.end();++it)
				mag_[it-arr.begin()] = *it;
		}

		/** Copy Constructor
	    * 
	    @param 
	    @return
		*/
		EuclideanVector::EuclideanVector(const EuclideanVector& ev)
			: dim_(ev.dim_)
			, mag_(new double[dim_])
		{
			for(unsigned int i = 0; i<dim_; ++i)
				mag_[i] = ev.mag_[i];
		}

		/** Move Constructor
	    * 
	    @param 
	    @return
		*/
		EuclideanVector::EuclideanVector(EuclideanVector&& ev)
			: dim_(std::move(ev.dim_))
			, mag_(std::move(ev.mag_))
		{
			ev.mag_ = nullptr;
			ev.dim_ = 0;
		}
		
		/** Destructor
	    * 
	    @param 
	    @return
		*/
		EuclideanVector::~EuclideanVector(){
			delete [] mag_;
		}
		

		/*  Operations
		 *  
		 * */

		/** Copy Assignment
	    * 
	    @param
	    @return
		*/
		EuclideanVector& EuclideanVector::operator=(const EuclideanVector& ev){
			if (this != &ev) {
				dim_ = ev.dim_;
				mag_ = new double[dim_];
				for(unsigned int i = 0; i<dim_; ++i)
					mag_[i] = ev.mag_[i];
			}
			return *this;
		}

		/** Move assignment
	    * 
	    @param 
	    @return
		*/
		EuclideanVector& EuclideanVector::operator=(EuclideanVector&& ev){
			if (this != &ev) {
				dim_ = ev.dim_;

				delete[] mag_;
        		mag_ = ev.mag_;
        		ev.mag_ = nullptr;
			}
			return *this;
		}

		/** get the number of dimension for the euclidean vector
	    * 
	    @param 
	    @return number of dimension
		*/
		unsigned int EuclideanVector::getNumDimensions() const {
			return dim_;
		}

		/** get the magnitude of the given dimension
	    * 
	    @param d index of the dimension
	    @return the magnitude in d th dimension
		*/
		double EuclideanVector::get(unsigned int d) const{
			return mag_[d];
		}

		/** get the norm value of euclidean vector
	    * 
	    @param 
	    @return the norm of euclidean vector
		*/
		double EuclideanVector::getEuclideanNorm() const{
			double norm = 0;
			for (unsigned int i=0;i<dim_; ++i)
				norm += mag_[i]*mag_[i];
			return sqrt(norm);
		}

		/** Create unit vector of the same dimension
	    * 
	    @param 
	    @return the created unit vector
		*/
		EuclideanVector& EuclideanVector::createUnitVector() const{
			std::vector<double> norMag(dim_);
			double norm = this->getEuclideanNorm();
			for (unsigned int i=0;i<dim_;++i) 
				norMag[i] = mag_[i] / norm;
			EuclideanVector *unitVector = new EuclideanVector(norMag.begin(), norMag.end());
			return *unitVector;
		}


		double& EuclideanVector::operator[] (int index){
			return mag_[index];
		}

		double EuclideanVector::operator[] (int index) const{
			return mag_[index];
		}

		EuclideanVector& EuclideanVector::operator+=(const EuclideanVector &ev){
			if (dim_ != ev.getNumDimensions()) {
				std::cerr << "Error: Incompatible dimensions." << std::endl;
				abort();
			}
			for (unsigned int i=0; i<dim_; ++i)
				mag_[i] += ev.get(i);
			return *this;
		}

		EuclideanVector& EuclideanVector::operator-=(const EuclideanVector &ev){
			if (dim_ != ev.getNumDimensions()) {
				std::cerr << "Error: Incompatible dimensions." << std::endl;
				abort();
			}
			for (unsigned int i=0; i<dim_; ++i)
				mag_[i] -= ev.get(i);
			return *this;
		}

		EuclideanVector& EuclideanVector::operator*=(const double rhs){
			for (unsigned int i=0; i<dim_; ++i)
				mag_[i] *= rhs;
			return *this;
		}

		EuclideanVector& EuclideanVector::operator/=(const double rhs){
			for (unsigned int i=0; i<dim_; ++i)
				mag_[i] /= rhs;
			return *this;
		}
	
		EuclideanVector::operator std::vector<double>() const{
			std::vector<double> arr(dim_);
			for (unsigned int i=0;i<dim_;++i) {
				arr[i] = mag_[i];
			}
			return arr;
		}

		EuclideanVector::operator std::list<double>() const{
			std::list<double> arr;
			for (unsigned int i=0;i<dim_;++i) {
				arr.push_back(mag_[i]);
			}
			return arr;
		}


		/*  Friends
		 *  
		 * */
		bool operator==(const EuclideanVector &lhs, const EuclideanVector &rhs){
			if (lhs.getNumDimensions() != rhs.getNumDimensions()) return false;
			std::vector<double> l = lhs;
			std::vector<double> r = rhs;
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				if (lhs.get(i)!=rhs.get(i)) return false;
			
			return true;
		}

		bool operator!=(const EuclideanVector &lhs, const EuclideanVector &rhs){
			return !(lhs==rhs);
		}

		EuclideanVector& operator+(const EuclideanVector &lhs, const EuclideanVector &rhs){
			if (lhs.getNumDimensions() != rhs.getNumDimensions()){
				std::cerr << "Error: Incompatible dimensions." << std::endl;
				abort();
			}
			std::vector<double> l = lhs;
			std::vector<double> r = rhs;
			std::vector<double> v(lhs.getNumDimensions());
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				v[i] = l[i] + r[i];
			
			EuclideanVector *ret = new EuclideanVector{v.begin(),v.end()};
			return *ret;
		}
		

		EuclideanVector& operator-(const EuclideanVector &lhs, const EuclideanVector &rhs){
			if (lhs.getNumDimensions() != rhs.getNumDimensions()){
				std::cerr << "Error: Incompatible dimensions." << std::endl;
				abort();
			}
			std::vector<double> l = lhs;
			std::vector<double> r = rhs;
			std::vector<double> v(lhs.getNumDimensions());
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				v[i] = l[i] - r[i];
			
			EuclideanVector *ret = new EuclideanVector{v.begin(),v.end()};
			return *ret;
		}


		double operator*(const EuclideanVector &lhs, const EuclideanVector &rhs){
			if (lhs.getNumDimensions() != rhs.getNumDimensions()){
				std::cerr << "Error: Incompatible dimensions." << std::endl;
				abort();
			}
			std::vector<double> l = lhs;
			std::vector<double> r = rhs;
			double ret = 0;
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				ret += l[i] * r[i];

			return ret;
		}

		EuclideanVector& operator*(const EuclideanVector &lhs, double rhs){
			std::vector<double> l = lhs;
			std::vector<double> v(lhs.getNumDimensions());
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				v[i] = l[i] * rhs;
			
			EuclideanVector *ret = new EuclideanVector{v.begin(),v.end()};
			return *ret;			
		}

		EuclideanVector& operator/(const EuclideanVector &lhs, double rhs){
			std::vector<double> l = lhs;
			std::vector<double> v(lhs.getNumDimensions());
			for (unsigned int i=0; i<lhs.getNumDimensions(); ++i)
				v[i] = l[i] / rhs;
			EuclideanVector *ret = new EuclideanVector{v.begin(),v.end()};
			return *ret;			
		}

		std::ostream& operator<<(std::ostream &os, const EuclideanVector &ev){
			std::stringstream ss;
			unsigned int i;
			if(ev.dim_ > 0){
				for(i=0; i<ev.dim_-1; ++i){
					ss << ev.mag_[i];
					ss << " ";
				}
				ss << ev.mag_[i];
			}
		   	os << "[" << ss.str() << "]";
			return os;
		}
}


