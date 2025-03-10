#include "dynamic_vector.hpp"


dynamic_vector::dynamic_vector(size_t N, double init[]) {
    this->ndim = N;
    this->coord = new double[N];
    for (size_t i=0; i<N; i++) {
        this->coord[i] = init[i];
    }
}
dynamic_vector::dynamic_vector(size_t N, double xx) {
    this->ndim = N;
    this->coord = new double[N];
    for (size_t i=0; i<N; i++) {
        this->coord[i] = xx;
    }
}

dynamic_vector::~dynamic_vector() {
    delete[] this->coord;
}

void dynamic_vector::set(size_t i, const double& a) {
    this->coord[i] = a;
}

double dynamic_vector::operator[](size_t i) const {
    return this->coord[i];
}

double& dynamic_vector::operator[](size_t i) {
    return this->coord[i];
}

dynamic_vector::dynamic_vector(const dynamic_vector &V) {
    this->ndim = V.ndim;
    this->coord = new double[V.ndim];
    for (size_t i=0; i<V.ndim; i++) {
        this->coord[i] = V.coord[i];
    }
}

dynamic_vector& dynamic_vector::operator=(const dynamic_vector& other) {
    this->ndim = other.ndim;
    delete this->coord;
    this->coord = other.coord;

    return *this;
}

dynamic_vector operator*(const double &a, const dynamic_vector &V) {
    dynamic_vector ret = dynamic_vector(V);
    for (size_t i=0; i<V.ndim; i++) {
        ret.coord[i] *= a;
    }
    return ret;
}

dynamic_vector operator*(const dynamic_vector &V, const double &a) {return a*V;}

const dynamic_vector dynamic_vector::operator=(const double& alpha) {
    for (size_t i=0; i<this->ndim; i++) {
        this->coord[i] = alpha;
    }
    return *(this);
}

dynamic_vector& dynamic_vector::operator+=(const dynamic_vector &V) {
    for (size_t i=0; i<this->ndim; i++) {
        this->coord[i] += V.coord[i];
    }
    return *(this);
}

dynamic_vector& dynamic_vector::operator-=(const dynamic_vector &V) {
    for (size_t i=0; i<this->ndim; i++) {
        this->coord[i] -= V.coord[i];
    }
    return *(this);
}

dynamic_vector dynamic_vector::operator+() const {
    return dynamic_vector(*(this));
}

dynamic_vector dynamic_vector::operator-() const {
    return -1 * (*(this));
}

dynamic_vector operator+(const dynamic_vector& V, const dynamic_vector& W) {
    dynamic_vector R = dynamic_vector(W);
    for (size_t i =0; i<V.ndim; i++ ) {
        R.coord[i] += V[i];
    }
    return R;
}

dynamic_vector operator-(const dynamic_vector& V, const dynamic_vector& W) {return V + -W;}


std::ostream& operator<<(std::ostream& o, const dynamic_vector &V) {
    o << "(";
    for (size_t i=0; i<V.ndim-1; i++) {
        o << V.coord[i] << ", ";
    }
    o << V.coord[V.ndim-1] << ")";
    return o;
}
