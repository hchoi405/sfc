#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <cassert>
#include <cstdint>

using namespace boost::multiprecision;

namespace sfc {
template <typename T, int nDimensions>
class Vector {
   public:
    Vector() {}

    Vector(T val) {
        for (auto& el : v) el = val;
    }

    // variadic template constructor for varying size of vector
    template <typename... Args,
              typename std::enable_if<std::conjunction<std::is_convertible<Args, T>...>::value &&
                                      sizeof...(Args) == nDimensions>::type>
    Vector(Args... values) : v{values...} {}

    Vector(std::array<T, nDimensions>& arr) : v(arr) {}
    Vector(std::array<T, nDimensions>&& arr) : v(arr) {}

    // Casting operator
    operator std::array<T, nDimensions>() { return v; }

#ifdef PBRT_CORE_PBRT_H
    Vector(pbrt::Point2<T> p) : v{p.x, p.y} {}
    template <typename U>
    Vector(pbrt::Point2<U> p) : v{T(p.x), T(p.y)} {}

    Vector(pbrt::Point3<T> p) : v{p.x, p.y, p.z} {}
    template <typename U>
    Vector(pbrt::Point3<U> p) : v{T(p.x), T(p.y), T(p.z)} {}

    Vector(pbrt::Point2<T> p1, pbrt::Point3<T> p2) : v{p1.x, p1.y, p2.x, p2.y, p2.z} {}
    template <typename U>
    Vector(pbrt::Point2<U> p1, pbrt::Point3<U> p2)
        : v{T(p1.x), T(p1.y), T(p2.x), T(p2.y), T(p2.z)} {}
#endif
    T& operator[](int ind) {
        assert(ind < nDimensions);
        return v[ind];
    }
    T operator[](int ind) const {
        assert(ind < nDimensions);
        return v[ind];
    }
    void set(int i, T val) { this->v[i] = val; }
    Vector& operator+=(const Vector& rhs) {
        for (int i = 0; i < nDimensions; ++i) {
            v[i] += rhs[i];
        }
        return *this;
    }
    Vector& operator-=(const Vector& rhs) {
        for (int i = 0; i < nDimensions; ++i) {
            v[i] -= rhs[i];
        }
        return *this;
    }
    Vector& operator*=(const Vector& rhs) {
        for (int i = 0; i < nDimensions; ++i) {
            v[i] *= rhs[i];
        }
        return *this;
    }
    Vector& operator/=(const Vector& rhs) {
        for (int i = 0; i < nDimensions; ++i) {
            v[i] /= rhs[i];
        }
        return *this;
    }
    bool operator!=(const Vector& rhs) {
        for (int i = 0; i < nDimensions; ++i)
            if (v[i] != rhs[i]) return true;
        return false;
    }

    Vector operator+(const Vector& rhs) const {
        Vector ret(*this);
        return ret += rhs;
    }
    Vector operator+(const T& rhs) const {
        Vector ret(*this);
        return ret += Vector(rhs);
    }
    Vector operator-(const Vector& rhs) const {
        Vector ret(*this);
        return ret -= rhs;
    }
    Vector operator-(const T rhs) const {
        Vector ret(*this);
        return ret -= Vector(rhs);
    }
    Vector operator*(const Vector& rhs) const {
        Vector ret(*this);
        return ret *= rhs;
    }
    Vector operator/(const Vector& rhs) const {
        Vector ret(*this);
        return ret /= rhs;
    }

    template <typename U>
    Vector operator*(U rhs) const {
        Vector ret(*this);
        for (int i = 0; i < nDimensions; ++i) {
            ret.set(i, v[i] * rhs);
        }
        return ret;
    }
    template <typename U>
    Vector& operator*=(const U& rhs) {
        for (int i = 0; i < nDimensions; ++i) {
            v[i] *= rhs;
        }
        return *this;
    }
    bool operator<(const Vector<T, nDimensions>& rhs) const {
        for (int i = 0; i < nDimensions; ++i) {
            if (v[i] >= rhs[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator<=(const Vector<T, nDimensions>& rhs) const {
        for (int i = 0; i < nDimensions; ++i) {
            if (v[i] > rhs[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator>(const Vector<T, nDimensions>& rhs) const {
        for (int i = 0; i < nDimensions; ++i) {
            if (v[i] <= rhs[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator>=(const Vector<T, nDimensions>& rhs) const {
        for (int i = 0; i < nDimensions; ++i) {
            if (v[i] < rhs[i]) {
                return false;
            }
        }
        return true;
    }

    std::string toString() const {
        auto str = std::to_string(v[0]);
        for (int i = 0; i < nDimensions; ++i) {
            str += ", " + v[i];
        }
        return str;
    }

    template <typename UNSIGNED_INTEGER>
    Vector<UNSIGNED_INTEGER, nDimensions> floor() const {
        Vector<UNSIGNED_INTEGER, nDimensions> ret;
        for (int i = 0; i < nDimensions; ++i) {
            ret.set(i, UNSIGNED_INTEGER(this->v[i]));
        }
        return ret;
    }

    template <typename U, int nD>
    friend std::ostream& operator<<(std::ostream&, const Vector<U, nD>&);

   private:
    std::array<T, nDimensions> v;
};

typedef Vector<float, 2> Vector2f;
typedef Vector<uint64_t, 2> Vector2ui;

typedef Vector<float, 3> Vector3f;
typedef Vector<uint64_t, 3> Vector3ui;

typedef Vector<float, 5> Vector5f;
typedef Vector<uint64_t, 5> Vector5ui;

template <typename T, int nDimensions>
std::ostream& operator<<(std::ostream& os, const sfc::Vector<T, nDimensions>& v) {
    os << "[" << v[0];
    for (int i = 1; i < nDimensions; ++i) {
        os << ", " << v[i];
    }
    os << "]" << std::endl;
    return os;
}
}  // namespace sfc