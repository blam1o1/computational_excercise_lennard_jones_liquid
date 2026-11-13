// header file for vec3.cpp
#pragma once
#include <ostream>
#include <cmath>

// no need to implement vec3.cpp, we can define everythign inline in vec3.h

// Vec3 struct to store 3D vectors and some operations 
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;


    Vec3& operator+=(const Vec3& b) {
        x += b.x;
        y += b.y;
        z += b.z;
        return *this;
    }
    Vec3& operator-=(const Vec3& b) {
        x -= b.x;
        y -= b.y;
        z -= b.z;
        return *this;
    }
    Vec3& operator*=(double s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    Vec3& operator/=(double s) {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }
};


// Vector operations

// Addition
inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

// Subtraction 
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

// Scalar multiplication
inline Vec3 operator*(const Vec3& a, double s) {
    return Vec3{a.x * s, a.y * s, a.z * s};
}
// Scalar multiplication reversed
inline Vec3 operator*(double s, const Vec3& a) {
    return Vec3{s * a.x, s * a.y, s * a.z};
}

// Scalar division
inline Vec3 operator/(const Vec3& a, double s) {
    return Vec3{a.x / s, a.y / s, a.z / s};
}

// Dot product
inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// Magnitude
inline double magnitude(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

//Magnitude squared
inline double magnitude_squared(const Vec3& a) {
    return dot(a, a);
}

//Printing the vector 
inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}