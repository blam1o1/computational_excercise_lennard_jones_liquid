#include <iostream>
#include "../include/vec3.h"

int main() {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};

    // Print a and b
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";

    // Test addition
    Vec3 c = a + b;
    std::cout << "a + b = " << c << "\n";
    // Test subtraction
    Vec3 d = a - b;
    std::cout << "a - b = " << d << "\n";
    std::cout << "b - a = " << b - a << "\n";
    // Test scalar multiplication
    Vec3 e = a * 2.0;
    std::cout << "a * 2.0 = " << e << "\n";
    // Test scalar division
    Vec3 f = a / 2.0;
    std::cout << "a / 2.0 = " << f << "\n";
    // Test dot product
    double g = dot(a, b);
    std::cout << "dot(a, b) = " << g << "\n";
    // Test cross product
    Vec3 h = cross(a, b);
    std::cout << "cross(a, b) = " << h << "\n";
    // Test magnitude
    double i = magnitude(a);
    std::cout << "magnitude(a) = " << i << "\n";
    // Test magnitude squared
    double j = magnitude_squared(a);
    std::cout << "magnitude_squared(a) = " << j << "\n";

    return 0;
}