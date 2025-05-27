#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include "../core/konal_constants.hpp"
#include <cmath> // For std::sqrt, std::abs, etc.
#include <array>
#include <numeric> // For std::inner_product

namespace Ckonal {
namespace utils {

template<typename T, size_t N>
T magnitude(const std::array<T, N>& vec) {
    T sum_sq = 0;
    for (size_t i = 0; i < N; ++i) {
        sum_sq += vec[i] * vec[i];
    }
    return std::sqrt(sum_sq);
}

template<typename T, size_t N>
std::array<T, N> normalize(const std::array<T, N>& vec) {
    std::array<T, N> result;
    T mag = magnitude(vec);
    if (std::abs(mag) < 1e-9) { // Avoid division by zero or very small numbers
        // Return zero vector or handle error as appropriate
        for (size_t i = 0; i < N; ++i) result[i] = 0;
        return result;
    }
    for (size_t i = 0; i < N; ++i) {
        result[i] = vec[i] / mag;
    }
    return result;
}

// Function to solve quadratic equation ax^2 + bx + c = 0
// Returns a vector of real roots.
inline std::vector<real_t> solve_quadratic(real_t a, real_t b, real_t c) {
    std::vector<real_t> roots;
    if (std::abs(a) < 1e-9) { // Linear equation
        if (std::abs(b) > 1e-9) {
            roots.push_back(-c / b);
        }
        // if b is also zero, it's either 0=c (no solution if c!=0) or 0=0 (infinite solutions)
        // For Eikonal, a=0 implies a different update rule, or no contribution from that dimension
        return roots;
    }

    real_t discriminant = b * b - 4 * a * c;

    if (discriminant > -1e-9) { // Allow for small negative due to precision
        if (discriminant < 0) discriminant = 0; // Per original code's hack for negative discriminant
        
        real_t sqrt_discriminant = std::sqrt(discriminant);
        // See Numerical Recipes for stable quadratic formula
        // Or use the direct approach as in the original code if it's simpler for now
        // roots.push_back((-b + sqrt_discriminant) / (2 * a));
        // roots.push_back((-b - sqrt_discriminant) / (2 * a));

        // Python code used: new = (-b + sqrt(b**2 - 4*a*c)) / (2*a)
        // This selects one root. For Eikonal, usually the larger, causally consistent root.
        // Let's return both for now, solver can pick.
        // More robust formula to avoid catastrophic cancellation:
        real_t r1, r2;
        if (b >= 0) {
            r1 = (-b - sqrt_discriminant) / (2 * a);
            if (std::abs(r1) > 1e-9) { // Avoid division by zero if r1 is zero
                 r2 = c / (a * r1);
            } else { // r1 is zero, means c must be zero for this root
                 r2 = -b / a; // The other root
            }
        } else {
            r1 = (-b + sqrt_discriminant) / (2 * a);
             if (std::abs(r1) > 1e-9) {
                 r2 = c / (a * r1);
            } else {
                 r2 = -b / a;
            }
        }
        roots.push_back(r1);
        roots.push_back(r2);
    }
    // No real roots if discriminant is significantly negative
    return roots;
}


} // namespace utils
} // namespace Ckonal

#endif // MATH_UTILS_HPP