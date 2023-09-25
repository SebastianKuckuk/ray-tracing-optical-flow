#pragma once

struct Vec3 {
    double x, y, z;

    // Constructor
    Vec3(double xVal, double yVal, double zVal)
            : x(xVal), y(yVal), z(zVal) {}

    explicit Vec3(double scalar)
            : x(scalar), y(scalar), z(scalar) {}

    // Addition operator
    Vec3 operator+(const Vec3 &other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    // Subtraction operator
    Vec3 operator-(const Vec3 &other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    // Scalar multiplication operator
    Vec3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    Vec3 operator*(const Vec3 &other) const {
        return {x * other.x, y * other.y, z * other.z};
    }

    // Compound assignment addition operator (+=)
    Vec3 &operator+=(const Vec3 &other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Compound assignment subtraction operator (-=)
    Vec3 &operator-=(const Vec3 &other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // Compound assignment scalar multiplication operator (*=)
    Vec3 &operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    // Scalar division operator
    Vec3 operator/(double scalar) const {
        if (scalar != 0.0) {
            return {x / scalar, y / scalar, z / scalar};
        } else {
            // Handle division by zero gracefully or throw an exception
            throw std::runtime_error("Division by zero");
        }
    }

    // Dot product
    [[nodiscard]] double dot(const Vec3 &other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Magnitude (length) of the vector
    [[nodiscard]] double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Normalize the vector (make it a unit vector)
    [[nodiscard]] Vec3 normalize() const {
        double mag = magnitude();
        if (mag != 0.0) {
            return *this / mag;
        } else {
            // Handle zero magnitude gracefully or throw an exception
            throw std::runtime_error("Normalization of zero vector");
        }
    }

    // Print the vector
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
};

double dot(const Vec3 &a, const Vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Scalar-vector multiplication as a free function
Vec3 operator*(double scalar, const Vec3 &vec) {
    return {scalar * vec.x, scalar * vec.y, scalar * vec.z};
}
