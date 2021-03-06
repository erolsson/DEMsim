/*
  This file is part of the DEM program DEMsim. The code is free to run and
  develop after getting permission from the original developers Erik Olsson
  and Per-Lennart Larsson, KTH Solid Mechanics, Stockholm, Sweden

  erolsson@kth.se
  pelle@kth.se
 */


/*
  Small 3D-vector class with common vector operations

*/

#ifndef DEMSIM6_0_VEC3_H
#define DEMSIM6_0_VEC3_H



#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

namespace DEM {
    class Vec3 {
    public:

        //Cartesian coordinates

        //Empty default constructor and a constructor initializing with the components;
        Vec3() = default;

        Vec3(double x, double y, double z)
                :data_{x, y, z} { }

        [[nodiscard]] const double& x() const {
            return data_[0];
        }

        [[nodiscard]] const double& y() const {
            return data_[1];
        }

        [[nodiscard]] const double& z() const {
            return data_[2];
        }

        double& x() {
            return data_[0];
        }

        double& y() {
            return data_[1];
        }

        double& z() {
            return data_[2];
        }

        const double& operator[](std::size_t axis) const {
            return data_[axis];
        }

        double& operator[](std::size_t axis) {
            return data_[axis];
        }

        inline Vec3& normalize()
        {
            double l = inv_length();
            x() *= l;
            y() *= l;
            z() *= l;
            return *this;
        }

        inline void flip()
        {
            x() = -x();
            y() = -y();
            z() = -z();
        }

        [[nodiscard]] inline double length() const
        {
            return sqrt(x()*x() + y()*y() + z()*z());
        }

        [[nodiscard]] inline Vec3 normal() const
        {
            double inv = inv_length();
            return *this*inv;
        }

        //Common vector operations, should be self-explaining
        inline Vec3 operator-() const
        {
            return Vec3(-x(), -y(), -z());
        }

        [[nodiscard]] inline bool is_zero() const
        {
            return x()==0.0 && y()==0.0 && z()==0.0;
        }

        inline Vec3 operator-(const Vec3& vec) const
        {
            return Vec3(x()-vec.x(), y()-vec.y(), z()-vec.z());
        }

        inline Vec3 operator+(const Vec3& vec) const
        {
            return Vec3(x()+vec.x(), y()+vec.y(), z()+vec.z());
        }

        inline Vec3 operator*(double factor) const
        {
            return Vec3(x()*factor, y()*factor, z()*factor);
        }

        inline Vec3 operator/(double factor) const
        {
            double inv = 1/factor;
            return *this*inv;
        }

        inline Vec3& operator*=(double factor)
        {
            x() *= factor;
            y() *= factor;
            z() *= factor;
            return *this;
        }

        inline Vec3& operator/=(double factor)
        {
            x() /= factor;
            y() /= factor;
            z() /= factor;
            return *this;
        }

        inline Vec3& operator+=(const Vec3& rhs)
        {
            x() += rhs.x();
            y() += rhs.y();
            z() += rhs.z();
            return *this;
        }

        inline Vec3& operator-=(const Vec3& rhs)
        {
            x() -= rhs.x();
            y() -= rhs.y();
            z() -= rhs.z();
            return *this;
        }

        inline void set_zero()
        {
            x() = 0.;
            y() = 0.;
            z() = 0.;
        }

        [[nodiscard]] std::string csv_string() const
        {
            std::stringstream ss;
            ss << x() << ", " << y() << ", " << z();
            return ss.str();
        }

    private:
        double data_[3]{0., 0., 0.};

        [[nodiscard]] double inv_length() const
        {
            return 1./length();
        }

    };

    //Vector functions that aren't part of the class
    inline double dot_product(const Vec3& a, const Vec3& b)
    {
        return a.x()*b.x() + a.y()*b.y()+ a.z()*b.z();
    }

    inline bool operator==(const Vec3& a, const Vec3& b)
    {
        return (a.x()==b.x() && a.y()==b.y() && a.z()==b.z());
    }

    inline Vec3 cross_product(const Vec3& a, const Vec3& b)
    {
        return Vec3(a.y()*b.z()-a.z()*b.y(), a.z()*b.x()-a.x()*b.z(), a.x()*b.y()-a.y()*b.x());
    }

    //Cute printing
    inline std::ostream& operator<<(std::ostream& os, const Vec3& vec)
    {
        os << '(' << vec.x() << ',' << vec.y() << ',' << vec.z() << ')';
        return os;
    }

    inline Vec3 operator*(double factor, const Vec3& a)
    {
        return Vec3(a.x()*factor, a.y()*factor, a.z()*factor);
    }



}

#endif //DEMSIM6_0_VEC3_H
