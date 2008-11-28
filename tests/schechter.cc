#include "milia/schechter.h"
#include <iostream>
#include <cmath>
#include <boost/lambda/lambda.hpp>

using namespace boost::lambda;
using milia::luminosity_functions::schechter;

class Per {
public:
Per(double a) : m_a(a) {}
double operator()(double z) {return m_a;}
private:
double m_a;
};

class Har {
public:
Har(double a) : m_a(a) {}
double operator()(double z) {return m_a + z;}
private:
double m_a;
};
class Har2 {
public:
Har2(double a,double b) : m_a(a),m_b(b) {}
double operator()(double z) {return m_a * pow((1. + z),m_b) ;}
private:
double m_a;
double m_b;
};

int main() {
double ls = 8.0;
double x1 = 1.455555555555;
double x2 = 1.455555555556;
schechter a(2.,ls,-1.3,0.);
std::cout << a.integrate(ls*x1,ls*x2) << " " << 
a.integrate2(ls*x1,ls*x2) << " " << a.integrate3(ls*x1, ls*x2) << std::endl;
}
