#include "gsl/gsl_sf.h"

#include <cmath>
#include <time.h>
#include "iostream"
#include "fstream"
#include "array"
#include "string"


class MCEV_benchmark {
public:
  double f_large_z(const double z) {
    constexpr double one = 1.;
    double f = one;
    double z_pow = one; 
    double theta_shft = one - theta;
    double omega_shft = omega - theta ;
    std::array<double, 64> d;
    std::array<double, 65> coef;
    d[0] = one;
    coef[0] = one;
    coef[1] = theta_shft * omega_shft;
    double res = one;
    double s_double = one;
    for (size_t s = 1 ; s < N && fabs(res) > epsilon; ++s, s_double += one) {
      double d_new = one;
      double rhs_term = 0.;
      double j_double = 0.;
      for (size_t j = 0; j < s; ++j, j_double += one) {
        d_new *= ((theta_shft + j_double) * (omega_shft + j_double))  / (j_double + one);
//       i = s - j
        rhs_term += d[j] * coef[s - j];
      }
      double tmp0 = (theta_shft + s_double) * (omega_shft + s_double);
      coef[s + 1] = d_new * tmp0 / (s_double + one);
      d_new *= tmp0 / (theta_shft * omega_shft);
      d_new -= rhs_term;
      d[s] = d_new;
      z_pow /= z;
      res = d_new * z_pow; 
      f += res;
    }
    return - f * theta_shft / z;
  }

  double f_small_z(const double z) {
    constexpr double one = 1.;
    double f = one;
    double z_pow = one; 
    std::array<double, 64> c;
    std::array<double, 65> coef;
    c[0] = one;
    coef[0] = one;
    coef[1] = theta / omega;
    double theta_shft = theta - one;
    double res = one;
    double s_double = one;
    for (size_t s = 1 ; s < N && fabs(res) > epsilon; ++s, s_double += one) {
      double c_new = one;
      double rhs_term = 0.;
      double j_double = 0.;
      for (size_t j = 0; j < s; ++j, j_double += one) {
        c_new *= (theta_shft + j_double) / ((j_double + one) * (omega + j_double));
//       i = s - j
        rhs_term += c[j] * coef[s - j];
      }
      coef[s + 1] = c_new * ((theta_shft + s_double) / ((s_double + one) * (omega + s_double))) 
                          * ((theta_shft + s_double + one) / theta_shft);
      c_new -= rhs_term;
      c[s] = c_new;
      z_pow *= z;
      res = c_new * z_pow;
      f += res;
    }
    return f;
  }

  double f_gsl(const double z) {
      double numerator  = gsl_sf_hyperg_1F1(theta - 1, omega, z);
      double denominator = gsl_sf_hyperg_1F1(theta, omega, z);
      return numerator / denominator;
  }
  double time_gsl(const double z) {
    clock_t tStart = clock();
    double f = 0.;
    for (size_t i = 0; i < it_num; ++i) {
      f = f_gsl(z);
    }
    return (double)(clock() - tStart)/CLOCKS_PER_SEC;
  }
  double time_small_z(const double z) {
    clock_t tStart = clock();
    double f = 0.;
    for (size_t i = 0; i < it_num; ++i) {
      f = f_small_z(z);
    }
    return (double)(clock() - tStart)/CLOCKS_PER_SEC;
  }
  double time_large_z(const double z) {
    clock_t tStart = clock();
    double f = 0.;
    for (size_t i = 0; i < it_num; ++i) {
      f = f_large_z(z);
    }
    return (double)(clock() - tStart)/CLOCKS_PER_SEC;
  }

  void comparison_test_small_z(const std::string& f_name, const double z_min, const double z_max) {
    std::ofstream dump(f_name);
    for (double z = z_min ; z<=z_max; z += 0.001) {
      dump << time_small_z(z) <<"," <<time_gsl(z) <<"," <<z<<std::endl; 
    }
  }

  void comparison_test_large_z(const std::string& f_name, const double z_min, const double z_max) {
    std::ofstream dump(f_name);
    for (double z = z_min ; z<=z_max; z += 1.) {
      dump << time_large_z(z) <<"," <<time_gsl(z) <<"," <<z<<std::endl; 
    }
  }


  MCEV_benchmark(double lambda, double eta) : theta(0.5 + eta - lambda), omega(1 + 2. * eta) {};
  MCEV_benchmark() = default;
private: 
  double theta;
  double omega;
  size_t N = 64;  
  double epsilon = 1e-10;
  size_t it_num = 10000;  
};

int main ()
{
  double lambda = 1.2;
  double eta = 2.12;
  double z = 0.07;
  MCEV_benchmark bench (lambda, eta);

  std::cout.precision(10);
  std::cout << "GSL time: " <<  bench.time_gsl(z) << std::endl;  
  std::cout << "Our time: " <<  bench.time_small_z(z) << std::endl;  
  std::cout<<"Our value: "<< bench.f_small_z(z) << std::endl;
  std::cout<<"GSL value: "<< bench.f_gsl(z) << std::endl;

  return 0;
}
