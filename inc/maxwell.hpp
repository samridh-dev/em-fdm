#ifndef MAXWELL_HPP
#define MAXWELL_HPP

#include <vector>
#include <array>
#include <string>

using scalar_t = double;
using vector_t = std::array<scalar_t, 2>;

template <typename T>
class field {

  public:

    field(const std::array<scalar_t, 2>& _min,
          const std::array<scalar_t, 2>& _max,
          const std::array<scalar_t, 2>& _h);

    size_t 
    size(const size_t n) const;

    T& 
    operator[](const size_t idx);

    const T& 
    operator[](const size_t idx) const;

  private:
    std::array<size_t, 2> N;
    std::vector<T> data;

};

using sfield = field<scalar_t>;
using vfield = field<std::array<scalar_t, 2>>;

class maxwell {
  
  public:

    maxwell(const std::array<scalar_t, 2>& _min,
            const std::array<scalar_t, 2>& _max,
            const std::array<scalar_t, 2>& _h,
            const scalar_t mu, const scalar_t eps,
            const std::string& fname);

    void
    assert__courant(const scalar_t dt);

    void
    init__point_charge_static(const scalar_t dt, const scalar_t q = 1e2);

    void
    init__point_charge_sinusoidal(const scalar_t t, 
                                  const scalar_t q = 1e2,
                                  const scalar_t r = 2e-1,
                                  const scalar_t omega = 2e1);
    void
    init__point_charge_dipole(const scalar_t t, 
                              const scalar_t q = 1e4,
                              const scalar_t r = 1e-2,
                              const scalar_t omega = 1e2);
    void
    setbounds__PML(const scalar_t dt,
                   const scalar_t sigma = 5e1,
                   const scalar_t alpha = 0.1,
                   const size_t width = 10);
    void
    setbounds__PEC(const scalar_t dt);

    void
    solve(const scalar_t dt);

    void
    save(const size_t i);

  private:

    vfield E;
    sfield B;

    sfield p;
    vfield J;

    sfield V0;
    sfield V1;
    sfield V2;

    vfield A0;
    vfield A1;
    vfield A2;

    const std::array<scalar_t, 2> min;
    const std::array<scalar_t, 2> max;
    const std::array<scalar_t, 2> h;

    const scalar_t mu;
    const scalar_t eps;

    const std::string fname;

};

struct metadata {
  double tmax;
  double dt;
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  double dx;
  double dy;
  double mu;
  double eps;
};

void 
h5concat(const std::string &fin, const std::string &fout,
         scalar_t t_max, scalar_t dt, struct metadata md);

#endif // MAXWELL_HPP
