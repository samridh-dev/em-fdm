#include <maxwell.hpp>

#include <omp.h>

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
///                                                                         ///
///                                                                         ///
/// CLASS                                                                   ///
///                                                                         ///
///                                                                         ///
///////////////////////////////////////////////////////////////////////////////

template<typename T>
field<T>::field(const std::array<scalar_t, 2>& _min,
                      const std::array<scalar_t, 2>& _max,
                      const std::array<scalar_t, 2>& _h)
  : N({
      static_cast<size_t>((_max[0] - _min[0]) / _h[0]), 
      static_cast<size_t>((_max[1] - _min[1]) / _h[1])
    }),
    data(N[0] * N[1]) {}

template<typename T>
size_t 
field<T>::size(size_t n) const {
  assert(n < 2 && "dimension out of range");
  return N[n];
}

template<typename T>
T& 
field<T>::operator[](const size_t idx) { 
  return data[idx]; 
}

template<typename T>
const T& 
field<T>::operator[](const size_t idx) const { 
  return data[idx]; 
}

template class field<scalar_t>;
template class field<vector_t>;

///////////////////////////////////////////////////////////////////////////////
///                                                                         ///
///                                                                         ///
/// IMPLEMENTATION                                                          ///
///                                                                         ///
///                                                                         ///
///////////////////////////////////////////////////////////////////////////////

maxwell::maxwell(const std::array<scalar_t, 2>& _min,
                 const std::array<scalar_t, 2>& _max,
                 const std::array<scalar_t, 2>& _h,
                 const scalar_t _mu, const scalar_t _eps,
                 const std::string& _fname)

  : E(_min, _max, _h), B(_min, _max, _h), 
    p(_min, _max, _h), J(_min, _max, _h), 
    V0(_min, _max, _h), V1(_min, _max, _h), V2(_min, _max, _h),
    A0(_min, _max, _h), A1(_min, _max, _h), A2(_min, _max, _h),
    min(_min), max(_max), h(_h),
    mu(_mu), eps(_eps), fname(_fname)

{}

void
maxwell::assert__courant(const scalar_t dt) {

  const scalar_t c = 1 / std::sqrt(mu * eps);

  const scalar_t courant = 1 / (c * std::sqrt( 1.0 / (h[0] * h[0]) + 
                                               1.0 / (h[1] * h[1]) ));

  if (dt >= courant) {
    std::cout << "Warning: "
              << "Time step dt = " << dt
              << " exceeds Courant condition value of " << courant << ".\n";
  } else {
    std::cout << " Courant condition value of " 
              << "t(" << dt << ") <= " << courant << " satisfied.\n";
  }

}
void
maxwell::init__point_charge_static(const scalar_t dt, const scalar_t q) {

  (void) dt;

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  auto delta = [&](const scalar_t x, const scalar_t y) -> scalar_t {
    const scalar_t tol = 2*h[0];
    return (std::abs(x) < tol) && (std::abs(y) < tol) ? 1.0 : 0.0;
  };

  #pragma omp parallel for 
  for (size_t i = 0; i < Nx; ++i) {
  for (size_t j = 0; j < Ny; ++j) {

    const scalar_t x = min[0] + static_cast<scalar_t>(i) * h[0];
    const scalar_t y = min[1] + static_cast<scalar_t>(j) * h[1];
    p[Ny * i + j] = q * delta(x, y);
    
  }}
  

}

void
maxwell::init__point_charge_sinusoidal(const scalar_t t, 
                                       const scalar_t q,
                                       const scalar_t r,
                                       const scalar_t omega) {

  if (t == 0.00) {
    init__point_charge_static(t, q);
  }

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  auto delta = [&](const scalar_t x, const scalar_t y) -> scalar_t {
    const scalar_t tol = 1e-1;
    return (std::abs(x) < tol) && (std::abs(y) < tol) ? 1.0 : 0.0;
  };

  #pragma omp parallel for 
  for (size_t i = 0; i < Nx; ++i) {
  for (size_t j = 0; j < Ny; ++j) {

    const scalar_t x = min[0] + static_cast<scalar_t>(i) * h[0];
    const scalar_t y = min[1] + static_cast<scalar_t>(j) * h[1];
    const scalar_t del = delta(x, y - r * std::sin(omega * t));

    /* p = q del(x, y) */
    p[Ny * i + j]  = q * del;
    
    /* J = q v del(x, y) */
    J[Ny * i + j][0] = 0.0;
    J[Ny * i + j][1] = q * r * omega * std::cos(omega * t) * del;
//    J[Ny * i + j][1] = 0.0;
    
  }}
  

}

void
maxwell::init__point_charge_dipole(const scalar_t t, 
                                   const scalar_t q,
                                   const scalar_t r,
                                   const scalar_t omega) {


  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  auto delta = [&](const scalar_t x, const scalar_t y) -> scalar_t {
    const scalar_t tol = 2 * h[0];
    return (std::abs(x) < tol) && (std::abs(y) < tol) ? 1.0 : 0.0;
  };

  #pragma omp parallel for 
  for (size_t i = 0; i < Nx; ++i) {
  for (size_t j = 0; j < Ny; ++j) {

    const scalar_t x = min[0] + static_cast<scalar_t>(i) * h[0];
    const scalar_t y = min[1] + static_cast<scalar_t>(j) * h[1];
    const scalar_t del1 = delta(x + (r/2) * std::cos(omega * t),
                                y - (r/2) * std::sin(omega * t));
    const scalar_t del2 = delta(x - (r/2) * std::cos(omega * t),
                                y + (r/2) * std::sin(omega * t));

    p[Ny * i + j]  = q * (del1 - del2);

    /* J = q v del(x, y) */
    J[Ny * i + j][0] = - q * r * omega * std::sin(omega * t) * (del2 - del1);
    J[Ny * i + j][1] =   q * r * omega * std::cos(omega * t) * (del1 - del2);
//    J[Ny * i + j][0] = 0.0; J[Ny * i + j][1] = 0.0;
    
  }}

#if 0
  if (t == 0.00) {

    auto idx = [&](const size_t x, const size_t y) -> size_t {
      return Ny * x + y;
    };
    const size_t ij  = idx(i+0, j+0);
    const size_t i0j = idx(i-1, j+0);
    const size_t ij0 = idx(i+0, j-1);
    const size_t i1j = idx(i+1, j+0);
    const size_t ij1 = idx(i+0, j+1);

    size_t iter = 0;
    size_t imax = 1e2;

    while (iter++ < imax) {
      for (size_t i = 1; i < Nx -1; ++i) {
      for (size_t j = 1; j < Ny -1; ++j) {

      V0[ij] = 0.25 (V[i0j] + V[i1j] + 
      A2[ij][0] = (c*c*dt*dt) * (laplacianv(A1, 0) - (J[ij][0]/mu)) 
                + 2*A1[ij][0] - A0[ij][0];
      A2[ij][1] = (c*c*dt*dt) * (laplacianv(A1, 1) - (J[ij][1]/mu)) 
                + 2*A1[ij][1] - A0[ij][1];


      }}
    }


  }
#endif


}

void
maxwell::setbounds__PEC(const scalar_t dt) {

  (void) dt;

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  #pragma omp parallel for 
  for (size_t j = 0; j < Ny; ++j) {
    V0[(     0) * Ny + j] = V0[(     1) * Ny + j];
    V0[(Nx - 1) * Ny + j] = V0[(Nx - 2) * Ny + j];
    V1[(     0) * Ny + j] = V1[(     1) * Ny + j];
    V1[(Nx - 1) * Ny + j] = V1[(Nx - 2) * Ny + j];
    V2[(     0) * Ny + j] = V2[(     1) * Ny + j];
    V2[(Nx - 1) * Ny + j] = V2[(Nx - 2) * Ny + j];
  }

  #pragma omp parallel for 
  for (size_t i = 0; i < Nx; ++i) {
    V0[i * Ny + (     0)] = V0[i * Ny + (     1)];
    V0[i * Ny + (Ny - 1)] = V0[i * Ny + (Ny - 2)];
    V1[i * Ny + (     0)] = V1[i * Ny + (     1)];
    V1[i * Ny + (Ny - 1)] = V1[i * Ny + (Ny - 2)];
    V2[i * Ny + (     0)] = V2[i * Ny + (     1)];
    V2[i * Ny + (Ny - 1)] = V2[i * Ny + (Ny - 2)];
  }

}

void maxwell::setbounds__PML(const scalar_t dt,
                             const scalar_t sigma,
                             const scalar_t alpha,
                             const size_t width) {
  (void) alpha;

  auto getsigma = [&](const size_t d, const size_t w) -> scalar_t {
    const scalar_t dnorm = static_cast<scalar_t>(d) / static_cast<scalar_t>(w);
    return sigma * std::pow(dnorm, 3);
  };

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  #pragma omp parallel for 
  for (size_t j = 0; j < Ny; ++j) {
  for (size_t d = 0; d < width; ++d) {
    const scalar_t sl = getsigma(width - d, width);
    const scalar_t sr = getsigma(        d, width);

    V0[d * Ny + j] *= std::exp(-sl * dt);
    V1[d * Ny + j] *= std::exp(-sl * dt);
    V2[d * Ny + j] *= std::exp(-sl * dt);

    V0[(Nx - 1 - d) * Ny + j] *= std::exp(-sr * dt);
    V1[(Nx - 1 - d) * Ny + j] *= std::exp(-sr * dt);
    V2[(Nx - 1 - d) * Ny + j] *= std::exp(-sr * dt);

    A0[d * Ny + j][0] *= std::exp(-sl * dt);
    A1[d * Ny + j][0] *= std::exp(-sl * dt);
    A2[d * Ny + j][0] *= std::exp(-sl * dt);

    A0[(Nx - 1 - d) * Ny + j][0] *= std::exp(-sr * dt);
    A1[(Nx - 1 - d) * Ny + j][0] *= std::exp(-sr * dt);
    A2[(Nx - 1 - d) * Ny + j][0] *= std::exp(-sr * dt);

    A1[d * Ny + j][1] *= std::exp(-sl * dt);
    A1[d * Ny + j][1] *= std::exp(-sl * dt);
    A2[d * Ny + j][1] *= std::exp(-sl * dt);

    A1[(Nx - 1 - d) * Ny + j][1] *= std::exp(-sr * dt);
    A1[(Nx - 1 - d) * Ny + j][1] *= std::exp(-sr * dt);
    A2[(Nx - 1 - d) * Ny + j][1] *= std::exp(-sr * dt);

  }}

  #pragma omp parallel for 
  for (size_t i = 0; i < Nx; ++i) {
  for (size_t d = 0; d < width; ++d) {
    const scalar_t sb = getsigma(width - d, width);
    const scalar_t su = getsigma(        d, width);

    V0[i * Ny + d] *= std::exp(-sb * dt);
    V1[i * Ny + d] *= std::exp(-sb * dt);
    V2[i * Ny + d] *= std::exp(-sb * dt);

    V0[i * Ny + (Ny - 1 - d)] *= std::exp(-su * dt);
    V1[i * Ny + (Ny - 1 - d)] *= std::exp(-su * dt);
    V2[i * Ny + (Ny - 1 - d)] *= std::exp(-su * dt);

    A0[i * Ny + d][0] *= std::exp(-sb * dt);
    A1[i * Ny + d][0] *= std::exp(-sb * dt);
    A2[i * Ny + d][0] *= std::exp(-sb * dt);

    A0[i * Ny + (Ny - 1 - d)][0] *= std::exp(-su * dt);
    A1[i * Ny + (Ny - 1 - d)][0] *= std::exp(-su * dt);
    A2[i * Ny + (Ny - 1 - d)][0] *= std::exp(-su * dt);

    A0[i * Ny + d][1] *= std::exp(-sb * dt);
    A1[i * Ny + d][1] *= std::exp(-sb * dt);
    A2[i * Ny + d][1] *= std::exp(-sb * dt);

    A0[i * Ny + (Ny - 1 - d)][1] *= std::exp(-su * dt);
    A1[i * Ny + (Ny - 1 - d)][1] *= std::exp(-su * dt);
    A2[i * Ny + (Ny - 1 - d)][1] *= std::exp(-su * dt);

  }}

}
void
maxwell::solve(const scalar_t dt) {

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  const scalar_t c = 1.0 / std::sqrt(mu * eps);

  #define SET_INDICES_VARIABLES                                               \
    auto idx = [&](const size_t x, const size_t y) -> size_t {                \
      return Ny * x + y;                                                      \
    };                                                                        \
    const size_t ij  = idx(i+0, j+0);                                         \
    const size_t i0j = idx(i-1, j+0);                                         \
    const size_t ij0 = idx(i+0, j-1);                                         \
    const size_t i1j = idx(i+1, j+0);                                         \
    const size_t ij1 = idx(i+0, j+1);                                         \

  /* obtain potential V and A */
  {

    #pragma omp parallel for 
    for (size_t i = 1; i < Nx - 1; ++i) {
    for (size_t j = 1; j < Ny - 1; ++j) {

      SET_INDICES_VARIABLES

      auto laplacians = [&](const sfield& s) -> scalar_t {
        return (s[i0j] - 2 * s[ij] + s[i1j]) / (h[0]*h[0])
             + (s[ij0] - 2 * s[ij] + s[ij1]) / (h[1]*h[1]);
      };

      auto laplacianv = [&](const vfield& F, const size_t d) -> scalar_t {
        return (F[i0j][d]- 2 * F[ij][d] + F[i1j][d]) / (h[0]*h[0])
             + (F[ij0][d]- 2 * F[ij][d] + F[ij1][d]) / (h[1]*h[1]);
      };

      V2[ij] = (c*c*dt*dt) * (laplacians(V1) + p[ij]/eps) + 2*V1[ij] - V0[ij];
      A2[ij][0] = (c*c*dt*dt) * (laplacianv(A1, 0) - (J[ij][0]/mu)) 
                + 2*A1[ij][0] - A0[ij][0];
      A2[ij][1] = (c*c*dt*dt) * (laplacianv(A1, 1) - (J[ij][1]/mu)) 
                + 2*A1[ij][1] - A0[ij][1];

    }}

  }

  /* Calculate Electric and Magnetic Fields */
  {
    #pragma omp parallel for 
    for (size_t i = 1; i < Nx - 1; ++i) {
    for (size_t j = 1; j < Ny - 1; ++j) {
      SET_INDICES_VARIABLES
      E[ij][0] = - (V1[i1j]   - V1[i0j]  ) / (h[0]) 
                 - (A2[ij][0] - A0[ij][0]) / (2*dt);
      E[ij][1] = - (V1[ij1]   - V1[ij0]  ) / (h[1]) 
                 - (A2[ij][1] - A0[ij][1]) / (2*dt);
      B[ij] = ((A1[ij1][0] - A1[ij0][0]) / h[1]) 
            - ((A1[i1j][1] - A1[i0j][1]) / h[0]);
    }}
  }

  V0 = V1;
  V1 = V2;
  A0 = A1;
  A1 = A2;

  #undef SET_INDICES_VARIABLES

}

