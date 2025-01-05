#include <iostream>
#include <maxwell.hpp>

constexpr scalar_t tmax = 1e+0;
constexpr scalar_t dt   = 1e-3;

constexpr scalar_t xmin = -1e+0;
constexpr scalar_t ymin = -1e+0;
constexpr scalar_t xmax =  1e+0;
constexpr scalar_t ymax =  1e+0;

constexpr scalar_t dx = 2e-3;
constexpr scalar_t dy = 2e-3;

constexpr scalar_t mu  = 1.0;
constexpr scalar_t eps = 1.0;

const std::string ftmp = "out/";
const std::string fout = "out/dat";

int
main(void) {

  class maxwell m({xmin, ymin}, {xmax, ymax}, {dx, dy}, mu, eps, ftmp);


  m.assert__courant(dt);

  for (size_t nt{0}; nt < static_cast<size_t>(tmax / dt); nt++) {

    const scalar_t t = static_cast<scalar_t>(nt) * dt;

    std::cout << "\rGenerating data : "
              << " timesteps : " << t << " / " << tmax
              << std::flush;

    m.init__point_charge_dipole(t);

    m.solve(dt);

    m.setbounds__PML(dt);

    m.save(nt);

  }

  std::cout << "\n"
            << "Compiling hdf5 data"
            << std::endl;

  h5concat(ftmp, fout, tmax, dt,
    { tmax, dt, xmin, ymin, xmax, ymax, dx, dy, mu, eps }
  );

  std::cout << "C++ program complete."
            << std::endl;

  return 0;

}
