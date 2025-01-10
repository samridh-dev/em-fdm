#include <maxwell.hpp>
#include <H5Cpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <thread>
#include <algorithm>
#include <cstdio>

static const std::string out_dir = "out/";

static struct metadata 
read_metadata(H5::H5File &file) {
  struct metadata md;
  file.openAttribute("tmax").read(H5::PredType::NATIVE_DOUBLE, &md.tmax);
  file.openAttribute("dt")  .read(H5::PredType::NATIVE_DOUBLE, &md.dt);
  file.openAttribute("xmin").read(H5::PredType::NATIVE_DOUBLE, &md.xmin);
  file.openAttribute("ymin").read(H5::PredType::NATIVE_DOUBLE, &md.ymin);
  file.openAttribute("xmax").read(H5::PredType::NATIVE_DOUBLE, &md.xmax);
  file.openAttribute("ymax").read(H5::PredType::NATIVE_DOUBLE, &md.ymax);
  file.openAttribute("dx")  .read(H5::PredType::NATIVE_DOUBLE, &md.dx);
  file.openAttribute("dy")  .read(H5::PredType::NATIVE_DOUBLE, &md.dy);
  file.openAttribute("mu")  .read(H5::PredType::NATIVE_DOUBLE, &md.mu);
  file.openAttribute("eps") .read(H5::PredType::NATIVE_DOUBLE, &md.eps);
  return md;
}

struct chunk {
  std::vector<scalar_t> E;
  std::vector<scalar_t> B;
  chunk(void) : E(), B() {}  
};

int main(int argc, char* argv[]) 
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <file.h5>\n";
    return 1;
  }
  const std::string fname = argv[1];

  try {

    H5::H5File file(fname, H5F_ACC_RDONLY);

    H5::DataSet datasetE   = file.openDataSet("E");
    H5::DataSet datasetB   = file.openDataSet("B");
    H5::DataSpace spaceE   = datasetE.getSpace();
    H5::DataSpace spaceB   = datasetB.getSpace();

    hsize_t dimsE[3], dimsB[3];
    spaceE.getSimpleExtentDims(dimsE, nullptr);
    spaceB.getSimpleExtentDims(dimsB, nullptr);

    const size_t n_t  = dimsE[0];
    const size_t n_x  = dimsE[1];
    const size_t n_y2 = dimsE[2];
    const size_t n_y  = dimsB[2];

    struct metadata md = read_metadata(file);

    if (n_y2 / 2 != n_y) {
      throw std::runtime_error("Inconsistent dimensions: n_y2/2 != n_y.");
    }

    size_t nthreads = 
      #if 1
        std::thread::hardware_concurrency() 
      #else
        1
      #endif
      ;
    nthreads = std::min<size_t>(nthreads, n_t);

    const size_t nbatches = n_t / nthreads;

    for (size_t batch_i = 0; batch_i < nbatches; ++batch_i) {

      std::vector<std::thread> threads;
      threads.reserve(nthreads);

      for (size_t slot = 0; slot < nthreads; ++slot) {
        size_t t_index = batch_i * nthreads + slot;

        chunk c;
        c.E.resize(n_x * n_y2);
        c.B.resize(n_x * n_y);
        
        {
          {
            hsize_t offs[3] = { t_index, 0, 0 };
            hsize_t cnt[3]  = { 1, n_x, n_y2 };
            spaceE.selectHyperslab(H5S_SELECT_SET, cnt, offs);
            H5::DataSpace memspaceE(3, cnt);
            datasetE.read(c.E.data(), H5::PredType::NATIVE_DOUBLE,
                           memspaceE, spaceE);
          }

          {
            hsize_t offs[3] = { t_index, 0, 0 };
            hsize_t cnt[3]  = { 1, n_x, n_y };
            spaceB.selectHyperslab(H5S_SELECT_SET, cnt, offs);

            H5::DataSpace memspaceB(3, cnt);
            datasetB.read(c.B.data(), H5::PredType::NATIVE_DOUBLE,
                           memspaceB, spaceB);
          }
        }

        threads.emplace_back([=, md_copy = md]() mutable {

          std::vector<scalar_t> E0(n_x * n_y);
          std::vector<scalar_t> E1(n_x * n_y);
          for (size_t ix = 0; ix < n_x; ++ix) {
            for (size_t iy = 0; iy < n_y; ++iy) {
              E0[ix*n_y + iy] = c.E[ix*n_y2 + 2*iy + 0];
              E1[ix*n_y + iy] = c.E[ix*n_y2 + 2*iy + 1];
            }
          }

          std::string dfname = "/dev/shm/t" + std::to_string(t_index) + ".tmp";

          {
            std::ofstream tmp(dfname, std::ios::out | std::ios::trunc);
            if (!tmp.is_open()) {
              std::cerr << "Failed to open " << dfname << " for writing.\n";
              return;
            }

            auto save_field = [&](const std::vector<scalar_t> &v) {
              for (size_t ix = 0; ix < n_x; ++ix) {
                for (size_t iy = 0; iy < n_y; ++iy) {
                  tmp << v[ix*n_y + iy];
                  if (iy < n_y - 1) tmp << " ";
                }
                tmp << "\n";
              }
              tmp << "\n";
            };

            save_field(c.B);
            save_field(E0);
            save_field(E1);

          }

          // Build a gnuplot command
          std::ostringstream cmd; cmd 
          << "gnuplot -persist << EOF\n"
          << "set terminal pngcairo size 3600,900 enhanced font 'Verdana,20'\n"
          << "set output '" << out_dir << "/t" << t_index << ".png'\n"
          << "set multiplot layout 1,3 title "
          << "'t = " << scalar_t(t_index) * md_copy.dt << " / " << md_copy.tmax
          << ", dx = " << md_copy.dx << ", dy = " << md_copy.dy
          << ", mu = " << md_copy.mu << ", eps = " << md_copy.eps
          << "' font ',36'\n"

          << "set palette defined (-5 'blue', 0 'white', 5 'red')\n"
          << "set cbrange [-5:5]\n"
          << "set zrange [-1:1]\n"

          << "set xrange [0:*]\n"
          << "unset xlabel\n"
          << "set format x ''\n"

          << "set yrange [0:*]\n"
          << "unset ylabel\n"
          << "set format y ''\n"
      
          << "unset colorbox\n"
        
          << "set title 'B'\n"
          << "plot '" << dfname << "' index 0 matrix with image notitle\n"
        
          << "set title 'Ex'\n"
          << "plot '" << dfname << "' index 1 matrix with image notitle\n"
        
          << "set title 'Ey'\n"
          << "plot '" << dfname << "' index 2 matrix with image notitle\n"
        
          << "unset multiplot\n"
          << "EOF\n";

          if (t_index == 0) {
            std::cout << cmd.str() << std::endl;
          }

          int ret = system(cmd.str().c_str());
          if (ret != 0) {
            std::cerr << "[t=" << t_index << "] Gnuplot command failed (exit "
                      << ret << ")\n";
          }

          std::remove(dfname.c_str());

          std::cout << "Generated timestep:" << t_index << "\n";

        });

      }

      for (auto &th : threads) {
        if (th.joinable()) {
          th.join();
        }
      }
      threads.clear();

    }

    file.close();

  } catch (const H5::Exception &e) {
    std::cerr << "HDF5 Error: " << e.getDetailMsg() << std::endl;
    return 1;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
