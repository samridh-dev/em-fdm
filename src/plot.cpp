#include <maxwell.hpp>

#include <H5Cpp.h>

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <algorithm>

const std::string out_dir = "out/";

size_t nthreads = 1;

static struct metadata 
read_metadata(H5::H5File &file) {
  struct metadata md;
  file.openAttribute("tmax").read(H5::PredType::NATIVE_DOUBLE, &md.tmax);
  file.openAttribute("dt").read(H5::PredType::NATIVE_DOUBLE, &md.dt);
  file.openAttribute("xmin").read(H5::PredType::NATIVE_DOUBLE, &md.xmin);
  file.openAttribute("ymin").read(H5::PredType::NATIVE_DOUBLE, &md.ymin);
  file.openAttribute("xmax").read(H5::PredType::NATIVE_DOUBLE, &md.xmax);
  file.openAttribute("ymax").read(H5::PredType::NATIVE_DOUBLE, &md.ymax);
  file.openAttribute("dx").read(H5::PredType::NATIVE_DOUBLE, &md.dx);
  file.openAttribute("dy").read(H5::PredType::NATIVE_DOUBLE, &md.dy);
  file.openAttribute("mu").read(H5::PredType::NATIVE_DOUBLE, &md.mu);
  file.openAttribute("eps").read(H5::PredType::NATIVE_DOUBLE, &md.eps);
  return md;
}

int 
main(int argc, char* argv[]) {

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <file.h5>" << std::endl;
    return 1;
  }

  const std::string fname = argv[1];

  try {

  H5::H5File file(fname, H5F_ACC_RDONLY);

  H5::DataSet datasetE = file.openDataSet("E");
  H5::DataSet datasetB = file.openDataSet("B");

  H5::DataSpace dataspaceE = datasetE.getSpace();
  H5::DataSpace dataspaceB = datasetB.getSpace();

  hsize_t dimsE[3];
  hsize_t dimsB[3];

  dataspaceE.getSimpleExtentDims(dimsE, nullptr);
  dataspaceB.getSimpleExtentDims(dimsB, nullptr);

  const size_t n_t = dimsE[0];
  const size_t n_x = dimsE[1];
  const size_t n_y2 = dimsE[2];
  const size_t n_y = dimsB[2];

  struct metadata md = read_metadata(file);

  if (n_y2 / 2 != n_y) {
    throw std::runtime_error("Inconsistent dimensions between E and B.");
  }

  std::vector<scalar_t> E(n_x * n_y2);
  std::vector<scalar_t> B(n_x * n_y);
  std::vector<scalar_t> E0(n_x * n_y);
  std::vector<scalar_t> E1(n_x * n_y);


  size_t beg = 0;

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  if (nthreads == 0) nthreads = std::thread::hardware_concurrency(); 
  nthreads = std::min<size_t>(nthreads, static_cast<unsigned int>(n_t));
  for (size_t thread_i{0}; thread_i < nthreads; ++thread_i) {

    const size_t chunk_size = n_t / nthreads;
    const size_t remainder = n_t % nthreads;
    const size_t end = beg + chunk_size + (thread_i < remainder ? 1 : 0);

    threads.emplace_back([&, beg, end, thread_i]() {
    for (size_t t_index{beg}; t_index < end; ++t_index) {

      { /* read E */
        hsize_t offs[3] = {t_index, 0, 0};
        hsize_t cnt[3] = {1, n_x, n_y2};
        dataspaceE.selectHyperslab(H5S_SELECT_SET, cnt, offs);
        H5::DataSpace mspace(3, cnt);
        datasetE.read(E.data(), H5::PredType::NATIVE_DOUBLE,
                      mspace, dataspaceE);
      }

      { /* read B */
        hsize_t offs[3] = {t_index, 0, 0};
        hsize_t cnt[3] = {1, n_x, n_y};
        dataspaceB.selectHyperslab(H5S_SELECT_SET, cnt, offs);
        H5::DataSpace mspace(3, cnt);
        datasetB.read(B.data(), H5::PredType::NATIVE_DOUBLE,
                      mspace, dataspaceB);
      }


      { /* read E0, E1 */
        for (size_t i{0}; i < n_x; ++i) {
        for (size_t j{0}; j < n_y; ++j) {
          E0[i * n_y + j] = E[2 * (i * n_y + j) + 0];
          E1[i * n_y + j] = E[2 * (i * n_y + j) + 1];
        }}
      }
      
      std::string dfname = "/dev/shm/t" + std::to_string(t_index) + ".tmp";

      {

        std::ofstream tmp(dfname, std::ios::out | std::ios::trunc);
        if (!tmp.is_open()) {
          std::cerr << "Failed to open tmp.dat for writing.\n";
          return; 
        }

        auto save_field = [&](const std::vector<scalar_t>& v) {
          for (size_t i{0}; i < n_x; ++i) {
            for (size_t j{0}; j < n_y; ++j) {
              tmp << v[i * n_y + j];
              if (j < n_y - 1) tmp << " ";
            } tmp << "\n";
          }   tmp << "\n";
        };

        save_field(B);
        save_field(E0);
        save_field(E1);
        
        tmp.close();
      }


      {
        std::ostringstream cmd; cmd 

    << "gnuplot -persist << EOF\n"
    << "set terminal pngcairo size 3600,900 enhanced font 'Verdana,20'\n"
    << "set output '" << out_dir << "/t" << t_index << ".png'\n"
    << "set multiplot layout 1,3 title "
    << "'t = " << static_cast<scalar_t>(t_index) * md.dt << " / " << md.tmax
    << ", dx = " << md.dx << ", dy = " << md.dy
    << ", mu = " << md.mu << ", eps = " << md.eps
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
          std::cerr << "Gnuplot command failed with exit code " << ret << "\n";
          return;
        }

      }
      std::remove(dfname.c_str());

      std::cout << "Generated timestep: " << t_index << "\r" << std::flush;

    }});

  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
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
