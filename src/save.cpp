#include <maxwell.hpp>
#include <H5Cpp.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

#include <fstream>
#include <sstream>
#include <cstdlib>

void 
maxwell::save(const size_t nt) {

  const std::string fileName = fname + "t" + std::to_string(nt) + ".h5";

  H5::H5File file(fileName, H5F_ACC_TRUNC);

  auto flatten_vfield = [&](const vfield &VF) {
    
    const size_t Nx = VF.size(0);
    const size_t Ny = VF.size(1);

    std::vector<scalar_t> v(Nx * Ny * 2);

    for (size_t i = 0; i < Nx; ++i) {
    for (size_t j = 0; j < Ny; ++j) {
      v[2 * (i * Ny + j) + 0] = VF[Ny * i + j][0];
      v[2 * (i * Ny + j) + 1] = VF[Ny * i + j][1];
    }}

    return v;

  };

  auto flatten_sfield = [&](const sfield &SF) {

    const size_t Nx = SF.size(0);
    const size_t Ny = SF.size(1);

    std::vector<scalar_t> v(Nx * Ny);
    for (size_t i = 0; i < Nx; ++i) {
    for (size_t j = 0; j < Ny; ++j) {
      v[i * Ny + j] = SF[Ny*i + j];
    }}

    return v;

  };

  std::vector<scalar_t> E_data = flatten_vfield(E);
  std::vector<scalar_t> B_data = flatten_sfield(B);

  const size_t Nx = E.size(0);
  const size_t Ny = E.size(1);

  {
    hsize_t dims[2];
    dims[0] = Nx;
    dims[1] = Ny * 2;
    H5::DataSpace dataspace(2, dims);
    H5::DataSet datasetE = file.createDataSet(
      "E", H5::PredType::NATIVE_DOUBLE, dataspace
    );
    datasetE.write(E_data.data(), H5::PredType::NATIVE_DOUBLE);
  }

  {
    hsize_t dims[2];
    dims[0] = Nx;
    dims[1] = Ny;
    H5::DataSpace dataspace(2, dims);
    H5::DataSet datasetB = file.createDataSet(
      "B", H5::PredType::NATIVE_DOUBLE, dataspace
    );
    datasetB.write(B_data.data(), H5::PredType::NATIVE_DOUBLE);
  }

}

void 
h5concat(const std::string &fin, const std::string &fout,
         scalar_t t_max, scalar_t dt, 
         struct metadata md) {

  const size_t n_t = static_cast<size_t>(t_max / dt);
  if (n_t == 0) {
    std::cerr << "No timesteps to concatenate.\n";
    return;
  }

  const std::string first_f_name = fin + "t0.h5";
  H5::H5File first_file(first_f_name, H5F_ACC_RDONLY);
  H5::DataSet dset_e_first = first_file.openDataSet("E");
  H5::DataSpace dspace_e_first = dset_e_first.getSpace();
  hsize_t dims_e[2];
  dspace_e_first.getSimpleExtentDims(dims_e, nullptr);
  const hsize_t n_x = dims_e[0];
  const hsize_t n_y2 = dims_e[1];
  if (n_y2 % 2 != 0) {
    throw std::runtime_error(
      "Dataset E dimension does not match expected (n_y*2)."
    );
  }
  const hsize_t n_y = n_y2 / 2;
  first_file.close();

  const std::string out_f_name = fout + ".h5";
  H5::H5File out_file(out_f_name, H5F_ACC_TRUNC);

  {
    hsize_t dims_e_out[3] = {n_t, n_x, n_y2};
    H5::DataSpace dspace_e_out(3, dims_e_out);
    out_file.createDataSet("E", H5::PredType::NATIVE_DOUBLE, dspace_e_out);
  }

  {
    hsize_t dims_b_out[3] = {n_t, n_x, n_y};
    H5::DataSpace dspace_b_out(3, dims_b_out);
    out_file.createDataSet("B", H5::PredType::NATIVE_DOUBLE, dspace_b_out);
  }

  for (hsize_t t_index = 0; t_index < n_t; ++t_index) {

    const std::string f0 = fin + "t" + std::to_string(t_index) + ".h5";
    H5::H5File in_file(f0, H5F_ACC_RDONLY);

    // ------------------ E dataset ------------------
    {
      H5::DataSet ds_e_in = in_file.openDataSet("E");
      std::vector<double> buf_e(n_x * n_y2);
      ds_e_in.read(buf_e.data(), H5::PredType::NATIVE_DOUBLE);
      H5::DataSet ds_e_out = out_file.openDataSet("E");
      H5::DataSpace fs_e = ds_e_out.getSpace();
      hsize_t offset_e[3] = {t_index, 0, 0};
      hsize_t count_e[3]  = {1, n_x, n_y2};
      fs_e.selectHyperslab(H5S_SELECT_SET, count_e, offset_e);
      H5::DataSpace ms__e(3, count_e);
      ds_e_out.write(buf_e.data(), H5::PredType::NATIVE_DOUBLE, ms__e, fs_e);
    }

    // ------------------ B dataset ------------------
    {
      H5::DataSet ds_b_in = in_file.openDataSet("B");
      std::vector<double> buf_b(n_x * n_y);
      ds_b_in.read(buf_b.data(), H5::PredType::NATIVE_DOUBLE);
      H5::DataSet ds_b_out = out_file.openDataSet("B");
      H5::DataSpace fs_b = ds_b_out.getSpace();
      hsize_t offset_b[3] = {t_index, 0, 0};
      hsize_t count_b[3]  = {1, n_x, n_y};
      fs_b.selectHyperslab(H5S_SELECT_SET, count_b, offset_b);
      H5::DataSpace ms__b(3, count_b);
      ds_b_out.write(buf_b.data(), H5::PredType::NATIVE_DOUBLE, ms__b, fs_b);
    }

    in_file.close();

  }

    try {

    auto add_metadata = [&](const std::string &attr_name, double value) {
      H5::Attribute attr = out_file.createAttribute(
        attr_name, H5::PredType::NATIVE_DOUBLE, H5S_SCALAR
      );
      attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    };

    add_metadata("tmax", md.tmax);
    add_metadata("dt", md.dt);
    add_metadata("xmin", md.xmin);
    add_metadata("ymin", md.ymin);
    add_metadata("xmax", md.xmax);
    add_metadata("ymax", md.ymax);
    add_metadata("dx", md.dx);
    add_metadata("dy", md.dy);
    add_metadata("mu", md.mu);
    add_metadata("eps", md.eps);

  } catch (H5::Exception &error) {
    std::cerr << "Error writing metadata attribute: " 
              << error.getCDetailMsg() << "\n";
  }

  out_file.close();

  for (hsize_t t_index = 0; t_index < n_t; ++t_index) {
    const std::string f0 = fin + "t" + std::to_string(t_index) + ".h5";
    std::remove(f0.c_str());
  }


}

