#ifndef EQUILIBRIUM_SRC_HDF5_OUTPUT_HPP_
#define EQUILIBRIUM_SRC_HDF5_OUTPUT_HPP_

#include "equilibrium/writer_helpers.hpp"

#include <cstddef>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <hdf5.h>

// Also covers callers of equilibrium_add_io_library with supplied targets.
#if !H5_VERSION_GE(1, 14, 6)
#error "equilibrium-io requires HDF5 >= 1.14.6"
#endif

namespace equilibrium::io_detail {

template <herr_t (*Close)(hid_t)>
class Hdf5Handle {
   public:
    explicit Hdf5Handle(hid_t id) : id_(id) {}
    ~Hdf5Handle() {
        if (id_ >= 0) Close(id_);
    }
    Hdf5Handle(const Hdf5Handle&) = delete;
    Hdf5Handle& operator=(const Hdf5Handle&) = delete;
    hid_t get() const { return id_; }
    herr_t close() { return Close(std::exchange(id_, -1)); }

   private:
    hid_t id_;
};

// HDF5 can leave its library state damaged after a filesystem write failure.
// Build an image without a backing file, close HDF5, then publish through the
// ordinary checked file writer. Extraction briefly holds two image buffers.
class Hdf5Output {
   public:
    explicit Hdf5Output(const std::string& path)
        : path_(path),
          temporary_(temp_path_for(path)),
          file_(create_image(temporary_)) {}
    ~Hdf5Output() {
        output_.reset();
        if (owns_temporary_) std::remove(temporary_.c_str());
    }
    hid_t get() const { return file_.get(); }

    std::string publish() {
        if (H5Fflush(file_.get(), H5F_SCOPE_LOCAL) < 0)
            return "H5Fflush failed";
        const auto size = H5Fget_file_image(file_.get(), nullptr, 0);
        if (size <= 0) return "H5Fget_file_image size failed";
        std::vector<std::byte> image(static_cast<std::size_t>(size));
        // Use the public image API, which also clears superblock status flags.
        if (H5Fget_file_image(file_.get(), image.data(), image.size()) != size)
            return "H5Fget_file_image failed";
        if (file_.close() < 0) return "H5Fclose failed";

        output_.reset(std::fopen(temporary_.c_str(), "wbx"));
        if (!output_) return "cannot create HDF5 output temporary file";
        owns_temporary_ = true;
        if (std::fwrite(image.data(), 1, image.size(), output_.get()) !=
            image.size())
            return "could not write HDF5 output image";
        const auto error = publish_atomic(output_.release(), temporary_, path_);
        owns_temporary_ = false;
        return error;
    }

   private:
    static hid_t create_image(const std::string& name) {
        Hdf5Handle<H5Pclose> access(H5Pcreate(H5P_FILE_ACCESS));
        if (access.get() < 0 ||
            H5Pset_fapl_core(access.get(), 64 * 1024, false) < 0)
            return -1;
        return H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         access.get());
    }

    std::string path_, temporary_;
    Hdf5Handle<H5Fclose> file_;
    std::unique_ptr<std::FILE, decltype(&std::fclose)> output_{nullptr,
                                                               std::fclose};
    bool owns_temporary_ = false;
};

}  // namespace equilibrium::io_detail

#endif
