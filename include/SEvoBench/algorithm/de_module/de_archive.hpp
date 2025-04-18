
#include "../../common/population.hpp"
#include "../evolutionary_algorithm.hpp"

namespace sevobench::de_module {
template <std::floating_point T> class de_archive {
protected:
  population<T> archive_;

public:
  virtual std::span<const solution<T>> get() = 0;

  virtual void add(solution<T> &) = 0;

  virtual void resize(const population<T> &) = 0;

  virtual void prepare(const population<T> &) = 0;

  virtual ~de_archive() = default;
};

template <std::floating_point T>
class fifo_archive final : public de_archive<T> {
  int replace_index = 0;
  int A_Size = 0;
  int Archive_Size = 0;
  const T r_Arc = T(2.6);

public:
  fifo_archive() = default;

  fifo_archive(T ratio_archive) : r_Arc(ratio_archive) {}

  void prepare(const population<T> &pop) override {
    replace_index = 0;
    A_Size = 0;
    Archive_Size = static_cast<int>(pop.pop_size() * r_Arc);
    this->archive_.resize(Archive_Size, solution<T>(pop.dim()));
    this->archive_.set_dim(pop.dim());
  }

  void add(solution<T> &sol) override {
    std::swap(this->archive_[replace_index++], sol);
    replace_index %= Archive_Size;
    A_Size = Archive_Size == A_Size ? Archive_Size : A_Size + 1;
  }

  void resize(const population<T> &pop) override {
    if (Archive_Size > static_cast<int>(r_Arc * pop.pop_size())) {
      Archive_Size = static_cast<int>(r_Arc * pop.pop_size());
      if (Archive_Size < A_Size) {
        int len = A_Size - Archive_Size;
        if ((replace_index < A_Size) && replace_index != 0) {
          std::rotate(this->archive_.begin(),
                      this->archive_.begin() + replace_index,
                      this->archive_.begin() + A_Size);
        }
        this->archive_.erase(this->archive_.begin(),
                             this->archive_.begin() + len);
        A_Size = Archive_Size;
        replace_index = 0;
      } else if (Archive_Size == A_Size) {
        replace_index %= Archive_Size;
      }
    }
  }

  std::span<const solution<T>> get() override {
    return std::span<const solution<T>>(this->archive_.data(), A_Size);
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class random_archive final : public de_archive<T> {
  int replace_index = 0;
  int A_Size = 0;
  int Archive_Size = 0;
  const T r_Arc = T(2.6);

public:
  sevobench::tool::rng RNG;

  random_archive() = default;

  random_archive(T ratio_archive) : r_Arc(ratio_archive) {}

  void prepare(const population<T> &pop) override {
    replace_index = 0;
    A_Size = 0;
    Archive_Size = static_cast<int>(pop.pop_size() * r_Arc);
    this->archive_.resize(Archive_Size + pop.pop_size(),
                          solution<T>(pop.dim()));
    this->archive_.set_dim(pop.dim());
  }

  void add(solution<T> &sol) override {
    std::swap(this->archive_[replace_index++], sol);
    A_Size++;
  }

  void resize(const population<T> &pop) override {
    Archive_Size = static_cast<int>(r_Arc * pop.pop_size());
    if (Archive_Size < A_Size) {
      int len = A_Size - Archive_Size;
      for (int i = 0; i < len; i++) {
        auto i1 = RNG.rand_int(A_Size--);
        std::swap(this->archive_[i1], this->archive_[A_Size]);
      }
      replace_index = Archive_Size;
    }
  }

  std::span<const solution<T>> get() override {
    return std::span<const solution<T>>(this->archive_.data(), A_Size);
  }
};
} // namespace sevobench::de_module
