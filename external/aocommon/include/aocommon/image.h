#ifndef AOCOMMON_IMAGE_H_
#define AOCOMMON_IMAGE_H_

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>
#include <aocommon/uvector.h>

#include <cassert>
#include <cstring>
#include <cmath>
#include <complex>
#include <memory>

namespace aocommon {
namespace detail {
/**
 * Because Image can be instantiated for complex values, a
 * helper function is used to implement ordering operations like
 * median that do not have complex implementations.
 */
template <typename NumT>
inline NumT MedianWithCopyImplementation(const NumT* data, size_t size,
                                         aocommon::UVector<NumT>& copy) {
  copy.reserve(size);
  for (const NumT* i = data; i != data + size; ++i) {
    if (std::isfinite(*i)) copy.push_back(*i);
  }
  if (copy.empty())
    return 0.0;
  else {
    bool even = (copy.size() % 2) == 0;
    typename aocommon::UVector<NumT>::iterator mid =
        copy.begin() + (copy.size() - 1) / 2;
    std::nth_element(copy.begin(), mid, copy.end());
    NumT median = *mid;
    if (even) {
      std::nth_element(mid, mid + 1, copy.end());
      median = (median + *(mid + 1)) * 0.5;
    }
    return median;
  }
}

template <typename NumT>
inline NumT MadImplementation(const NumT* data, size_t size) {
  aocommon::UVector<NumT> copy;
  NumT median = MedianWithCopyImplementation(data, size, copy);
  if (copy.empty()) return 0.0;

  // Replace all values by the difference from the mean
  typename aocommon::UVector<NumT>::iterator mid =
      copy.begin() + (copy.size() - 1) / 2;
  for (typename aocommon::UVector<NumT>::iterator i = copy.begin();
       i != mid + 1; ++i)
    *i = median - *i;
  for (typename aocommon::UVector<NumT>::iterator i = mid + 1; i != copy.end();
       ++i)
    *i = *i - median;

  std::nth_element(copy.begin(), mid, copy.end());
  median = *mid;
  bool even = (copy.size() % 2) == 0;
  if (even) {
    std::nth_element(mid, mid + 1, copy.end());
    median = (median + *(mid + 1)) * 0.5;
  }
  return median;
}
}  // namespace detail

template <typename NumT>
class ImageBase {
 public:
  using value_type = NumT;
  using num_t = value_type;

  using iterator = value_type*;
  using const_iterator = const value_type*;

  constexpr ImageBase() noexcept : data_(nullptr), width_(0), height_(0) {}

  /**
   * @brief Construct new ImageBase object of given width and
   * height, but with uninitialized values.
   */
  ImageBase(size_t width, size_t height)
      : data_(new value_type[width * height]), width_(width), height_(height) {}

  /**
   * @brief Construct an new ImageBase object with given width and height, with
   * all values initialized to \c initial_value.
   */
  ImageBase(size_t width, size_t height, value_type initial_value)
      : data_(new value_type[width * height]), width_(width), height_(height) {
    std::fill(data_, data_ + width_ * height_, initial_value);
  }

  /**
   * @brief Construct an new ImageBase object with given width and height, with
   * values initialized from an initializer list.
   */
  ImageBase(size_t width, size_t height,
            std::initializer_list<NumT> initial_values)
      : data_(new value_type[width * height]), width_(width), height_(height) {
    assert(initial_values.size() == width * height);
    std::copy_n(initial_values.begin(), width_ * height_, data_);
  }

  ~ImageBase() noexcept { delete[] data_; }

  ImageBase(const ImageBase<NumT>& source)
      : data_(new value_type[source.width_ * source.height_]),
        width_(source.width_),
        height_(source.height_) {
    std::copy(source.data_, source.data_ + width_ * height_, data_);
  }

  ImageBase<NumT>& operator=(const ImageBase<NumT>& source) {
    if (width_ * height_ != source.width_ * source.height_) {
      delete[] data_;
      // Make data_ robust against failure of "new" (e.g. due to out of mem
      // failures)
      data_ = nullptr;
      data_ = new value_type[source.width_ * source.height_];
    }
    width_ = source.width_;
    height_ = source.height_;
    std::copy(source.data_, source.data_ + width_ * height_, data_);
    return *this;
  }

  ImageBase<NumT>& operator=(value_type value) noexcept {
    std::fill_n(data_, width_ * height_, value);
    return *this;
  }

  ImageBase<NumT>& Assign(const NumT* begin, const NumT* end) {
    if (size_t(end - begin) != width_ * height_)
      throw std::runtime_error("Invalid Assign()");
    std::copy(begin, end, data_);
    return *this;
  }

  ImageBase(ImageBase<NumT>&& source) noexcept
      : data_(source.data_), width_(source.width_), height_(source.height_) {
    source.width_ = 0;
    source.height_ = 0;
    source.data_ = nullptr;
  }

  ImageBase<NumT>& operator=(ImageBase<NumT>&& source) noexcept {
    std::swap(data_, source.data_);
    std::swap(width_, source.width_);
    std::swap(height_, source.height_);
    return *this;
  }

  static std::unique_ptr<ImageBase<NumT>> Make(size_t width, size_t height) {
    return std::unique_ptr<ImageBase<NumT>>(new ImageBase<NumT>(width, height));
  }
  static std::unique_ptr<ImageBase<NumT>> Make(size_t width, size_t height,
                                               value_type initial_value) {
    return std::unique_ptr<ImageBase<NumT>>(
        new ImageBase<NumT>(width, height, initial_value));
  }

  value_type* Data() { return data_; }
  const value_type* Data() const { return data_; }

  size_t Width() const { return width_; }
  size_t Height() const { return height_; }
  size_t Size() const { return width_ * height_; }
  bool Empty() const { return width_ == 0 || height_ == 0; }

  iterator begin() { return data_; }
  const_iterator begin() const { return data_; }

  iterator end() { return data_ + width_ * height_; }
  const_iterator end() const { return data_ + width_ * height_; }

  const value_type& operator[](size_t index) const { return data_[index]; }
  value_type& operator[](size_t index) { return data_[index]; }

  ImageBase<NumT>& operator+=(const ImageBase<NumT>& other) {
    assert(Size() == other.Size());
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] += other[i];
    return *this;
  }

  ImageBase<NumT>& operator-=(const ImageBase<NumT>& other) {
    assert(Size() == other.Size());
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] -= other[i];
    return *this;
  }

  ImageBase<NumT>& operator*=(value_type factor) {
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] *= factor;
    return *this;
  }

  ImageBase<NumT>& operator*=(const ImageBase<NumT>& other) {
    assert(Size() == other.Size());
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] *= other[i];
    return *this;
  }

  ImageBase<NumT>& operator/=(value_type factor) {
    return (*this) *= value_type(1.0) / factor;
  }

  ImageBase<NumT>& Sqrt() {
    for (size_t i = 0; i != width_ * height_; ++i)
      data_[i] = std::sqrt(data_[i]);
    return *this;
  }

  ImageBase<NumT>& SqrtWithFactor(NumT factor) {
    for (size_t i = 0; i != width_ * height_; ++i)
      data_[i] = std::sqrt(data_[i]) * factor;
    return *this;
  }

  ImageBase<NumT>& Square() {
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] *= data_[i];
    return *this;
  }

  ImageBase<NumT>& SquareWithFactor(NumT factor) {
    for (size_t i = 0; i != width_ * height_; ++i)
      data_[i] *= data_[i] * factor;
    return *this;
  }

  ImageBase<NumT>& AddWithFactor(const ImageBase<NumT>& rhs, NumT factor) {
    assert(Size() == rhs.Size());
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] += rhs[i] * factor;
    return *this;
  }

  ImageBase<NumT>& AddSquared(const ImageBase<NumT>& rhs) {
    assert(Size() == rhs.Size());
    for (size_t i = 0; i != width_ * height_; ++i) data_[i] += rhs[i] * rhs[i];
    return *this;
  }

  ImageBase<NumT>& AddSquared(const ImageBase<NumT>& rhs, NumT factor) {
    assert(Size() == rhs.Size());
    for (size_t i = 0; i != width_ * height_; ++i)
      data_[i] += rhs[i] * rhs[i] * factor;
    return *this;
  }

  void Reset() {
    delete[] data_;
    data_ = nullptr;
    width_ = 0;
    height_ = 0;
  }

  /** Cut-off the borders of an image.
   * @param out_width Should be &lt;= in_width.
   * @param out_height Should be &lt;= in_height.
   */
  static void Trim(value_type* output, size_t out_width, size_t out_height,
                   const value_type* input, size_t in_width, size_t in_height) {
    const size_t start_x = (in_width - out_width) / 2;
    const size_t start_y = (in_height - out_height) / 2;
    const size_t end_y = (in_height + out_height) / 2;
    for (size_t y = start_y; y != end_y; ++y) {
      std::copy_n(&input[y * in_width + start_x], out_width,
                  &output[(y - start_y) * out_width]);
    }
  }

  ImageBase<NumT> Trim(size_t out_width, size_t out_height) const {
    ImageBase<NumT> image(out_width, out_height);
    Trim(image.Data(), out_width, out_height, Data(), Width(), Height());
    return image;
  }

  /**
   * @brief Trim image to box, given a lower left (x,y)-coordinate
   * and box width and height.
   *
   * @param x1 Lower left x-coordinate
   * @param y1 Lower left y-coordinate
   * @param box_width Box width
   * @param box_height Box height
   * @return ImageBase<NumT> Trimmed image
   */
  ImageBase<NumT> TrimBox(size_t x1, size_t y1, size_t box_width,
                          size_t box_height) const {
    ImageBase<NumT> image(box_width, box_height);
    TrimBox(image.Data(), x1, y1, box_width, box_height, Data(), Width(),
            Height());
    return image;
  }

  template <typename T>
  static void TrimBox(T* output, size_t x1, size_t y1, size_t box_width,
                      size_t box_height, const T* input, size_t in_width,
                      size_t in_height) {
    size_t end_y = y1 + box_height;
    for (size_t y = y1; y != end_y; ++y) {
      std::copy_n(&input[y * in_width + x1], box_width,
                  &output[(y - y1) * box_width]);
    }
  }

  template <typename T>
  static void CopyMasked(T* to, size_t to_x, size_t to_y, size_t to_width,
                         const T* from, size_t from_width, size_t from_height,
                         const bool* from_mask) {
    for (size_t y = 0; y != from_height; ++y) {
      for (size_t x = 0; x != from_width; ++x) {
        if (from_mask[y * from_width + x])
          to[to_x + (to_y + y) * to_width + x] = from[y * from_width + x];
      }
    }
  }

  /**
   * Add the contents of an image to a bigger image at a specific place.
   * This function is templated such that it also works with images stored
   * as (one-dimensional) arrays of size width x height.
   * The height of the target image is not necessary as parameter, but it
   * is assumed that the subimage fits completely in the target image at the
   * specified position.
   *
   * @param to Target image to which the subimage is added.
   * @param to_x Horizontal size of the target image where the subimage is
   * placed
   * @param to_y Vertical size of target image
   * @param to_width Width of target image.
   * @param from Subimage to be placed on the target image
   * @param from_width Horizontal size of the @c from image
   * @param from_height Vertical size of the @c from image
   */
  template <typename ValueTypeTo, typename ValueTypeFrom>
  static void AddSubImage(ValueTypeTo* to, size_t to_x, size_t to_y,
                          size_t to_width, const ValueTypeFrom* from,
                          size_t from_width, size_t from_height) {
    assert(to_x + from_width <= to_width);
    for (size_t y = 0; y != from_height; ++y) {
      ValueTypeTo* to_line = &to[to_x + (to_y + y) * to_width];
      const ValueTypeFrom* from_line = &from[y * from_width];
      for (size_t x = 0; x != from_width; ++x) {
        to_line[x] += from_line[x];
      }
    }
  }

  /** Extend an image with zeros, complement of Trim.
   * @param out_width Should be &gt;= in_width.
   * @param out_height Should be &gt;= in_height.
   */
  static void Untrim(value_type* output, size_t out_width, size_t out_height,
                     const value_type* input, size_t in_width,
                     size_t in_height) {
    const size_t start_x = (out_width - in_width) / 2;
    const size_t end_x = (out_width + in_width) / 2;
    const size_t start_y = (out_height - in_height) / 2;
    const size_t end_y = (out_height + in_height) / 2;
    for (size_t y = 0; y != start_y; ++y) {
      value_type* ptr = &output[y * out_width];
      std::fill_n(ptr, out_width, 0.0);
    }
    for (size_t y = start_y; y != end_y; ++y) {
      value_type* ptr = &output[y * out_width];
      std::fill_n(ptr, start_x, 0.0);
      std::copy_n(&input[(y - start_y) * in_width], in_width,
                  &output[y * out_width + start_x]);
      std::fill_n(ptr + end_x, out_width, 0.0);
    }
    for (size_t y = end_y; y != out_height; ++y) {
      value_type* ptr = &output[y * out_width];
      std::fill_n(ptr, out_width, 0.0);
    }
  }

  ImageBase<NumT> Untrim(size_t out_width, size_t out_height) const {
    ImageBase<NumT> image(out_width, out_height);
    Untrim(image.Data(), out_width, out_height, Data(), Width(), Height());
    return image;
  }

  static value_type Median(const value_type* data, size_t size);

  static value_type MAD(const value_type* data, size_t size);

  value_type Sum() const {
    value_type sum = 0.0;
    for (const value_type& v : *this) sum += v;
    return sum;
  }

  value_type Average() const { return Sum() / NumT(Size()); }

  value_type Min() const;
  value_type Max() const;

  value_type StdDevFromMAD() const {
    return StdDevFromMAD(data_, width_ * height_);
  }
  static value_type StdDevFromMAD(const value_type* data, size_t size) {
    // norminv(0.75) x MAD
    return value_type(1.48260221850560) * MAD(data, size);
  }

  value_type RMS() const { return RMS(data_, width_ * height_); }

  static value_type RMS(const value_type* data, size_t size) {
    value_type sum = 0.0;
    for (size_t i = 0; i != size; ++i) sum += data[i] * data[i];
    return std::sqrt(sum / value_type(size));
  }

  void Negate() {
    for (value_type& d : *this) d = -d;
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt64(width_).UInt64(height_);
    size_t n = sizeof(NumT) * width_ * height_;
    std::copy_n(reinterpret_cast<const unsigned char*>(data_), n,
                stream.Chunk(n));
  }
  void Unserialize(aocommon::SerialIStream& stream) {
    delete[] data_;
    // Make robust against stream failing
    data_ = nullptr;
    stream.UInt64(width_).UInt64(height_);
    if (width_ * height_ != 0) {
      data_ = new value_type[width_ * height_];
    }
    const size_t n = sizeof(NumT) * width_ * height_;
    std::copy_n(stream.Chunk(n), n, reinterpret_cast<unsigned char*>(data_));
  }

 private:
  value_type* data_;
  size_t width_;
  size_t height_;

  static value_type MedianWithCopy(const value_type* data, size_t size,
                                   aocommon::UVector<value_type>& copy);
};

/**
 * A single-precision two-dimensional image.
 */
using Image = ImageBase<float>;

/**
 * A double-precision two-dimensional image.
 * This type is (currently) not used inside WSClean, but it is used in
 * some other projects that share the same Image class, hence it is here.
 */
using DImage = ImageBase<double>;

using ImageCF = ImageBase<std::complex<float>>;

template <class T>
using ComplexImageBase = ImageBase<std::complex<T>>;

template <>
inline typename ImageBase<double>::value_type ImageBase<double>::Min() const {
  return *std::min_element(begin(), end());
}

template <>
inline typename ImageBase<float>::value_type ImageBase<float>::Min() const {
  return *std::min_element(begin(), end());
}

template <>
inline typename ImageBase<double>::value_type ImageBase<double>::Max() const {
  return *std::max_element(begin(), end());
}

template <>
inline typename ImageBase<float>::value_type ImageBase<float>::Max() const {
  return *std::max_element(begin(), end());
}

template <>
inline typename ImageBase<float>::value_type ImageBase<float>::MedianWithCopy(
    const value_type* data, size_t size, aocommon::UVector<value_type>& copy) {
  return detail::MedianWithCopyImplementation(data, size, copy);
}

template <>
inline typename ImageBase<double>::value_type ImageBase<double>::MedianWithCopy(
    const value_type* data, size_t size, aocommon::UVector<value_type>& copy) {
  return detail::MedianWithCopyImplementation(data, size, copy);
}

template <typename NumT>
typename ImageBase<NumT>::value_type ImageBase<NumT>::MedianWithCopy(
    const value_type*, size_t, aocommon::UVector<value_type>&) {
  throw std::runtime_error("not implemented");
}

template <typename NumT>
typename ImageBase<NumT>::value_type ImageBase<NumT>::Median(
    const value_type* data, size_t size) {
  aocommon::UVector<value_type> copy;
  return MedianWithCopy(data, size, copy);
}

template <>
inline typename ImageBase<float>::value_type ImageBase<float>::MAD(
    const value_type* data, size_t size) {
  return detail::MadImplementation(data, size);
}

template <>
inline typename ImageBase<double>::value_type ImageBase<double>::MAD(
    const value_type* data, size_t size) {
  return detail::MadImplementation(data, size);
}

template <typename NumT>
typename ImageBase<NumT>::value_type ImageBase<NumT>::MAD(const value_type*,
                                                          size_t) {
  throw std::runtime_error("not implemented");
}

}  // namespace aocommon
#endif
