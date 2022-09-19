#ifndef AOCOMMON_IMAGE_ACCESSOR_H_
#define AOCOMMON_IMAGE_ACCESSOR_H_

#include "image.h"

namespace aocommon {

/**
 * @brief Abstract interface for loading and storing an image.
 *
 * An ImageAccessor object knows how to store a single image and load it back
 * later. When used in an interface, the caller can implement storing the image
 * in different ways (e.g. on disk, in memory, etc)."
 */
class ImageAccessor {
 public:
  virtual ~ImageAccessor() {}

  /**
   * @brief Load the image.
   *
   * @param image Location where the image will be stored.
   */
  virtual void Load(Image& image) const = 0;

  /**
   * @brief Store the image, so it can be loaded back later.
   *
   * @param image The image that must be stored.
   */
  virtual void Store(const Image& image) = 0;
};

}  // namespace aocommon

#endif