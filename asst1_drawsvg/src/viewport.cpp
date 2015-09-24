#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 4 (part 2):
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvas coordinates.
  this->x = x;
  this->y = y;
  this->span = span;
  float scaleBy = 1.0 / (span * 2);
  float translateX = x - span;
  float translateY = y - span;
  Matrix3x3 scaleTransform;
  scaleTransform(0, 0) = scaleBy;
  scaleTransform(0, 1) = 0;
  scaleTransform(0, 2) = 0;
  scaleTransform(1, 0) = 0;
  scaleTransform(1, 1) = scaleBy;
  scaleTransform(1, 2) = 0;
  scaleTransform(2, 0) = 0;
  scaleTransform(2, 1) = 0;
  scaleTransform(2, 2) = 1;
  Matrix3x3 translateTransform;
  translateTransform(0, 0) = 1;
  translateTransform(0, 1) = 0;
  translateTransform(0, 2) = -translateX;
  translateTransform(1, 0) = 0;
  translateTransform(1, 1) = 1;
  translateTransform(1, 2) = -translateY;
  translateTransform(2, 0) = 0;
  translateTransform(2, 1) = 0;
  translateTransform(2, 2) = 1;
  this->svg_2_norm = scaleTransform * translateTransform;
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) {

  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CMU462
