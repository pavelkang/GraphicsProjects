#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

  void SoftwareRendererImp::draw_svg( SVG& svg ) {
    clear_target();


    supersample_buffer = vector<unsigned char> (4*w*h, 255);
    // set top level transformation
    transformation = canvas_to_screen;

    // draw all elements
    for ( size_t i = 0; i < svg.elements.size(); ++i ) {
      draw_element(svg.elements[i]);
    }

    // draw canvas outline
    Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y++;
    Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y++;
    Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y--;
    Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y--;

    rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
    rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
    rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
    rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

    // resolve and send to render target
    resolve();

  }

  void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {
    this->sample_rate = sample_rate;
    this->w = this->target_w * sample_rate;
    this->h = this->target_h * sample_rate;
    this->supersample_buffer = vector<unsigned char> (4 * w * h, 255);
  }

  void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                               size_t width, size_t height ) {

    // Task 3:
    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;
    // modify supersample_buffer
    this->w = this->target_w * sample_rate;
    this->h = this->target_h * sample_rate;
    this->supersample_buffer = vector<unsigned char> (4 * w * h, 255);
  }

  void SoftwareRendererImp::draw_element( SVGElement* element ) {

    // Task 4 (part 1):
    // Modify this to implement the transformation stack
    Matrix3x3 original = transformation;
    transformation = original * element->transform;
    switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
    }
    transformation = original;
  }


// Primitive Drawing //

  void SoftwareRendererImp::draw_point( Point& point ) {

    Vector2D p = transform(point.position);
    rasterize_point( p.x, p.y, point.style.fillColor );

  }

  void SoftwareRendererImp::draw_line( Line& line ) {

    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

  }

  void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

    Color c = polyline.style.strokeColor;

    if( c.a != 0 ) {
      int nPoints = polyline.points.size();
      for( int i = 0; i < nPoints - 1; i++ ) {
        Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
        Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
        rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
      }
    }
  }

  void SoftwareRendererImp::draw_rect( Rect& rect ) {

    Color c;

    // draw as two triangles
    float x = rect.position.x;
    float y = rect.position.y;
    float w = rect.dimension.x;
    float h = rect.dimension.y;

    Vector2D p0 = transform(Vector2D(   x   ,   y   ));
    Vector2D p1 = transform(Vector2D( x + w ,   y   ));
    Vector2D p2 = transform(Vector2D(   x   , y + h ));
    Vector2D p3 = transform(Vector2D( x + w , y + h ));


    // draw fill
    c = rect.style.fillColor;
    if (c.a != 0 ) {
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
      rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
    }

    // draw outline
    c = rect.style.strokeColor;
    if( c.a != 0 ) {
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
      rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
      rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
      rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
    }

  }

  void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if( c.a != 0 ) {

      // triangulate
      vector<Vector2D> triangles;
      triangulate( polygon, triangles );

      // draw as triangles
      for (size_t i = 0; i < triangles.size(); i += 3) {
        Vector2D p0 = transform(triangles[i + 0]);
        Vector2D p1 = transform(triangles[i + 1]);
        Vector2D p2 = transform(triangles[i + 2]);
        rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
      }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if( c.a != 0 ) {
      int nPoints = polygon.points.size();
      for( int i = 0; i < nPoints; i++ ) {
        Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
        Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
        rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
      }
    }
  }

  void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {
    // Extra credit

  }

  void SoftwareRendererImp::draw_image( Image& image ) {

    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
  }

  void SoftwareRendererImp::draw_group( Group& group ) {

    for ( size_t i = 0; i < group.elements.size(); ++i ) {
      draw_element(group.elements[i]);
    }

  }

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

  float minf(float a, float b) {
    return a > b ? b : a;
  }


  void SoftwareRendererImp::rasterize_point_for_line (float x, float y, Color color) {

    int sx = (int) floor(x) * sample_rate;
    int sy = (int) floor(y) * sample_rate;

    if ( sx < 0 || sx >= w ) return;
    if ( sy < 0 || sy >= h ) return;

    float r = color.r;
    float g = color.g;
    float b = color.b;
    float a = color.a;

    float cr, cg, cb, ca;
    float ca_after, cr_after, cg_after, cb_after;
    int i = sx;
    int j = sy;
    cr = (supersample_buffer[4 * (i + j * w)    ]) / 255.0;
    cg = (supersample_buffer[4 * (i + j * w) + 1]) / 255.0;
    cb = (supersample_buffer[4 * (i + j * w) + 2]) / 255.0;
    ca = (supersample_buffer[4 * (i + j * w) + 3]) / 255.0;

    for (int i=sx; i<sx+sample_rate; i++) {
        for (int j=sy; j<sy+sample_rate; j++) {
            cr = (supersample_buffer[4 * (i + j * w)    ]) / 255.0;
            cg = (supersample_buffer[4 * (i + j * w) + 1]) / 255.0;
            cb = (supersample_buffer[4 * (i + j * w) + 2]) / 255.0;
            ca = (supersample_buffer[4 * (i + j * w) + 3]) / 255.0;

            ca_after = a + ca - a * ca;

            if (ca_after == 0) {
              cr_after = cg_after = cb_after = 0;
            } else {
              cr_after = (r*a + cr * ca * (1.0-a)) / ca_after ;
              cg_after = (g*a + cg * ca * (1.0-a)) / ca_after ;
              cb_after = (b*a + cb * ca * (1.0-a)) / ca_after ;
            }

            supersample_buffer[4 * (i + j * w)    ] = (uint8_t) (cr_after*255);
            supersample_buffer[4 * (i + j * w) + 1] = (uint8_t) (cg_after*255);
            supersample_buffer[4 * (i + j * w) + 2] = (uint8_t) (cb_after*255);
            supersample_buffer[4 * (i + j * w) + 3] = (uint8_t) (ca_after*255);
          }
    }

  }

  void SoftwareRendererImp::super_sample_plot(float x, float y, Color color) {
    // source color values
    float r, g, b, a;
    // canvas color values
    float cr, cg, cb, ca;
    // color values after alpha compositing
    float cr_after, cg_after, cb_after, ca_after;

    int sx = (int) floor(x);
    int sy = (int) floor(y);
    int index = 4 * (sx + sy * w);

    // get source colors
    r = color.r;
    g = color.g;
    b = color.b;
    a = color.a;

    // get canvas colors
    cr = (supersample_buffer[index]) / 255.0;
    cg = (supersample_buffer[index + 1]) / 255.0;
    cb = (supersample_buffer[index + 2]) / 255.0;
    ca = (supersample_buffer[index + 3]) / 255.0;


    // computing alpha composition
    ca_after = a + ca - a * ca;
    if (ca_after == 0) {
      cr_after = cg_after = cb_after = 0;
    } else {
      cr_after = (r * a + cr * ca * (1 - a)) / ca_after;
      cg_after = (g * a + cg * ca * (1 - a)) / ca_after;
      cb_after = (b * a + cb * ca * (1 - a)) / ca_after;
    }

    // write the color values to super render target
    supersample_buffer[index] = (uint8_t)(cr_after * 255);
    supersample_buffer[index + 1] = (uint8_t)(cg_after * 255);
    supersample_buffer[index + 2] = (uint8_t)(cb_after * 255);
    supersample_buffer[index + 3] = (uint8_t)(ca_after * 255);
  }


  void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {
    int i, j, super_x, super_y, index;
    // source color values
    float r, g, b, a;
    // canvas color values
    float cr, cg, cb, ca;
    // color values after alpha compositing
    float cr_after, cg_after, cb_after, ca_after;

    int sx = (int) floor(x);
    int sy = (int) floor(y);

    // check bounds
    if (sx < 0 || sx >= target_w)
      return;
    if (sy < 0 || sy >= target_h)
      return;

    // get source color values
    r = color.r;
    g = color.g;
    b = color.b;
    a = color.a;

    // loop through the super sample grid for each pixel
    for (i = 0; i < sample_rate; i++) {
      for (j = 0; j < sample_rate; j++) {
        super_x = sample_rate * sx + i;
        super_y = sample_rate * sy + j;
        super_sample_plot(super_x, super_y, color);
      }
    }
    /*
    int sx = (int) floor(x * sample_rate);
    int sy = (int) floor(y * sample_rate);

    if ( sx < 0 || sx >= w ) return;
    if ( sy < 0 || sy >= h ) return;

    float ca_after, cr_after, cg_after, cb_after;
    float r = color.r;
    float g = color.g;
    float b = color.b;
    float a = color.a;

    float cr = (supersample_buffer[4 * (sx + sy * w)    ]) / 255.0;
    float cg = (supersample_buffer[4 * (sx + sy * w) + 1]) / 255.0;
    float cb = (supersample_buffer[4 * (sx + sy * w) + 2]) / 255.0;
    float ca = (supersample_buffer[4 * (sx + sy * w) + 3]) / 255.0;

    ca_after = a + ca - a * ca;

    if (ca_after == 0) {
      cr_after = cg_after = cb_after = 0;
    } else {
      cr_after = (r*a + cr * ca * (1.0-a)) / ca_after ;
      cg_after = (g*a + cg * ca * (1.0-a)) / ca_after ;
      cb_after = (b*a + cb * ca * (1.0-a)) / ca_after ;
    }

    this->supersample_buffer[4 * (sx + sy * w)    ] = (uint8_t) (cr_after * 255.0);
    this->supersample_buffer[4 * (sx + sy * w) + 1] = (uint8_t) (cg_after * 255.0);
    this->supersample_buffer[4 * (sx + sy * w) + 2] = (uint8_t) (cb_after * 255.0);
    this->supersample_buffer[4 * (sx + sy * w) + 3] = (uint8_t) (ca_after * 255.0);*/

  }


  void swap(float &x0, float &x1) {
    int temp = x1;
    x1 = x0;
    x0 = temp;
  }

  float min3(float x, float y, float z) {
    return min(min(x, y), z);
  }

  float max3(float x, float y, float z) {
    return max(max(x, y), z);
  }

  /**
   * @brief find nearest sample point(.5) to the "dir"
   * @param x
   * @param dir: 0 for larger, 1 for smaller
   * @return
   **/
  float SoftwareRendererImp::findNearestSamplePoint(float x) {
    return floor(x) + (.5 / this->sample_rate);
  }

  void SoftwareRendererImp::bresenham_line( float x0, float y0,
                                            float x1, float y1,
                                            Color color, int thickness) {
    bool steep = fabs(y1-y0) > fabs(x1-x0);
    if (steep) { // steep slope
      swap(x0, y0);
      swap(x1, y1);
    }
    if (x0 > x1) { // left to right
      swap(x0, x1);
      swap(y0, y1);
    }
    float dx = x1 - x0;
    float dy = fabs(y1-y0);
    float ystep = (y0 < y1) ? 1.0f : -1.0f;
    float y  = y0;
    float eps= 0;
    for (int x=x0; x<=x1; x++) {
      if (steep) {
        for (int i=0; i<thickness; i++) {
          for (int j=0; j<thickness; j++) {
            rasterize_point_for_line(y+i, x+j, color);
          }
        }
      } else {
        for (int i=0; i<thickness; i++) {
          for (int j=0; j<thickness; j++) {
            rasterize_point_for_line(x+i, y+j, color);
          }
        }
      }
      eps += dy;
      if ((eps *2)  >= dx) {
        y += ystep;
        eps -= dx;
      }
    }
  }

  float rfpart(float x) {
    return 1 - (x - floor(x));
  }
  float fpart(float x) {
    return x - floor(x);
  }

  void SoftwareRendererImp::xiaolinwu_line( float x0, float y0,
                                            float x1, float y1,
                                            Color color) {
    bool steep = fabs(y1-y0) > fabs(x1-x0);
    if (steep) { // steep slope
      swap(x0, y0);
      swap(x1, y1);
    }
    if (x0 > x1) { // left to right
      swap(x0, x1);
      swap(y0, y1);
    }
    float dx = x1 - x0;
    float dy = y1 - y0;
    float ystep = dy > 0 ? 1 : -1;
    // vertical line
    if (dx == 0) {
      if (y0 < y1) {
        while (y0 < y1) {
          rasterize_point_for_line(x0, y0, color);
          y0++;
        }
      } else {
        while (y0 > y1) {
          rasterize_point_for_line(x0, y0, color);
          y0--;
        }
      }
      return ;
    }
    float gradient = dy / dx;
    float intery = y0 + gradient*(round(x0+1.0)-x0);
    // first end point
    float xend = round(x0);
    float yend = y0 + gradient*(xend - x0);
    float xgap = rfpart(x0 + 0.5);
    float ypxl1 = floor(yend+1.0);
    if (steep) {

      rasterize_point_for_line(ypxl1, xend, Color(color.r, color.g, color.b,
                                         color.a*(rfpart(yend))*xgap));
      rasterize_point_for_line(ypxl1+1, xend, Color(color.r, color.g, color.b,
                                           color.a*(fpart(yend))*xgap));

    } else {
      rasterize_point_for_line(xend, ypxl1, Color(color.r, color.g, color.b,
                                         color.a*(rfpart(yend)*xgap)));
      rasterize_point_for_line(xend, ypxl1+1, Color(color.r, color.g, color.b,
                                           color.a*(fpart(yend)*xgap)));
    }
    // second end point
    xend = round(x1);
    yend = y1 + gradient*(xend - x1);
    xgap = rfpart(x1 + 0.5);
    ypxl1 = floor(yend+1.0);
    if (steep) {
      rasterize_point_for_line(ypxl1, xend, Color(color.r, color.g, color.b,
                                         color.a*(rfpart(yend))*xgap));
      rasterize_point_for_line(ypxl1+1, xend, Color(color.r, color.g, color.b,
                                           color.a*(fpart(yend))*xgap));
    } else {
      rasterize_point_for_line(xend, ypxl1, Color(color.r, color.g, color.b,
                                         color.a*(rfpart(yend)*xgap)));
      rasterize_point_for_line(xend, ypxl1+1, Color(color.r, color.g, color.b,
                                           color.a*(fpart(yend)*xgap)));
    }
    // main loop
    for (float x=round(x0); x<round(x1); x++) {
      if (steep) {

        rasterize_point_for_line(int(intery), x, Color(color.r, color.g, color.b,
                                         color.a*rfpart(intery)));
        rasterize_point_for_line(int(intery)+1, x, Color(color.r, color.g, color.b,
                                         color.a*(intery - floor(intery))));

      } else {

        rasterize_point_for_line(x, int(intery), Color(color.r, color.g, color.b,
                                           color.a*rfpart(intery)));
        rasterize_point_for_line(x, int(intery)+1, Color(color.r, color.g, color.b,
                                           color.a*(intery-floor(intery))));
      }
      intery += gradient;
    }
  }



  void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                            float x1, float y1,
                                            Color color) {
    // Extra credit 1: support line thickness
    int thickness = 1;
    // Uncomment the following to use Bresenham's algorithm
    bresenham_line(x0, y0, x1, y1, color, thickness);
    // Extra credit 2: xiaolin wu's algorithm
    // I implemented Xiaolin wu's algorithm according to wikipedia, which is off by 1.
    // xiaolinwu_line(x0, y0, x1, y1, color);
  }

  void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                                float x1, float y1,
                                                float x2, float y2,
                                                Color color ) {
    float sample_step = 1.0 / this->sample_rate;
    // find bounding box
    float boxTopLeftX = min3(x0, x1, x2);
    float boxTopLeftY = min3(y0, y1, y2);
    float boxBotRightX= max3(x0, x1, x2);
    float boxBotRightY= max3(y0, y1, y2);
    // find top-left, bot-right sample points
    float sx0 = floor(boxTopLeftX);
    float sy0 = floor(boxTopLeftY);
    float sx1 = ceil(boxBotRightX);
    float sy1 = ceil(boxBotRightY);
    float dx0 = x1 - x0;
    float dy0 = y1 - y0;
    float dx1 = x2 - x1;
    float dy1 = y2 - y1;
    float dx2 = x0 - x2;
    float dy2 = y0 - y2;
    int clockwise = 1;
    if ((dx0 * dy1 -  dx1 * dy0) > 0) {
      clockwise = 1;
    } else {
      clockwise = -1;
    }
    // sampleStep is the distance between two sample points

    // loop and color points
    for (float i=sx0; i<sx1; i+=sample_step) {
      for (float j=sy0; j<sy1; j+=sample_step) {
        if ((((i-x0)*dy0 - (j-y0)*dx0)*clockwise <= 0) &&
            (((i-x1)*dy1 - (j-y1)*dx1)*clockwise <= 0) &&
            (((i-x2)*dy2 - (j-y2)*dx2)*clockwise <= 0)) {
          super_sample_plot(sample_rate*i, sample_rate*j, color);
          //rasterize_point(i, j, color);
        }
      }
    }
  }

  void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                             float x1, float y1,
                                             Texture& tex ) {
    // Implement image rasterization
    float u, v;
    float imageW = x1 - x0;
    float imageH = y1 - y0;
    float zoomFactorU = this->transformation(0, 0);
    float zoomFactorV = this->transformation(1, 1);
    for (float x=floor(x0)+.5; x<floor(x1); x += 1) {
      for (float y=floor(y0)+.5; y<floor(y1); y += 1) {
        u = (x - floor(x0)) / imageW;
        v = (y - floor(y0)) / imageH;
        rasterize_point(x, y, sampler->sample_trilinear(tex, u, v, 1/zoomFactorU, 1/zoomFactorV));
      }
    }
  }

// resolve samples to render target
  void SoftwareRendererImp::resolve( void ) {
    for (int x = 0; x < target_w; x++) {
      for (int y = 0; y < target_h; y++) {
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        float a = 0.0f;
        // Find Box filter color
        for (int dx=0; dx<sample_rate; dx++) {
          for (int dy=0; dy<sample_rate; dy++) {
            int bigX = (x*sample_rate+dx);
            int bigY = (y*sample_rate+dy);
            r += supersample_buffer[4*(bigX+bigY*w)] / 255.0;
            g += supersample_buffer[4*(bigX+bigY*w) + 1] / 255.0;
            b += supersample_buffer[4*(bigX+bigY*w) + 2] / 255.0;
            a += supersample_buffer[4*(bigX+bigY*w) + 3] / 255.0;
          }
        }
        float sq = sample_rate * sample_rate;
        render_target[4 * (x + y * target_w)    ] = (uint8_t) ((r/sq) * 255);
        render_target[4 * (x + y * target_w) + 1] = (uint8_t) ((g/sq) * 255);
        render_target[4 * (x + y * target_w) + 2] = (uint8_t) ((b/sq) * 255);
        render_target[4 * (x + y * target_w) + 3] = (uint8_t) ((a/sq) * 255);
      }
    }
  }
} // namespace CMU462
