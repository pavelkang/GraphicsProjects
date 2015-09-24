#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE(sky):
  // The starter code allocates the mip levels and generates a level
  // map simply fills each level with a color that differs from its
  // neighbours'. The reference solution uses trilinear filtering
  // and it will only work when you have mipmaps.

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors

  int row, col, rLow, rHigh, cLow, cHigh;
  float r, g, b, a;


  for (size_t i=1; i<tex.mipmap.size(); ++i) {
    MipLevel& mip = tex.mipmap[i];
    MipLevel& prevMip = tex.mipmap[i-1];

    for (size_t j = 0; j < 4*mip.width*mip.height; j += 4) {
      rLow = 2 * ((j/4) / mip.width);
      rHigh = 2 * ((j/4) / mip.width + 1);
      cLow = 2 * ((j/4) % mip.width);
      cHigh = 2 * ((j/4) % mip.width + 1);
      r = 0.0;
      g = 0.0;
      b = 0.0;
      a = 0.0;
      for (size_t k = rLow; k < rHigh; k++) {
        for (size_t l = cLow; l < cHigh; l++) {
          r += (float) prevMip.texels[4 * (k * prevMip.width + l)];
          g += (float) prevMip.texels[4 * (k * prevMip.width + l) + 1];
          b += (float) prevMip.texels[4 * (k * prevMip.width + l) + 2];
          a += (float) prevMip.texels[4 * (k * prevMip.width + l) + 3];
        }
      }
      mip.texels[j] = r / 4.0;
      mip.texels[j+1] = g / 4.0;
      mip.texels[j+2] = b / 4.0;
      mip.texels[j+3] = a / 4.0;

    }

  }


  /*
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {
    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }*/

}

Color Sampler2DImp::sample_nearest(Texture& tex,
                                   float u, float v,
                                   int level) {

  // Task ?: Implement nearest neighbour interpolation

  int x = round(u*tex.height);
  int y = round(v*tex.width);
  MipLevel& ml = tex.mipmap[level];
  float r = ml.texels[4*(ml.width*x+y)] / 255.0;
  float g = ml.texels[4*(ml.width*x+y)+1] / 255.0;
  float b = ml.texels[4*(ml.width*x+y)+2] / 255.0;
  float a = ml.texels[4*(ml.width*y+y)+3] / 255.0;
  // return magenta for invalid level
  Color c(r, g, b, a);
  return c;

}

void get_color(Texture &tex, int x, int y, int level, Color& c) {
  MipLevel& ml = tex.mipmap[level];
  c.r = ml.texels[4*(ml.width*y+x)] / 255.0;
  c.g = ml.texels[4*(ml.width*y+x)+1] / 255.0;
  c.b = ml.texels[4*(ml.width*y+x)+2] / 255.0;
  c.a = ml.texels[4*(ml.width*y+x)+3] / 255.0;
}

float getDeci(float x) {
  return x - floor(x);
}
  float getRf(float x) {
    return 1 - (x - floor(x));
  }

Color Sampler2DImp::sample_bilinear(Texture& tex,
                                    float u, float v,
                                    int level) {

  // Task ?: Implement bilinear filtering
  MipLevel& ml = tex.mipmap[level];
  float r = 0.0;
  float g = 0.0;
  float b = 0.0;
  float a = 0.0;
  Color bl, br, tl, tr;
  u = u * ml.width;
  v = v * ml.height;

  get_color(tex, floor(u), floor(v), level, tl);
  get_color(tex, floor(u), floor(v+1), level, bl);
  get_color(tex, round(u+.5), floor(v), level, tr);
  get_color(tex, round(u+.5), floor(v+1), level, br);
  r = tl.r*getRf(u)*getRf(v)+tr.r*getDeci(u)*getRf(v)+bl.r*getRf(u)*getDeci(v)+br.r*getDeci(u)*getDeci(v);
  g = tl.g*getRf(u)*getRf(v)+tr.g*getDeci(u)*getRf(v)+bl.g*getRf(u)*getDeci(v)+br.g*getDeci(u)*getDeci(v);
  b = tl.b*getRf(u)*getRf(v)+tr.b*getDeci(u)*getRf(v)+bl.b*getRf(u)*getDeci(v)+br.b*getDeci(u)*getDeci(v);
  a = tl.a*getRf(u)*getRf(v)+tr.a*getDeci(u)*getRf(v)+bl.a*getRf(u)*getDeci(v)+br.a*getDeci(u)*getDeci(v);
  return Color(r, g, b, a);
}

Color Sampler2DImp::sample_trilinear(Texture& tex,
                                     float u, float v,
                                     float u_scale, float v_scale) {

  float d = log2f(sqrt(v_scale*v_scale + u_scale*u_scale));
  if (d < 0) {
    return sample_bilinear(tex, u, v, 0);
  } else {
    Color lowColor = sample_bilinear(tex, u, v, floor(d));
    Color highColor = sample_bilinear(tex, u, v, floor(d)+1);
    return Color(lowColor.r*getRf(d)+highColor.r*getDeci(d),
                 lowColor.g*getRf(d)+highColor.g*getDeci(d),
                 lowColor.b*getRf(d)+highColor.b*getDeci(d),
                 lowColor.a*getRf(d)+highColor.a*getDeci(d));
  }


}

} // namespace CMU462
