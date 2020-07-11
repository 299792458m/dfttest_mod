/*
**                    dfttest for Avisynth+
**
**   2D/3D frequency domain denoiser.
**
**   Copyright (C) 2007-2010 Kevin Stone, 2017 (C) DJATOM
**             (C) 2020 pinterf
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "dfttest_avx2.h"
#include "avs/config.h"
#include <stdint.h>
#include "fmath.h"
#include <windows.h>

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
// 8 bit data -> float, 8 pixels at a time
void proc0_uint8_to_float_AVX2_8pixels(const unsigned char* s0, const float* s1, float* d,
  const int p0, const int p1, const int /*offset_lsb*/)
{
  for (int u = 0; u < p1; ++u)
  {
    for (int v = 0; v < p1; v += 8)
    {
      auto src = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s0 + p0 * u + v));
      // 8x 8bit->8x32 bits
      auto int8x32 = _mm256_cvtepu8_epi32(src);
      // int -> float
      auto float8x32 = _mm256_cvtepi32_ps(int8x32);
      // mul by weight
      auto w = _mm256_loadu_ps(s1 + p1 * u + v);
      auto result = _mm256_mul_ps(float8x32, w);
      _mm256_storeu_ps(d + p1 * u + v, result);
    }
  }
  _mm256_zeroupper();
}

template<int bits_per_pixel>
#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void proc0_uint16_to_float_AVX2_8pixels(const unsigned char* s0, const float* s1, float* d,
  const int p0, const int p1, const int /*offset_lsb*/)
{
  constexpr float scaler1 = 1.0f / (1 << (bits_per_pixel - 8)); // back to the 0.0 - 255.0 range
  auto scaler = _mm256_set1_ps(scaler1);
  for (int u = 0; u < p1; ++u)
  {
    for (int v = 0; v < p1; v += 8)
    {
      // 8x16
      auto src = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s0 + v * sizeof(uint16_t)));
      //8x16->8x32
      auto int8x32 = _mm256_cvtepu16_epi32(src);
      // int->float
      auto float8x32 = _mm256_cvtepi32_ps(int8x32);
      // weight
      auto w = _mm256_loadu_ps(s1 + v);
      auto result = _mm256_mul_ps(float8x32, w);
      // scale back
      result = _mm256_mul_ps(result, scaler);
      _mm256_storeu_ps(d + v, result);
    }
    s0 += p0;
    s1 += p1;
    d += p1;
  }
  _mm256_zeroupper();
}

// instantiate
template void proc0_uint16_to_float_AVX2_8pixels<10>(const unsigned char* s0, const float* s1, float* d, const int p0, const int p1, const int /*offset_lsb*/);
template void proc0_uint16_to_float_AVX2_8pixels<12>(const unsigned char* s0, const float* s1, float* d, const int p0, const int p1, const int /*offset_lsb*/);
template void proc0_uint16_to_float_AVX2_8pixels<14>(const unsigned char* s0, const float* s1, float* d, const int p0, const int p1, const int /*offset_lsb*/);
template void proc0_uint16_to_float_AVX2_8pixels<16>(const unsigned char* s0, const float* s1, float* d, const int p0, const int p1, const int /*offset_lsb*/);

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void proc0_float_to_float_AVX2_8pixels(const unsigned char* s0, const float* s1, float* d,
  const int p0, const int p1, const int /*offset_lsb*/)
{
  constexpr float scaler1 = 255.0f; // from 0..1 =-0.5..+0.5 back to the *256 range
  auto scaler = _mm256_set1_ps(scaler1);
  for (int u = 0; u < p1; ++u)
  {
    for (int v = 0; v < p1; v += 8)
    {
      // 8x float
      auto src = _mm256_loadu_ps(reinterpret_cast<const float *>(s0 + p0 * u + v * sizeof(float)));
      // weight
      auto w = _mm256_loadu_ps(s1 + p1 * u + v);
      auto result = _mm256_mul_ps(src, w);
      // scale back
      result = _mm256_mul_ps(result, scaler);
      _mm256_storeu_ps(d + p1 * u + v, result);
    }
  }
  _mm256_zeroupper();
}


#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void proc1_AVX2_8(const float* s0, const float* s1, float* d,
  const int p0, const int p1)
{
  for (int u = 0; u < p0; ++u)
  {
    for (int v = 0; v < p0; v += 8)
    {
      auto s0_loop = _mm256_loadu_ps(s0 + u * p0 + v);
      auto s1_loop = _mm256_loadu_ps(s1 + u * p0 + v);
      auto d_loop = _mm256_loadu_ps(d + u * p1 + v);
      auto result = _mm256_fmadd_ps(s0_loop, s1_loop, d_loop);
      _mm256_storeu_ps(d + u * p1 + v, result);
    }
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_0_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  auto zero = _mm256_setzero_ps();
  auto sse_1em15 = _mm256_set1_ps(1e-15f);
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq

    auto diff = _mm256_sub_ps(sumofsquares, sigmas_loop); // SQ-S SQ-S sq-s sq-s // . is undefined
    auto div_as_rcp_of_ss = _mm256_rcp_ps(_mm256_add_ps(sumofsquares, sse_1em15)); // 1 / (sq+1e-15)

    auto coeff = _mm256_max_ps(zero, _mm256_mul_ps(diff, div_as_rcp_of_ss));
    // max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f);

    auto result = _mm256_mul_ps(dftc_loop, coeff);
    _mm256_storeu_ps(dftc + h, result);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_1_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq

    auto cmple1_loop = _mm256_cmp_ps(sigmas_loop, sumofsquares, _CMP_LE_OQ);
    auto and1_loop = _mm256_and_ps(cmple1_loop, dftc_loop);
    _mm256_storeu_ps(dftc + h, and1_loop);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_2_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h);
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h);
    auto mul_loop = _mm256_mul_ps(dftc_loop, sigmas_loop);
    _mm256_storeu_ps(dftc + h, mul_loop);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_3_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto sigmas2_loop = _mm256_loadu_ps(sigmas2 + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq

    auto pmin_loop = _mm256_loadu_ps(pmin + h);
    auto pmax_loop = _mm256_loadu_ps(pmax + h);

    //if (psd >= pmin[h] && psd <= pmax[h])
    auto cmp1_loop = _mm256_cmp_ps(pmin_loop, sumofsquares, _CMP_LE_OQ);
    auto cmp2_loop = _mm256_cmp_ps(sumofsquares, pmax_loop, _CMP_LE_OQ);
    auto and1_loop = _mm256_and_ps(cmp1_loop, cmp2_loop);
    auto sigma_chosen = _mm256_blendv_ps(sigmas2_loop, sigmas_loop, and1_loop);
    auto mul2_loop = _mm256_mul_ps(sigma_chosen, dftc_loop);
    _mm256_storeu_ps(dftc + h, mul2_loop);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_4_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  auto sse_1em15 = _mm256_set1_ps(1e-15f);
  for (int h = 0; h < ccnt; h += 4)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq
    auto psd = _mm256_add_ps(sse_1em15, sumofsquares);

    auto pmin_loop = _mm256_loadu_ps(pmin + h);
    auto pmax_loop = _mm256_loadu_ps(pmax + h);

    auto add3_loop = _mm256_add_ps(psd, pmin_loop);
    auto add4_loop = _mm256_add_ps(psd, pmax_loop);
    auto mul2_loop = _mm256_mul_ps(psd, pmax_loop);
    auto mul3_loop = _mm256_mul_ps(add4_loop, add3_loop);
    auto rcp1_loop = _mm256_rcp_ps(mul3_loop);
    auto mul4_loop = _mm256_mul_ps(rcp1_loop, mul2_loop);
    auto sqrt_loop = _mm256_sqrt_ps(mul4_loop);
    auto mul5_loop = _mm256_mul_ps(sigmas_loop, sqrt_loop);
    auto mul6_loop = _mm256_mul_ps(mul5_loop, dftc_loop);
    _mm256_storeu_ps(dftc + h, mul6_loop);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_5_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  auto zero = _mm256_setzero_ps();
  auto sse_1em15 = _mm256_set1_ps(1e-15f);
  auto pmin_zero = _mm256_set1_ps(pmin[0]);
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq

    auto sub1_loop = _mm256_sub_ps(sumofsquares, sigmas_loop);
    auto add2_loop = _mm256_add_ps(sse_1em15, sumofsquares);
    auto rcp1_loop = _mm256_rcp_ps(add2_loop);
    auto mul3_loop = _mm256_mul_ps(rcp1_loop, sub1_loop);
    auto max1_loop = _mm256_max_ps(zero, mul3_loop);
    auto pow1_loop = fmath::pow_ps256(max1_loop, pmin_zero);
    auto mul4_loop = _mm256_mul_ps(dftc_loop, pow1_loop);
    _mm256_storeu_ps(dftc + h, mul4_loop);
  }
  _mm256_zeroupper();
}

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("fma,avx2")))
#endif
void filter_6_AVX2_8(float* dftc, const float* sigmas, const int ccnt,
  const float* pmin, const float* pmax, const float* sigmas2)
{
  auto zero = _mm256_setzero_ps();
  auto sse_1em15 = _mm256_set1_ps(1e-15f);
  for (int h = 0; h < ccnt; h += 8)
  {
    auto dftc_loop = _mm256_loadu_ps(dftc + h); // R I r i
    auto sigmas_loop = _mm256_loadu_ps(sigmas + h); // S S s s adjacent Sigma values should be the same!
    auto squares = _mm256_mul_ps(dftc_loop, dftc_loop); //R*R I*I r*r i*i

    auto sumofsquares = _mm256_add_ps(squares, _mm256_shuffle_ps(squares, squares, _MM_SHUFFLE(2, 3, 0, 1)));
    // R*R+I*I R*R+I*I r*r+i*i r*r+i*i from now SQ SQ sq sq

    auto sub1_loop = _mm256_sub_ps(sumofsquares, sigmas_loop);
    auto add2_loop = _mm256_add_ps(sse_1em15, sumofsquares);
    auto rcp1_loop = _mm256_rcp_ps(add2_loop);
    auto mul3_loop = _mm256_mul_ps(rcp1_loop, sub1_loop);
    //auto max1_loop = _mm256_max_ps(sse_1em15, mul3_loop); // bug fixed in 1.9.6: sse_1em15 instead of zero in max
    auto max1_loop = _mm256_max_ps(zero, mul3_loop);
#if 0
    auto rcp2_loop = _mm256_rcp_ps(_mm256_rsqrt_ps(max1_loop)); // 1.9.6: replacing rsqrt + rcp
#else
    auto rcp2_loop = _mm256_sqrt_ps(max1_loop); // 1.9.6: replacing rsqrt + rcp
#endif
    auto mul4_loop = _mm256_mul_ps(dftc_loop, rcp2_loop);
    _mm256_storeu_ps(dftc + h, mul4_loop);
  }
  _mm256_zeroupper();
}

void intcast_float_to_stacked16_AVX2(const float* p, unsigned char* dst, unsigned char* dst_lsb, const int src_height,
    const int src_width, const int dst_pitch, const int width)
{
    //const int saved_rounding_mode = ::_controlfp(RC_NEAR, _MCW_RC);
    for (int y = 0; y < src_height; ++y)
    {
        int x = 0;
        for (; x < src_width - 15; x += 16)
        {
            //const float vf = p[x] * 256;
            //int v;
            //v = int(vf + 0.5f);
            //v = min(max(v, 0), 65535);
            auto vf1 = _mm256_loadu_ps(p + x);
            auto vf2 = _mm256_loadu_ps(p + x + 8);
            vf1 = _mm256_mul_ps(vf1, _mm256_set1_ps(256.0));
            vf2 = _mm256_mul_ps(vf2, _mm256_set1_ps(256.0));
            auto v1 = _mm256_cvtps_epi32(vf1);    //cvtpsはMXCSRに従うが、デフォルトで一番近い整数なので四捨五入しなくてよい
            auto v2 = _mm256_cvtps_epi32(vf2);
            //vf1 = _mm256_fmadd_ps(vf1, _mm256_set1_ps(256.0),  _mm256_set1_ps(0.5));  //
            //vf2 = _mm256_fmadd_ps(vf2, _mm256_set1_ps(256.0),  _mm256_set1_ps(0.5));  //
            //auto v1 = _mm256_cvttps_epi32(vf1);
            //auto v2 = _mm256_cvttps_epi32(vf2);    //cvttpsは切り捨てなので、元の式からすると0.5足してこっちになるが・・・

            auto v = _mm256_packus_epi32(v1, v2);    //飽和付きで16bitへ なんか31|20のような順番になるので注意
            //ymmレーンが32bit単位に集めないと超えれないのでとりあえず集める
            v = _mm256_shuffle_epi8(v, _mm256_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0, 15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0));    //3H 1H 3L 1L 2H 0H 2L 0L
            v = _mm256_permutevar8x32_epi32(v, _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0));  //3H 2H 1H 0H 3L 2L 1L 0L  ymmレーン越え


            //dst[x] = static_cast <unsigned char> (v >> 8);
            //dst_lsb[x] = static_cast <unsigned char> (v);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + x), _mm256_extracti128_si256(v, 1));   //上位
            _mm_storeu_si128(reinterpret_cast<__m128i*>(dst_lsb + x), _mm256_extracti128_si256(v, 0));//下位

        }


        for (; x < src_width; ++x)
        {
            const float vf = p[x] * 256;
            int v;
            v = int(vf + 0.5f);
            v = min(max(v, 0), 65535);
            dst[x] = static_cast <unsigned char> (v >> 8);
            dst_lsb[x] = static_cast <unsigned char> (v);
        }
        p += width;
        dst += dst_pitch;
        dst_lsb += dst_pitch;
    }

    //::_controlfp(saved_rounding_mode, _MCW_RC);
}
