/*
**                    dfttest v1.9.4.3 for Avisynth+
**
**   2D/3D frequency domain denoiser.
**
**   Copyright (C) 2007-2010 Kevin Stone, 2017 (C) DJATOM
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

#ifdef AVX_BUILD
#include "dfttest_avx.h"

void removeMean_AVX_8(float *dftc, const float *dftgc, const int ccnt, float *dftc2)
{
	auto gf_asm = _mm256_set1_ps(dftc[0] / dftgc[0]);

	for (int h=0; h < ccnt; h += 8)
	{
		auto dftgc_loop  = _mm256_loadu_ps(dftgc + h);
		auto dftc_loop   = _mm256_loadu_ps(dftc + h);

		auto dftc2_result= _mm256_mul_ps(gf_asm, dftgc_loop);
		auto dftc_result = _mm256_sub_ps(dftc_loop, dftc2_result);

		_mm256_storeu_ps(dftc2 + h, dftc2_result);
		_mm256_storeu_ps(dftc  + h, dftc_result);
	}
	_mm256_zeroupper();

}

void addMean_AVX_8(float *dftc, const int ccnt, const float *dftc2)
{
	for (int h=0; h < ccnt; h += 8)
	{
		auto dftc_loop = _mm256_loadu_ps(dftc + h);
		auto dftc2_loop = _mm256_loadu_ps(dftc2 + h);
		auto dftc_result = _mm256_add_ps(dftc2_loop, dftc_loop);
		_mm256_storeu_ps(dftc + h, dftc_result);
	}
	_mm256_zeroupper();

}

void proc0_AVX_8(const unsigned char *s0, const float *s1, float *d,
	const int p0, const int p1, const int /*offset_lsb*/)
{
	auto maskl = _mm_set_epi8(-1,-1,-1,3,-1,-1,-1,2,-1,-1,-1,1,-1,-1,-1,0);
	auto maskh = _mm_set_epi8(-1,-1,-1,7,-1,-1,-1,6,-1,-1,-1,5,-1,-1,-1,4);

	for (int u = 0; u < p1; ++u)
	{
		for (int v = 0; v < p1; v += 8)
		{
			auto s1f = _mm256_loadu_ps(s1 + v);
			auto s064= _mm_loadu_si64(s0 + v);	//8å¬ÇÃÇ›
			auto s0il= _mm_shuffle_epi8(s064,maskl);	//char->int â∫à 8å¬(64bit)ÇÃïœä∑ (1)
			auto s0ih= _mm_shuffle_epi8(s064,maskh);
			auto s0i = _mm256_castsi128_si256(s0il);		//â∫à ÇymmÇ…ÉäÉlÅ[ÉÄ(0)
			     s0i = _mm256_insertf128_si256(s0i,s0ih,1);	//è„à Ç…xmmÇë}ì¸(3latency)
			auto s0f = _mm256_cvtepi32_ps(s0i);				//int->float
			auto d_res = _mm256_mul_ps(s0f, s1f);
			_mm256_storeu_ps(d + v, d_res);
		}
		s0 += p0;
		s1 += p1;
		d += p1;
	}
	_mm256_zeroupper();
}

void proc0_AVX2_8(const unsigned char *s0, const float *s1, float *d,
	const int p0, const int p1, const int /*offset_lsb*/)
{
	for (int u = 0; u < p1; ++u)
	{
		for (int v = 0; v < p1; v += 8)
		{
			auto s064 = _mm_loadu_si64(s0 + v);	//8å¬ÇÃÇ›
			auto s1f = _mm256_loadu_ps(s1 + v);

			auto s0i = _mm256_cvtepu8_epi32(s064);	//char->int â∫à 8å¬(64bit)ÇÃïœä∑	AVX2ñΩóﬂ
			auto s0f = _mm256_cvtepi32_ps(s0i);			//int->float
			auto d_res = _mm256_mul_ps(s0f, s1f);
			_mm256_storeu_ps(d + v, d_res);
		}
		s0 += p0;
		s1 += p1;
		d += p1;
	}
	_mm256_zeroupper();
}

void proc1_AVX_8(const float *s0, const float *s1, float *d,
	const int p0, const int p1)
{
#ifdef _WIN64	//Ç»Ç∫Ç©IntrinsicsÇæÇ∆asmÇÊÇËíxÇ¢ loopèàóùÇÃç∑ÅH
	for (int u = 0; u < p0; u++)
	{
		for (int v = 0; v < p0; v += 8)
		{
			auto ymm0 = _mm256_loadu_ps(s0 + v);
			auto ymm1 = _mm256_loadu_ps(s1 + v);
			auto ymm2 = _mm256_loadu_ps(d + v);

			ymm0 = _mm256_mul_ps(ymm0,ymm1);
			ymm0 = _mm256_add_ps(ymm0,ymm2);
			_mm256_storeu_ps(d + v, ymm0);
		}
		s0 += p0;
		s1 += p0;
		d += p1;
	}
#else
	__asm
	{
		mov esi,s0
		mov edi,s1
		mov edx,d
		mov eax,p0

		mov ecx,p0
		xor ebx,ebx
	uloop:
	vloop:
		vmovups ymm0,[esi+ebx]
		vmovups ymm1,[edi+ebx]
		vmovups ymm2,[edx+ebx]

		vmulps ymm0,ymm0,ymm1
		vaddps ymm0,ymm0,ymm2
		vmovups [edx+ebx],ymm0

		add ebx,32
		sub ecx,8
		jg vloop
		mov ecx,p1
		lea edx,[edx+ecx*4]
		mov ecx,p0
		lea esi,[esi+ecx*4]
		lea edi,[edi+ecx*4]
		xor ebx,ebx
		sub eax,1
		jg uloop
	}
#endif
	_mm256_zeroupper();
}

void proc1_AVX2_8(const float *s0, const float *s1, float *d,
	const int p0, const int p1)
{
#ifdef _WIN64
	for (int u = 0; u < p0; u++)
	{
		for (int v = 0; v < p0; v += 8)
		{
			auto ymm0 = _mm256_loadu_ps(s0 + v);
			auto ymm1 = _mm256_loadu_ps(s1 + v);
			auto ymm2 = _mm256_loadu_ps(d + v);

			ymm0 = _mm256_fmadd_ps(ymm0, ymm1, ymm2);	//å¬ï Ç…åvéZÇ∑ÇÈÇÊÇËê∏ìxÇÕçÇÇ≠Ç»ÇÈ(åãâ ÇÕè]óàÇ©ÇÁïœÇÌÇÈ)ÇÃÇ≈íçà” FMA3(ïÇìÆè¨êîÇÕHaswellà»ç~Å‡AVX2 êÆêîÇ»ÇÁSandyà»ç~)Ç™égÇ¶ÇÈÇ»ÇÁ2ñΩóﬂÇ1Ç¬Ç…Ç≈Ç´ÇÈ 
			_mm256_storeu_ps(d + v, ymm0);
		}
		s0 += p0;
		s1 += p0;
		d += p1;
	}
#else
	__asm
	{
		mov esi,s0
		mov edi,s1
		mov edx,d
		mov eax,p0

		mov ecx,p0
		xor ebx,ebx
	uloop:
	vloop:
		vmovups ymm0,[esi+ebx]
		vmovups ymm1,[edi+ebx]
		vmovups ymm2,[edx+ebx]

		//vmulps ymm0,ymm0,ymm1
		//vaddps ymm0,ymm0,ymm2
		vfmadd213ps ymm0,ymm1,ymm2
		vmovups [edx+ebx],ymm0

		add ebx,32
		sub ecx,8
		jg vloop
		mov ecx,p1
		lea edx,[edx+ecx*4]
		mov ecx,p0
		lea esi,[esi+ecx*4]
		lea edi,[edi+ecx*4]
		xor ebx,ebx
		sub eax,1
		jg uloop
	}
#endif
	_mm256_zeroupper();
}


void filter_0_AVX(float *dftc, const float *sigmas, const int ccnt,
	const float *pmin, const float *pmax, const float *sigmas2)
{
	register auto avx_1em15 = _mm256_set1_ps(1e-15f);
#ifdef _WIN64
	register auto zero = _mm256_setzero_ps();
	for (int h = 0; h<ccnt; h += 8)
	{
		auto dft = _mm256_loadu_ps(dftc + h);		//dftc[h+ 3,2,1,0]
		auto psd = _mm256_mul_ps(dft, dft);			//dftc[h+ 3,2,1,0].^2
		auto psd2= _mm256_permute_ps(psd, 177);		//dftc[h+ 2,3,0,1].^2
		     psd = _mm256_add_ps(psd, psd2);		//psd=dftc[h+ 3,2,1,0].^2

		auto den = _mm256_add_ps(psd, avx_1em15);	// psd+1e-15f
		     den = _mm256_rcp_ps(den);				// 1/(psd+1e-15f)
		auto sig = _mm256_loadu_ps(sigmas + h);
		auto num = _mm256_sub_ps(psd, sig);			// psd-sigmas
		auto res = _mm256_mul_ps(num, den);			// (psd-sigmas[h])/(psd+1e-15f)
		     res = _mm256_max_ps(res,zero);			// max(res,0)
		     res = _mm256_mul_ps(res,dft);			// res *= dftc
		_mm256_storeu_ps(dftc + h, res);
	}
#else
	__asm
	{
		mov esi,dftc
		mov edi,sigmas
		mov eax,ccnt
		xor ecx,ecx
		vxorps ymm4, ymm4, ymm4	//zero
		vmovaps ymm5,avx_1em15
	hloop:
		vmovups ymm1,[esi+ecx]	//dftc
		vmulps	ymm2, ymm1, ymm1	//dftc[h+ 3,2,1,0].^2
		vpermilps ymm3, ymm2, 177
		vaddps	ymm2, ymm2, ymm3	//psd

		vaddps	ymm3, ymm2, ymm5	//den=psd+1e-15f
		vrcpps	ymm3, ymm3

		vmovups ymm0,[edi+ecx]	//sigmas
		vsubps	ymm2, ymm2, ymm0	//num=psd-sigmas

		vmulps	ymm0, ymm2, ymm3	//num/den
		vmaxps	ymm0, ymm0, ymm4	//max num,zero

		vmulps	ymm0, ymm0, ymm1	//dftc*=res

		vmovups [esi+ecx],ymm0

		add ecx,32
		sub eax,8
		jg hloop

	}
#endif
	_mm256_zeroupper();
}

void intcast_AVX2(const float *p, unsigned char *dst, const int src_height,
	const int src_width, const int dst_pitch, const int width)
{
	auto avx_05 = _mm256_set1_ps(0.5f);
	for (int y = 0; y<src_height; ++y)
	{
		int x=0;
		for (; x < src_width; x += 16)
		{
			auto p_loop = _mm256_loadu_ps(p + x);
			auto p_loop2 = _mm256_loadu_ps(p + x + 8);
			auto add_loop = _mm256_add_ps(avx_05, p_loop);
			auto add_loop2 = _mm256_add_ps(avx_05, p_loop2);
			auto int1_loop = _mm256_cvttps_epi32(add_loop);
			auto int1_loop2 = _mm256_cvttps_epi32(add_loop2);
			auto packs_loop = _mm256_packs_epi32(int1_loop, int1_loop2);
			auto result = _mm256_packus_epi16(packs_loop, packs_loop);
			result = _mm256_permutevar8x32_epi32(result,_mm256_set_epi32(0,0,0,0,5,1,4,0));
			_mm_storeu_si128(reinterpret_cast<__m128i *>(dst + x), _mm256_castsi256_si128(result));
		}
		for (; x<src_width; ++x){
			//dst[x] = min(max((int)(p[x]+0.5f),0),255);
			int tmp = (int)(p[x]+0.5f);
			dst[x] = (tmp <0) ? 0 : (tmp >255) ? 255 : tmp;
		}
		p += width;
		dst += dst_pitch;
	}
	_mm256_zeroupper();
}


#endif

