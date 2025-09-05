// VS 2013
#if _MSC_VER >= 1800

#include <immintrin.h>
#include "JincRessizeMT.h"

template <typename T>
void resize_plane_avx2(const MT_Data_Info_JincResizeMT *MT_DataGF, const uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[])
{
	const T* src = reinterpret_cast<const T*>(MT_DataGF->src[idxPlane]);
	T* JincMT_RESTRICT dst = reinterpret_cast<T*>(MT_DataGF->dst[idxPlane]);

	const ptrdiff_t src_pitch = (ptrdiff_t)MT_DataGF->src_pitch[idxPlane] / sizeof(T);
	const ptrdiff_t dst_pitch = (ptrdiff_t)MT_DataGF->dst_pitch[idxPlane] / sizeof(T);

	const int Y_Min = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF->dst_UV_h_min : MT_DataGF->dst_Y_h_min;
	const int Y_Max = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF->dst_UV_h_max : MT_DataGF->dst_Y_h_max;
	const int dst_width = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF->dst_UV_w : MT_DataGF->dst_Y_w;

	const T val_min = static_cast<T>(Val_Min[idxPlane]);
	const T val_max = static_cast<T>(Val_Max[idxPlane]);
	const __m256 min_val = _mm256_set1_ps(Val_Min[idxPlane]);

	EWAPixelCoeffMeta *meta_y = coeff->meta + (Y_Min*dst_width);

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

    for (int y = Y_Min; y < Y_Max; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

        for (int x = 0; x < dst_width; ++x)
        {
            const T *src_ptr = src + (meta->start_y * src_pitch + meta->start_x);
            const float *coeff_ptr = coeff->factor + meta->coeff_meta;
            __m256 result = _mm256_setzero_ps();

            if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
                        const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
                        const __m256 coeff = _mm256_load_ps(coeff_ptr + lx);
                        result = _mm256_fmadd_ps(src_ps, coeff, result);
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()), _mm_setzero_si128()));
				dst[x] = clamp(final_res,val_min,val_max);
            }
            else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
                        const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
                        const __m256 coeff = _mm256_load_ps(coeff_ptr + lx);
                        result = _mm256_fmadd_ps(src_ps, coeff, result);
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()));
				dst[x] = clamp(final_res,val_min,val_max);
            }
            else
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
                        const __m256 src_ps = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr + lx)), min_val);
                        const __m256 coeff = _mm256_load_ps(coeff_ptr + lx);
                        result = _mm256_fmadd_ps(src_ps, coeff, result);
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));
				dst[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum)));
            }
			meta++;
        } // for (x)
		meta_y += dst_width;
        dst += dst_pitch;
	} // for (y)
}

template void resize_plane_avx2<uint8_t>(const MT_Data_Info_JincResizeMT *MT_DataGF, const  uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2<uint16_t>(const MT_Data_Info_JincResizeMT *MT_DataGF, const uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2<float>(const MT_Data_Info_JincResizeMT *MT_DataGF, const uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);

#endif