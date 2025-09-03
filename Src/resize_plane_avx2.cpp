#include <immintrin.h>

// VS 2013
#if _MSC_VER >= 1800

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_avx2(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch,
	const float ValMin, const float ValMax)
{
    const T *src = reinterpret_cast<const T*>(src_);
    T* VS_RESTRICT dst = reinterpret_cast<T*>(dst_);

    src_pitch /= sizeof(T);
    dst_pitch /= sizeof(T);

	const T val_min = static_cast<T>(ValMin);
	const T val_max = static_cast<T>(ValMax);
	const __m256 min_val = _mm256_set1_ps(ValMin);

	EWAPixelCoeffMeta *meta_y = coeff->meta;

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

    for (int y = 0; y < dst_height; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

        for (int x = 0; x < dst_width; ++x)
        {
            //EWAPixelCoeffMeta* meta = coeff->meta + static_cast<int64_t>(y) * dst_width + x;
            const T *src_ptr = src + (meta->start_y * static_cast<int64_t>(src_pitch)) + meta->start_x;
            const float *coeff_ptr = coeff->factor + meta->coeff_meta;
            __m256 result = _mm256_setzero_ps();

            if VS_CONSTEXPR (std::is_same<T, uint8_t>::value)
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
            else if VS_CONSTEXPR (std::is_same<T, uint16_t>::value)
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

template void resize_plane_avx2<uint8_t>(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch, const float ValMin, const float ValMax);
template void resize_plane_avx2<uint16_t>(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch, const float ValMin, const float ValMax);
template void resize_plane_avx2<float>(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch, const float ValMin, const float ValMax);

#endif