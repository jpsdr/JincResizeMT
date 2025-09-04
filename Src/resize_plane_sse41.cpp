#include <smmintrin.h>
#include "JincRessizeMT.h"

template <typename T>
void resize_plane_sse41(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax)
{
	const T* src = reinterpret_cast<const T*>(MT_DataGF.src[idxPlane]);
	T* JincMT_RESTRICT dst = reinterpret_cast<T*>(MT_DataGF.dst[idxPlane]);

	const int src_pitch = MT_DataGF.src_pitch[idxPlane] / sizeof(T);
	const int dst_pitch = MT_DataGF.dst_pitch[idxPlane] / sizeof(T);

	const int Y_Min = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF.dst_UV_h_min : MT_DataGF.dst_Y_h_min;
	const int Y_Max = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF.dst_UV_h_max : MT_DataGF.dst_Y_h_max;
	const int dst_width = ((idxPlane == 1) || (idxPlane == 2)) ? MT_DataGF.dst_UV_w : MT_DataGF.dst_Y_w;

	const T val_min = static_cast<T>(ValMin);
	const T val_max = static_cast<T>(ValMax);
	const __m128 min_val = _mm_set_ps1(ValMin);
	
	EWAPixelCoeffMeta *meta_y = coeff->meta + Y_Min*dst_width;

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

    for (int y = Y_Min; y < Y_Max; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;
		
        for (int x = 0; x < dst_width; x++)
        {
			const T *src_ptr = src + (meta->start_y * static_cast<int64_t>(src_pitch)) + meta->start_x;
            const float* coeff_ptr = coeff->factor + meta->coeff_meta; 
			__m128 result = _mm_setzero_ps();			

            if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 4)
                    {
                        const __m128 src_ps = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(reinterpret_cast<const int32_t*>(src_ptr + lx)))));
                        const __m128 coeff = _mm_load_ps(coeff_ptr + lx);
                        result = _mm_add_ps(result, _mm_mul_ps(src_ps, coeff));
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                const __m128 hsum = _mm_hadd_ps(_mm_hadd_ps(result, result), _mm_hadd_ps(result, result));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()), _mm_setzero_si128()));
                dst[x] = clamp(final_res,val_min,val_max);
            }
            else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 4)
                    {
                        const __m128 src_ps = _mm_cvtepi32_ps(_mm_cvtepu16_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr + lx))));
                        const __m128 coeff = _mm_load_ps(coeff_ptr + lx);
                        result = _mm_add_ps(result, _mm_mul_ps(src_ps, coeff));
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                const __m128 hsum = _mm_hadd_ps(_mm_hadd_ps(result, result), _mm_hadd_ps(result, result));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()));
                dst[x] = clamp(final_res,val_min,val_max);
            }
            else
            {
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; lx += 4)
                    {
                        const __m128 src_ps = _mm_max_ps(_mm_loadu_ps(reinterpret_cast<const float*>(src_ptr + lx)), min_val);
                        const __m128 coeff = _mm_load_ps(coeff_ptr + lx);
                        result = _mm_add_ps(result, _mm_mul_ps(src_ps, coeff));
                    }

                    coeff_ptr += coeff_stride;
                    src_ptr += src_pitch;
                }

                dst[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(result, result), _mm_hadd_ps(result, result)));
            }
			
			meta++;
        } // for (x)
		meta_y += dst_width;
        dst += dst_pitch;
    } // for (y)
}

template void resize_plane_sse41<uint8_t>(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);
template void resize_plane_sse41<uint16_t>(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);
template void resize_plane_sse41<float>(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);
