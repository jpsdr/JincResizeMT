#ifndef __JincResize_AVX2_H__
#define __JincResize_AVX2_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_avx2(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch,
	const float ValMin, const float ValMax);

#endif
