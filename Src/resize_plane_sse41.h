#ifndef __JincResize_SSE41_H__
#define __JincResize_SSE41_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_sse41(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch,
	const float ValMin, const float ValMax);

#endif
