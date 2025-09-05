#ifndef __JincResize_AVX512_H__
#define __JincResize_AVX512_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_avx512(const MT_Data_Info_JincResizeMT *MT_DataGF, const uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);

#endif
