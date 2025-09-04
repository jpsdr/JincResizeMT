#ifndef __JincResize_AVX2_H__
#define __JincResize_AVX2_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_avx2(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);

#endif
