#ifndef __JincResize_SSE41_H__
#define __JincResize_SSE41_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_sse41(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);

#endif
