#ifndef __JincResize_SSE41_H__
#define __JincResize_SSE41_H__

#include "JincRessizeMT.h"

template <typename T>
void resize_plane_sse41(const MT_Data_Info_JincResizeMT *MT_DataGF, const uint8_t idxPlane, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);

#endif
