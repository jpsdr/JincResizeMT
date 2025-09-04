#ifndef __JINCRESIZEMT_H__
#define __JINCRESIZEMT_H__

#include <string>
#include <vector>
#include <immintrin.h>

#include "avisynth.h"
#include "avs/minmax.h"
#include "ThreadPoolInterface.h"

#define JINCRESIZEMT_VERSION "JincResizeMT 1.0.0 JPSDR"

#define JincMT_RESTRICT __restrict

// VS 2017 v15.3
#if _MSC_VER >= 1911
#define JincMT_CONSTEXPR constexpr
#else
#define JincMT_CONSTEXPR
#endif


typedef struct _MT_Data_Info_JincResizeMT
{
	const BYTE* src[4];
	BYTE *JincMT_RESTRICT dst[4];
	int src_pitch[4];
	int dst_pitch[4];
	int32_t src_Y_h_min, src_Y_h_max, src_Y_w;
	int32_t src_UV_h_min, src_UV_h_max, src_UV_w;
	int32_t dst_Y_h_min, dst_Y_h_max, dst_Y_w;
	int32_t dst_UV_h_min, dst_UV_h_max, dst_UV_w;
	bool top, bottom;
} MT_Data_Info_JincResizeMT;

struct EWAPixelCoeffMeta
{
    int start_x;
    int start_y;
    int coeff_meta;
};

struct EWAPixelCoeff
{
    float* factor;
    EWAPixelCoeffMeta* meta;
    int* factor_map;
    int filter_size;
    int coeff_stride;
	
	EWAPixelCoeff() : factor(nullptr), meta(nullptr), factor_map(nullptr) {}
};

#define LUT_SIZE_VALUE 1024

class Lut
{
    int lut_size;

public:
    Lut();
    void InitLut(int lut_size, double radius, double blur);
    float GetFactor(int index);

    double* lut;
};


class JincResizeMT : public GenericVideoFilter
{
    Lut *init_lut;
	std::vector<EWAPixelCoeff*> out;
    bool avx512,avx2,sse41;
    int planecount;
    bool has_at_least_v8,has_at_least_v11;
	bool grey,isRGBPfamily,isAlphaChannel;
	uint8_t bits_per_pixel;
	bool subsampled;
    float ValMin[4],ValMax[4];

    void(*process_frame)(MT_Data_Info_JincResizeMT MT_DataGF, uint8_t idxPlane, EWAPixelCoeff *coeff, const float ValMin, const float ValMax);

	Public_MT_Data_Thread MT_Thread[MAX_MT_THREADS];
	MT_Data_Info_JincResizeMT MT_Data[MAX_MT_THREADS];
	uint8_t threads,threads_number;
	bool sleep;
	uint32_t UserId;
	
	ThreadPoolFunction Jinc_MT;

	static void StaticThreadpool(void *ptr);

	void FreeData(void);

public:
	JincResizeMT(PClip _child, int target_width, int target_height, double crop_left, double crop_top, double crop_width, double crop_height,
		int quant_x, int quant_y, int tap, double blur, int opt, int range, uint8_t _threads, bool _sleep, bool negativePrefetch,
		IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    ~JincResizeMT();
	int __stdcall SetCacheHints(int cachehints, int frame_range);
};

#endif
