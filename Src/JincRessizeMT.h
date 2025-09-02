#ifndef __JINCRESIZEMT_H__
#define __JINCRESIZEMT_H__

#include <string>
#include <vector>
#include <immintrin.h>

#include "avisynth.h"
#include "avs/minmax.h"
#include "ThreadPoolInterface.h"

#define JINCRESIZEMT_VERSION "JincResizeMT 1.0.0 JPSDR"

#define myfree(ptr) if (ptr!=nullptr) { free(ptr); ptr=nullptr;}
#define myalignedfree(ptr) if (ptr!=nullptr) { _aligned_free(ptr); ptr=nullptr;}
#define mydeleteT(ptr) if (ptr!=nullptr) { delete[] ptr; ptr=nullptr;}
#define mydelete(ptr) if (ptr!=nullptr) { delete ptr; ptr=nullptr;}

#define VS_RESTRICT __restrict

// VS 2013
#if _MSC_VER >= 1800
#define AVX2_BUILD_POSSIBLE
#endif

// VS 2017 v15.3
#if _MSC_VER >= 1911
#define AVX512_BUILD_POSSIBLE
#define C17_ENABLE
#define VS_CONSTEXPR constexpr
#else
#define VS_CONSTEXPR
#endif


typedef struct _MT_Data_Info_JincResizeMT
{
	const BYTE*src1, *src2, *src3, *src4;
	BYTE *dst1, *dst2, *dst3, *dst4;
	int src_pitch1, src_pitch2, src_pitch3, src_pitch4;
	int dst_pitch1, dst_pitch2, dst_pitch3, dst_pitch4;
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



/*template <typename T>
void resize_plane_sse41(EWAPixelCoeff* coeff, const void* src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch);
template <typename T>
void resize_plane_avx2(EWAPixelCoeff* coeff, const void* src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch);
template <typename T>
void resize_plane_avx512(EWAPixelCoeff* coeff, const void* src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch);
*/

class JincResize : public GenericVideoFilter
{
    int w, h;
    int _opt;
    Lut *init_lut;
	std::vector<EWAPixelCoeff*> out;
    bool avx512,avx2,sse41;
    int planecount;
    bool has_at_least_v8,has_at_least_v11;
	bool grey,isRGBPfamily,isAlphaChannel;
	uint8_t bits_per_pixel;
	bool subsampled;
    float ValMin[4],ValMax[4];

    void(*process_frame)(EWAPixelCoeff* coeff, const void* src_, void* VS_RESTRICT dst_, int dst_width, int dst_height, int src_pitch, int dst_pitch,
		const float ValMin, const float ValMax);

	Public_MT_Data_Thread MT_Thread[MAX_MT_THREADS];
	MT_Data_Info_JincResizeMT MT_Data[MAX_MT_THREADS];
	uint8_t threads,threads_number;
	uint32_t UserId;
	
	ThreadPoolFunction Jinc_MT;

	static void StaticThreadpool(void *ptr);

	uint8_t CreateMTData(uint8_t max_threads, int32_t src_size_x, int32_t src_size_y, int32_t dst_size_x, int32_t dst_size_y, int UV_w, int UV_h);

	void FreeData(void);

public:
    JincResize(PClip _child, int target_width, int target_height, double crop_left, double crop_top, double crop_width, double crop_height,
		int quant_x, int quant_y, int tap, double blur, int opt, int range, uint8_t _threads, bool negativePrefetch, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    ~JincResize();
	int __stdcall SetCacheHints(int cachehints, int frame_range);
};

/*
class Arguments
{
    AVSValue _args[12];
    const char* _arg_names[12];
    int _idx;

public:
    Arguments() : _args{}, _arg_names{}, _idx{} {}

    void add(AVSValue arg, const char* arg_name = nullptr)
    {
        _args[_idx] = arg;
        _arg_names[_idx] = arg_name;
        ++_idx;
    }

    AVSValue args() const { return{ _args, _idx }; }

    const char* const* arg_names() const { return _arg_names; }
};
*/
#endif
