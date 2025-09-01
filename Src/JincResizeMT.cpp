#define __STDCPP_WANT_MATH_SPEC_FUNCS__

#include <cmath>

#include "JincRessizeMT.h"
#include "resize_plane_sse41.h"

// VS 2015
#if _MSC_VER >= 1900
#define AVX2_BUILD_POSSIBLE
#define C17_ENABLE
#endif

// VS 2019 v16.3
#if _MSC_VER >= 1923
#define AVX512_BUILD_POSSIBLE
#endif

static ThreadPoolInterface *poolInterface;

static uint8_t CreateMTData(MT_Data_Info_JincResizeMT MT_Data[],int output,uint8_t threads_number,uint8_t max_threads,int32_t size_x,int32_t size_y)
{
	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=size_y;
		MT_Data[0].dst_Y_h_max=size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		MT_Data[0].src_UV_h_max=size_y >> 1;
		MT_Data[0].dst_UV_h_max=size_y;
		MT_Data[0].src_Y_w=size_x;
		if (output==0) MT_Data[0].dst_Y_w=size_x << 1;
		else MT_Data[0].dst_Y_w=size_x;
		MT_Data[0].src_UV_w=size_x >> 1;
		MT_Data[0].dst_UV_w=size_x >> 1;
		return(1);
	}

	int32_t dh_Y,dh_UV,h_y;
	uint8_t i,max=1;

	dh_Y=(size_y+(int32_t)max_threads-1)/(int32_t)max_threads;
	if (dh_Y<16) dh_Y=16;
	if ((dh_Y & 3)!=0) dh_Y=((dh_Y+3) >> 2) << 2;

	h_y=dh_Y;
	while (h_y<(size_y-16))
	{
		max++;
		h_y+=dh_Y;
	}

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=size_y;
		MT_Data[0].dst_Y_h_max=size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		MT_Data[0].src_UV_h_max=size_y >> 1;
		MT_Data[0].dst_UV_h_max=size_y;
		MT_Data[0].src_Y_w=size_x;
		if (output==0) MT_Data[0].dst_Y_w=size_x << 1;
		else MT_Data[0].dst_Y_w=size_x;
		MT_Data[0].src_UV_w=size_x >> 1;
		MT_Data[0].dst_UV_w=size_x >> 1;
		return(1);
	}

	dh_UV=dh_Y>>1; 

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dh_Y;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dh_Y;
		i++;
	}
	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=size_y;
	MT_Data[max-1].dst_Y_h_max=size_y;
	MT_Data[max-1].src_UV_h_max=size_y >> 1;
	MT_Data[max-1].dst_UV_h_max=size_y;
	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=size_x;
		if (output==0) MT_Data[i].dst_Y_w=size_x << 1;
		else MT_Data[i].dst_Y_w=size_x;
		MT_Data[i].src_UV_w=size_x >> 1;
		MT_Data[i].dst_UV_w=size_x >> 1;
	}
	return(max);
}

static AVS_FORCEINLINE unsigned portable_clz(size_t x)
{
    unsigned long index;
    return (_BitScanReverse(&index, static_cast<unsigned long>(x))) ? (31 - index) : 32;
}


#ifndef M_PI // GCC seems to have it
static double M_PI = 3.14159265358979323846;
#endif

// Taylor series coefficients of 2*BesselJ1(pi*x)/(pi*x) as (x^2) -> 0
static double jinc_taylor_series[31] =
{
    1.0,
    -1.23370055013616982735431137,
    0.507339015802096027273126733,
    -0.104317403816764804365258186,
    0.0128696438477519721233840271,
    -0.00105848577966854543020422691,
    6.21835470803998638484476598e-05,
    -2.73985272294670461142756204e-06,
    9.38932725442064547796003405e-08,
    -2.57413737759717407304931036e-09,
    5.77402672521402031756429343e-11,
    -1.07930605263598241754572977e-12,
    1.70710316782347356046974552e-14,
    -2.31434518382749184406648762e-16,
    2.71924659665997312120515390e-18,
    -2.79561335187943028518083529e-20,
    2.53599244866299622352138464e-22,
    -2.04487273140961494085786452e-24,
    1.47529860450204338866792475e-26,
    -9.57935105257523453155043307e-29,
    5.62764317309979254140393917e-31,
    -3.00555258814860366342363867e-33,
    1.46559362903641161989338221e-35,
    -6.55110024064596600335624426e-38,
    2.69403199029404093412381643e-40,
    -1.02265499954159964097119923e-42,
    3.59444454568084324694180635e-45,
    -1.17313973900539982313119019e-47,
    3.56478606255557746426034301e-50,
    -1.01100655781438313239513538e-52,
    2.68232117541264485328658605e-55
};

static double jinc_zeros[16] =
{
    1.2196698912665045,
    2.2331305943815286,
    3.2383154841662362,
    4.2410628637960699,
    5.2427643768701817,
    6.2439216898644877,
    7.2447598687199570,
    8.2453949139520427,
    9.2458926849494673,
    10.246293348754916,
    11.246622794877883,
    12.246898461138105,
    13.247132522181061,
    14.247333735806849,
    15.247508563037300,
    16.247661874700962
};

//  Modified from boost package math/tools/`rational.hpp`
//
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double evaluate_rational(const double* num, const double* denom, double z, int count)
{
    double s1, s2;
    if (z <= 1.0)
    {
        s1 = num[count - 1];
        s2 = denom[count - 1];
        for (auto i = count - 2; i >= 0; --i)
        {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    }
    else
    {
        z = 1.0f / z;
        s1 = num[0];
        s2 = denom[0];
        for (auto i = 1; i < count; ++i)
        {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    }

    return s1 / s2;
}

//  Modified from boost package `BesselJ1.hpp`
//
//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double jinc_sqr_boost_l(double x2)
{
    double bPC[7] =
    {
        -4.4357578167941278571e+06,
        -9.9422465050776411957e+06,
        -6.6033732483649391093e+06,
        -1.5235293511811373833e+06,
        -1.0982405543459346727e+05,
        -1.6116166443246101165e+03,
        0.0
    };
    double bQC[7] =
    {
        -4.4357578167941278568e+06,
        -9.9341243899345856590e+06,
        -6.5853394797230870728e+06,
        -1.5118095066341608816e+06,
        -1.0726385991103820119e+05,
        -1.4550094401904961825e+03,
        1.0
    };
    double bPS[7] =
    {
        3.3220913409857223519e+04,
        8.5145160675335701966e+04,
        6.6178836581270835179e+04,
        1.8494262873223866797e+04,
        1.7063754290207680021e+03,
        3.5265133846636032186e+01,
        0.0
    };
    double bQS[7] =
    {
        7.0871281941028743574e+05,
        1.8194580422439972989e+06,
        1.4194606696037208929e+06,
        4.0029443582266975117e+05,
        3.7890229745772202641e+04,
        8.6383677696049909675e+02,
        1.0
    };

    const auto y2 = M_PI * M_PI * x2;
    const auto xp = sqrt(y2);
    const auto y2p = 64.0 / y2;
    const auto sx = sin(xp);
    const auto cx = cos(xp);

    return (sqrt(xp / M_PI) * 2.0 / y2) * (evaluate_rational(bPC, bQC, y2p, 7) * (sx - cx) + (8.0 / xp) * evaluate_rational(bPS, bQS, y2p, 7) * (sx + cx));
}


#define MAX_TERMS 50
static double bessel_j1(double x)
{
    const double EPS = 1e-15;
	double TabD[MAX_TERMS];

    double term = x / 2.0; // premier terme m=0

	TabD[0] = term;

    double x2 = (x * x) / 4.0;

    unsigned long m=1;
	while ((std::fabs(term) >= EPS) && (m<MAX_TERMS))
	{
		term *= -x2 / (m * (m + 1)); // récurrence
		TabD[m++] = term;
	}

	double sum = 0.0;

	for (int i=m-1; i>=0; i--)
		sum += TabD[i];
	
    return sum;
}


// jinc(sqrt(x2))
static double jinc_sqr(double x2)
{
    if (x2 < 1.49)        // the 1-tap radius
    {
        double res = 0.0;
        for (auto j = 16; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 4.97)   // the 2-tap radius
    {
        double res = 0.0;
        for (auto j = 21; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 10.49)  // the 3-tap radius
    {
        double res = 0.0;
        for (auto j = 26; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 17.99)  // the 4-tap radius
    {
        double res = 0.0;
        for (auto j = 31; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 52.57)  // the 5~7-tap radius
    {
        const auto x = M_PI * sqrt(x2);
#ifdef C17_ENABLE
        return 2.0 * std::cyl_bessel_j(1, x) / x;
#else
		return 2.0 * bessel_j1(x)/x;
#endif
    }
    else if (x2 < 68.07)  // the 8-tap radius // Modify from pull request #4
    {
        return jinc_sqr_boost_l(x2);
    }
    else                  // the 9~16-tap radius
    {
        const auto x = M_PI * sqrt(x2);
#ifdef C17_ENABLE
        return 2.0 * std::cyl_bessel_j(1, x) / x;
#else
		return 2.0 * bessel_j1(x)/x;
#endif
    }
}


static double sample_sqr(double (*filter)(double), double x2, double blur2, double radius2)
{
    if (blur2 > 0.0)
        x2 /= blur2;

    if (x2 < radius2)
        return filter(x2);

    return 0.0;
}

double JINC_ZERO_SQR = 1.48759464366204680005356;

Lut::Lut() : lut_size(LUT_SIZE_VALUE)
{
    lut = new double[lut_size];
}

void Lut::InitLut(int lut_size, double radius, double blur)
{
    const auto radius2 = radius * radius;
    const auto blur2 = blur * blur;

    for (auto i = 0; i < lut_size; ++i)
    {
        const auto t2 = i / (lut_size - 1.0);
        lut[i] = sample_sqr(jinc_sqr, radius2 * t2, blur2, radius2) * sample_sqr(jinc_sqr, JINC_ZERO_SQR * t2, 1.0, radius2);
    }
}

float Lut::GetFactor(int index)
{
    if (index >= lut_size)
        return 0.f;
    return static_cast<float>(lut[index]);
}

double DOUBLE_ROUND_MAGIC_NUMBER = 6755399441055744.0;

static bool init_coeff_table(EWAPixelCoeff* out, int quantize_x, int quantize_y,
    int filter_size, int dst_width, int dst_height)
{
    out->filter_size = filter_size;
    out->coeff_stride = (filter_size + 15) & ~15;

    // This will be reserved to exact size in coff generating procedure
    out->factor = nullptr;

    // Allocate metadata
    out->meta = new EWAPixelCoeffMeta[static_cast<int64_t>(dst_width) * dst_height];
	if (out->meta == nullptr) return(false);

    // Alocate factor map
    out->factor_map = new int[static_cast<int64_t>(quantize_x) * quantize_y];

    // Zeroed memory
    if (out->factor_map != nullptr)
        memset(out->factor_map, 0, static_cast<int64_t>(quantize_x) * quantize_y * sizeof(int));
	else return(false);

    memset(out->meta, 0, static_cast<int64_t>(dst_width) * dst_height * sizeof(EWAPixelCoeffMeta));
	
	return(true);
}

static void delete_coeff_table(EWAPixelCoeff* out)
{
	myalignedfree(out->factor);
	mydeleteT(out->factor_map);
    mydeleteT(out->meta);
}

struct generate_coeff_params
{
    Lut *func;
    EWAPixelCoeff *out;
    int quantize_x;
    int quantize_y;
    int samples;
    int src_width;
    int src_height;
    int dst_width;
    int dst_height;
    double radius;
    double crop_left;
    double crop_top;
    double crop_width;
    double crop_height;
    int initial_capacity;
    double initial_factor;
};

#ifndef C17_ENABLE
#define round(x) floor(x + 0.5)
#endif

/* Coefficient table generation */
static bool generate_coeff_table_c(const generate_coeff_params &params)
{
    Lut *func = params.func;
    EWAPixelCoeff *out = params.out;
    int quantize_x = params.quantize_x;
    int quantize_y = params.quantize_y;
    int samples = params.samples;
    int src_width = params.src_width;
    int src_height = params.src_height;
    int dst_width = params.dst_width;
    int dst_height = params.dst_height;
    double radius = params.radius;

    const double filter_step_x = min(static_cast<double>(dst_width) / params.crop_width, 1.0);
    const double filter_step_y = min(static_cast<double>(dst_height) / params.crop_height, 1.0);

    const float filter_support_x = static_cast<float>(radius / filter_step_x);
    const float filter_support_y = static_cast<float>(radius / filter_step_y);

    const float filter_support = max(filter_support_x, filter_support_y);
    const int filter_size = max(static_cast<int>(ceil(filter_support_x * 2.0)), static_cast<int>(ceil(filter_support_y * 2.0)));

    const float start_x = static_cast<float>(params.crop_left + (params.crop_width / dst_width - 1.0) / 2.0);

    const float x_step = static_cast<float>(params.crop_width / dst_width);
    const float y_step = static_cast<float>(params.crop_height / dst_height);

    float xpos = start_x;
    float ypos = static_cast<float>(params.crop_top + (params.crop_height - dst_height) / (dst_height * static_cast<int64_t>(2)));

    // Initialize EWAPixelCoeff data structure
    if (!init_coeff_table(out, quantize_x, quantize_y, filter_size, dst_width, dst_height)) return(false);

    size_t tmp_array_capacity = params.initial_capacity;
    float* tmp_array = static_cast<float*>(_aligned_malloc(tmp_array_capacity * sizeof(float), 64));
    if (tmp_array==nullptr) return(false);
    size_t tmp_array_size = 0;
    int tmp_array_top = 0;
    unsigned base_clz = portable_clz(tmp_array_capacity);
    const double initial_growth_factor = params.initial_factor;
    const double radius2 = radius * radius;

    // Use to advance the coeff pointer
    const int coeff_per_pixel = out->coeff_stride * filter_size;

    for (int y = 0; y < dst_height; ++y)
    {
        for (int x = 0; x < dst_width; ++x)
        {
            bool is_border = false;

            EWAPixelCoeffMeta* meta = &out->meta[y * dst_width + x];

            // Here, the window_*** variable specified a begin/size/end
            // of EWA window to process.
            int window_end_x = static_cast<int>(xpos + filter_support);
            int window_end_y = static_cast<int>(ypos + filter_support);

            if (window_end_x >= src_width)
            {
                window_end_x = src_width - 1;
                is_border = true;
            }
            if (window_end_y >= src_height)
            {
                window_end_y = src_height - 1;
                is_border = true;
            }

            int window_begin_x = window_end_x - filter_size + 1;
            int window_begin_y = window_end_y - filter_size + 1;

            if (window_begin_x < 0)
            {
                window_begin_x = 0;
                is_border = true;
            }
            if (window_begin_y < 0)
            {
                window_begin_y = 0;
                is_border = true;
            }

            meta->start_x = window_begin_x;
            meta->start_y = window_begin_y;

            // Quantize xpos and ypos
            const int quantized_x_int = static_cast<int>(xpos * quantize_x);
            const int quantized_y_int = static_cast<int>(ypos * quantize_y);
            const int quantized_x_value = quantized_x_int % quantize_x;
            const int quantized_y_value = quantized_y_int % quantize_y;
            const float quantized_xpos = static_cast<float>(quantized_x_int) / quantize_x;
            const float quantized_ypos = static_cast<float>(quantized_y_int) / quantize_y;

            if (!is_border && out->factor_map[quantized_y_value * quantize_x + quantized_x_value] != 0)
            {
                // Not border pixel and already have coefficient calculated at this quantized position
                meta->coeff_meta = out->factor_map[quantized_y_value * quantize_x + quantized_x_value] - 1;
            }
            else
            {
                // then need computation
                float divider = 0.f;

                // This is the location of current target pixel in source pixel
                // Quantized
                //const float current_x = clamp(is_border ? xpos : quantized_xpos, 0.f, src_width - 1.f);
                //const float current_y = clamp(is_border ? ypos : quantized_ypos, 0.f, src_height - 1.f);

                if (!is_border)
                {
                    // Change window position to quantized position
                    window_begin_x = static_cast<int>(quantized_xpos + filter_support) - filter_size + 1;
                    window_begin_y = static_cast<int>(quantized_ypos + filter_support) - filter_size + 1;
                }

                // Windowing positon
                int window_x = window_begin_x;
                int window_y = window_begin_y;

                // First loop calcuate coeff
                const size_t new_size = tmp_array_size + coeff_per_pixel;
                if (new_size > tmp_array_capacity)
                {
                    size_t new_capacity = tmp_array_capacity * (1.0 + (initial_growth_factor - 1.0)
                        * (1.0 - static_cast<double>(max(0, static_cast<int>(base_clz - portable_clz(tmp_array_capacity)))) / 32.0));
                    if (new_capacity < new_size)
                        new_capacity = new_size;
                    float* new_tmp = static_cast<float*>(_aligned_malloc(new_capacity * sizeof(float), 64));
                    if (new_tmp==nullptr)
                    {
                        myalignedfree(tmp_array);
						return(false);
                    }
                    memcpy(new_tmp, tmp_array, tmp_array_size * sizeof(float));
                    myalignedfree(tmp_array);
                    tmp_array = new_tmp;
                    tmp_array_capacity = new_capacity;
                }
                memset(tmp_array + tmp_array_size, 0, coeff_per_pixel * sizeof(float));
                int curr_factor_ptr = tmp_array_top;
                tmp_array_size = new_size;

                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; ++lx)
                    {
                        // Euclidean distance to sampling pixel
                        const double dx = (clamp(is_border ? xpos : quantized_xpos, 0.f, static_cast<float>(src_width - 1)) - window_x) * filter_step_x;
                        const double dy = (clamp(is_border ? ypos : quantized_ypos, 0.f, static_cast<float>(src_height - 1)) - window_y) * filter_step_y;

                        //int index = static_cast<int>(llround((samples - 1) * (dx * dx + dy * dy) / radius2 + DOUBLE_ROUND_MAGIC_NUMBER));
						int index = static_cast<int>(round((samples - 1) * (dx * dx + dy * dy) / radius2));

                        const float factor = func->GetFactor(index);

                        tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] = factor;
                        divider += factor;

                        ++window_x;
                    }

                    curr_factor_ptr += out->coeff_stride;

                    window_x = window_begin_x;
                    ++window_y;
                }

                // Second loop to divide the coeff
                curr_factor_ptr = tmp_array_top;
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; ++lx)
                    {
                        tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] /= divider;
                    }

                    curr_factor_ptr += out->coeff_stride;
                }

                // Save factor to table
                if (!is_border)
                    out->factor_map[quantized_y_value * quantize_x + quantized_x_value] = tmp_array_top + 1;

                meta->coeff_meta = tmp_array_top;
                tmp_array_top += coeff_per_pixel;
            }

            xpos += x_step;
        }

        ypos += y_step;
        xpos = start_x;
    }

    // Copy from tmp_array to real array
    out->factor = tmp_array;
	
	return(true);
}


/* Planar resampling with coeff table */
/* 8-16 bit */
//#pragma intel optimization_parameter target_arch=sse
template<typename T>
static void resize_plane_c(EWAPixelCoeff *coeff, const void *src_, void* VS_RESTRICT dst_,
    int dst_width, int dst_height, int src_stride, int dst_stride, const float ValMin, const float ValMax)
{
    const T* srcp = reinterpret_cast<const T*>(src_);
    T* VS_RESTRICT dstp = reinterpret_cast<T*>(dst_);

    src_stride /= sizeof(T);
    dst_stride /= sizeof(T);
	
	EWAPixelCoeffMeta *meta_y = coeff->meta;

    for (int y = 0; y < dst_height; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

        for (int x = 0; x < dst_width; x++)
        {
            //EWAPixelCoeffMeta* meta = coeff->meta + static_cast<int64_t>(y) * dst_width + x;
			const T *src_ptr = srcp + meta->start_y * static_cast<int64_t>(src_stride) + meta->start_x;
            const float *coeff_ptr = coeff->factor + meta->coeff_meta;

            float result = 0.0f;

            for (int ly = 0; ly < coeff->filter_size; ly++)
            {
                for (int lx = 0; lx < coeff->filter_size; lx++)
                {
                    result += src_ptr[lx] * coeff_ptr[lx];
                }
                coeff_ptr += coeff->coeff_stride;
                src_ptr += src_stride;
            }

            if (!(std::is_same<T, float>::value))
                dstp[x] = static_cast<T>(round(clamp(result, ValMin, ValMax)));
            else
                dstp[x] = result;

            meta++;
        }
		
		meta_y += dst_width;
        dstp += dst_stride;
    }
}

void JincResize::FreeData(void)
{
	for (int i=0; i<static_cast<int>(out.size()); ++i)
	{
		if (out[i] != nullptr)
		{
			delete_coeff_table(out[i]);
			mydelete(out[i]);
		}
	}
	
	if (init_lut!=nullptr)
	{
		mydeleteT(init_lut->lut);
		mydelete(init_lut);
	}
}

JincResize::JincResize(PClip _child, int target_width, int target_height, double crop_left, double crop_top, double crop_width, double crop_height, int quant_x, int quant_y, int tap, double blur, int opt, IScriptEnvironment* env)
    : GenericVideoFilter(_child), w(target_width), h(target_height), _opt(opt), init_lut(nullptr), has_at_least_v8(false), has_at_least_v11(false), avx512(false), avx2(false), sse41(false), subsampled(false)
{
    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; };

    has_at_least_v11 = true;
    try { env->CheckVersion(11); }
    catch (const AvisynthError&) { has_at_least_v11 = false; };

    if (!vi.IsPlanar())
        env->ThrowError("JincResizeMT: clip must be in planar format.");

    if (tap < 1 || tap > 16)
        env->ThrowError("JincResizeMT: tap must be between 1..16.");

    if (quant_x < 1 || quant_x > 256)
        env->ThrowError("JincResizeMT: quant_x must be between 1..256.");

    if (quant_y < 1 || quant_y > 256)
        env->ThrowError("JincResizeMT: quant_y must be between 1..256.");

    if ((_opt > 3) || (_opt < -1))
        env->ThrowError("JincResizeMT: opt must be between -1..3.");

    if (blur < 0.0 || blur > 10.0)
        env->ThrowError("JincResizeMT: blur must be between 0.0..10.0.");

	if ((!(env->GetCPUFlags() & CPUF_AVX512F) || !has_at_least_v8) && (_opt == 3))
		env->ThrowError("JincResizeMT: opt=3 requires AVX-512F and AVS+.");

	if ((!(env->GetCPUFlags() & CPUF_AVX2) || !has_at_least_v8) && (_opt == 2))
        env->ThrowError("JincResizeMT: opt=2 requires AVX2 and AVS+.");

    if (!(env->GetCPUFlags() & CPUF_SSE4_1) && (_opt == 1))
        env->ThrowError("JincResizeMT: opt=1 requires SSE4.1.");

    int src_width = vi.width;
    int src_height = vi.height;

	if ( vi.Is420() && ( ((target_width%2)!=0) || ((target_height%2)!=0) ) )
	{
		FreeData();
		env->ThrowError("JincResizeMT: width and height must be multiple of 2 for 4:2:0 chroma subsampling.");
	}
	if ( vi.Is422() && ((target_width%2)!=0) )
	{
		FreeData();
		env->ThrowError("JincResizeMT: width must be multiple of 2 for 4:2:2 chroma subsampling.");
	}
	if (vi.IsYV411() && ((target_width%4)!=0) )
	{
		FreeData();
		env->ThrowError("JincResizeMT: width must be multiple of 4 for 4:1:1 chroma subsampling.");
	}

    if (crop_width <= 0.0)
        crop_width = vi.width - crop_left + crop_width;

    if (crop_height <= 0.0)
        crop_height = vi.height - crop_top + crop_height;

    vi.width = w;
    vi.height = h;
	
	const int initial_capacity = max(target_width * target_height, src_width * src_height);
	
	 std::string cplace = "mpeg2";
	
	const double radius = jinc_zeros[tap - 1];
	const int samples = 1024;  // should be a multiple of 4
	
	init_lut = new Lut();
	if (init_lut == nullptr)
	{
		FreeData();
		env->ThrowError("JincResizeMT: Error creating lut.");
	}
	if (init_lut->lut == nullptr)
	{
		FreeData();
		env->ThrowError("JincResizeMT: Error allocating lut.");
	}
	init_lut->InitLut(samples, radius, blur);
	planecount = vi.NumComponents();
	
	const double initial_factor = 1.50;
	 
    out.emplace_back(new EWAPixelCoeff());
    generate_coeff_params params =
    {
        init_lut,
        out[0],
        quant_x,
        quant_y,
        samples,
        src_width,
        src_height,
        target_width,
        target_height,
        radius,
        crop_left,
        crop_top,
        crop_width,
        crop_height,
        initial_capacity,
        initial_factor
    };

	if (!generate_coeff_table_c(params))
	{
		FreeData();
		env->ThrowError("JincResizeMT: Error generating coeff table [0].");
	}
	
	if ((planecount > 1) && !(vi.Is444() || vi.IsRGB()))
	{

        out.emplace_back(new EWAPixelCoeff());
        subsampled = true;
        const int sub_w = vi.GetPlaneWidthSubsampling(PLANAR_U);
        const int sub_h = vi.GetPlaneHeightSubsampling(PLANAR_U);
        const double div_w = 1 << sub_w;
        const double div_h = 1 << sub_h;

        const double crop_left_uv = (cplace == "mpeg2" || cplace == "topleft") ?
            (0.5 * (1.0 - static_cast<double>(src_width) / target_width) + crop_left) / div_w : crop_left / div_w;
        const double crop_top_uv = (cplace == "topleft") ?
            (0.5 * (1.0 - static_cast<double>(src_height) / target_height) + crop_top) / div_h : crop_top / div_h;

        generate_coeff_params params1 = {
            init_lut,
            out[1],
            quant_x,
            quant_y,
            samples,
            src_width >> sub_w,
            src_height >> sub_h,
            target_width >> sub_w,
            target_height >> sub_h,
            radius,
            crop_left_uv,
            crop_top_uv,
            crop_width / div_w,
            crop_height / div_h,
            initial_capacity / (static_cast<int>(div_w) * static_cast<int>(div_h)),
            initial_factor
        };
        if (!generate_coeff_table_c(params1))
		{
			FreeData();
			env->ThrowError("JincResizeMT: Error generating coeff table [1].");
		}
	}

	avx512 = (!!(env->GetCPUFlags() & CPUF_AVX512F) && (_opt < 0) && has_at_least_v8) || (_opt == 3);
    avx2 = (!!(env->GetCPUFlags() & CPUF_AVX2) && (_opt < 0) && has_at_least_v8) || (_opt == 2);
    sse41 = (!!(env->GetCPUFlags() & CPUF_SSE4_1) && (_opt < 0)) || (_opt == 1);
	
	for (unsigned char i=0; i<4; i++)
	{
		 ValMin[i] = 0.0f;
		 ValMax[i] = static_cast<float>((1 << vi.BitsPerPixel()) - 1);
	}

	threads_number = 1;
	threads = 1;


    if (vi.ComponentSize() == 1)
    {
/*        if (avx512)
            process_frame = resize_plane_avx512<uint8_t>;
        else if (avx2)
            process_frame = resize_plane_avx2<uint8_t>;
        else*/
            if (sse41)
				process_frame = resize_plane_sse41<uint8_t>;
			else
				process_frame = resize_plane_c<uint8_t>;
    }
    else if (vi.ComponentSize() == 2)
    {
/*        if (avx512)
            process_frame = resize_plane_avx512<uint16_t>;
        else if (avx2)
            process_frame = resize_plane_avx2<uint16_t>;
        else*/
            if (sse41)
				process_frame = resize_plane_sse41<uint16_t>;
			else
				process_frame = resize_plane_c<uint16_t>;
     }
    else
    {
/*        if (avx512)
            process_frame = resize_plane_avx512<float>;
        else if (avx2)
            process_frame = resize_plane_avx2<float>;
        else*/
            if (sse41)
				process_frame = resize_plane_sse41<float>;
			else
				process_frame = resize_plane_c<float>;
    }
}

JincResize::~JincResize()
{
	if (threads_number>1) poolInterface->RemoveUserId(UserId);
	
	FreeData();

	if (threads>1) poolInterface->DeAllocateAllThreads(true);
}

PVideoFrame JincResize::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src = child->GetFrame(n, env);
    PVideoFrame dst;
    if (has_at_least_v8) dst = env->NewVideoFrameP(vi, &src); else dst = env->NewVideoFrame(vi);

    int planes_y[4] = { PLANAR_Y, PLANAR_U, PLANAR_V, PLANAR_A };
    int planes_r[4] = { PLANAR_G, PLANAR_B, PLANAR_R, PLANAR_A };
    const int* current_planes = (vi.IsYUV() || vi.IsYUVA()) ? planes_y : planes_r;
    for (int i = 0; i < planecount; i++)
    {
        const int plane = current_planes[i];

        int src_stride = src->GetPitch(plane);
        int dst_stride = dst->GetPitch(plane);
        int dst_width = dst->GetRowSize(plane) / vi.ComponentSize();
        int dst_height = dst->GetHeight(plane);
        const uint8_t* srcp = src->GetReadPtr(plane);
        uint8_t* VS_RESTRICT dstp = dst->GetWritePtr(plane);
		
		int i_coeff;
		
		if (subsampled)
		{
			if ((planecount==1) || (planecount==2)) i_coeff = 1;
			else i_coeff = 0;
		}
		else i_coeff = 0;

		process_frame(out[i_coeff],srcp,dstp,dst_width,dst_height,src_stride,dst_stride,ValMin[i],ValMax[i]);
		
/*
		if (!sse41 && !avx2 && !avx512)
        {
            if (vi.ComponentSize() == 1)
                resize_plane_c<uint8_t>(out[i_coeff],srcp,dstp,dst_width,dst_height,src_stride,dst_stride,ValMin[i],ValMax[i]);
            else if (vi.ComponentSize() == 2)
                resize_plane_c<uint16_t>(out[i_coeff],srcp,dstp,dst_width,dst_height,src_stride,dst_stride,ValMin[i],ValMax[i]);
            else
                resize_plane_c<float>(out[i_coeff],srcp,dstp,dst_width,dst_height,src_stride,dst_stride,ValMin[i],ValMax[i]);
        }
        else
            process_frame(out[i], srcp, dstp, dst_width, dst_height, src_stride, dst_stride);*/
    }

    return dst;
}

AVSValue __cdecl Create_JincResize(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    const VideoInfo& vi = args[0].AsClip()->GetVideoInfo();

    return new JincResize(
        args[0].AsClip(),
        args[1].AsInt(),
        args[2].AsInt(),
        args[3].AsFloat(0),
        args[4].AsFloat(0),
        args[5].AsFloat(static_cast<float>(vi.width)),
        args[6].AsFloat(static_cast<float>(vi.height)),
        args[7].AsInt(256),
        args[8].AsInt(256),
        args[9].AsInt(3),
        args[10].AsFloat(1.0),
        args[11].AsInt(-1),
        env);
}

/*
template <int taps>
AVSValue __cdecl resizer_jinc36resize(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    Arguments mapped_args;

    resizer(args, &mapped_args);
    mapped_args.add(args[9].AsInt(taps), "tap");

    return env->Invoke("JincResizeMT", mapped_args.args(), mapped_args.arg_names()).AsClip();
}

const AVS_Linkage* AVS_linkage;
*/

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

	poolInterface=ThreadPoolInterface::Init(0);

	if (!poolInterface->GetThreadPoolInterfaceStatus()) env->ThrowError("JincResizeMT: Error with the TheadPool status!");

    env->AddFunction("JincResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i[tap]i[blur]f[opt]i" \
		"[range]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_JincResize, 0);
/*
    env->AddFunction("Jinc36ResizeMT", "cii[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i" \
		"[range]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", resizer_jinc36resize<3>, 0);
    env->AddFunction("Jinc64ResizeMT", "cii[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i" \
		"[range]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", resizer_jinc36resize<4>, 0);
    env->AddFunction("Jinc144ResizeMT", "cii[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i" \
		"[range]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", resizer_jinc36resize<6>, 0);
    env->AddFunction("Jinc256ResizeMT", "cii[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i" \
		"[range]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", resizer_jinc36resize<8>, 0);
*/
    return JINCRESIZEMT_VERSION;
}
