Internaly multi-threaded JincResize, based on v2.1.4.

===========================================================

The parameters of original JincResize v2.1.4 are kept with the same name, the same order
and with the same feature, except "threads" which changes just a little to match the behavior
of my others MT plugins, but "0" and "1" values still have the same effect, so scripts using
JincResize (v2.1.4) are fully compatible with JincResizeMT.

Several new parameters are added at the end of all the parameters.

JincResizeMT(clip, int target_width, int target_height, float src_left, float src_top, float src_width,
  float src_height, int quant_x, int quant_y, int tap, float blur, string cplace, int threads,
  int opt, int initial_capacity, float initial_factor,
  int range, bool logicalCores, bool MaxPhysCore, bool SetAffinity, bool sleep, int prefetch, int ThreadLevel)
)

------------------------------

##### There are 4 additional functions:
    Jinc36ResizeMT is an alias for JincResizeMT(tap=3).
    Jinc64ResizeMT is an alias for JincResizeMT(tap=4).
    Jinc144ResizeMT is an alias for JincResizeMT(tap=6).
    Jinc256ResizeMT is an alias for JincResizeMT(tap=8).

Jinc36ResizeMT/Jinc64ResizeMT/Jinc144ResizeMT/Jinc256ResizeMT(clip, int target_width, int target_height,
  float src_left, float src_top, float src_width, float src_height, int quant_x, int quant_y,
  string cplace, int threads,
  int range, bool logicalCores, bool MaxPhysCore, bool SetAffinity, bool sleep, int prefetch, int ThreadLevel)

------------------------------

   clip -
      A clip to process. All planar (and only planar) formats are supported.

   target_width -
      The width of the output.

   target_height -
      The height of the output.

   src_left -
      Cropping of the left edge.

      Default: 0.0.

   src_top -
      Cropping of the top edge.

      Default: 0.0.

   src_width -
      If > 0.0 it sets the width of the clip before resizing.
      If <= 0.0 it sets the cropping of the right edges before resizing.

      Default: Source width.

   src_height -
      If > 0.0 it sets the height of the clip before resizing.
      If <= 0.0 it sets the cropping of the bottom edges before resizing.

      Default: Source height.

   quant_x, quant_y -
      Controls the sub-pixel quantization.
      Must be between 1 and 256.
      Default: 256.

   tap (JincResizeMT only) -
      Corresponding to different zero points of Jinc function.
      Must be between 1 and 16.

      Default: 3.

   blur (JincResizeMT only) -
      Blur processing, it can reduce side effects.
      To achieve blur, the value should be less than 1.0.

      Default: 1.0.

   threads -
      Controls how many threads will be used for processing. If set to 0, threads will
      be set equal to the number of detected logical or physical cores,according logicalCores parameter.

      Default:  0  (int)

   cplace -
      The location of the chroma samples. Checked and used only on subsampled formats.
      "MPEG1": Chroma samples are located on the center of each group of 4 pixels.
      "MPEG2": Chroma samples are located on the left pixel column of the group.
      "topleft" or "top_left": Chroma samples are located on the left pixel column
        and the first row of the group.
      Empty or "auto": Use frame property. If frame properties are supported and frame property
        "_ChromaLocation" exists - "_ChromaLocation" value of the first frame is used.
        If frame properties aren't supported or there is no property "_ChromaLocation" - "MPEG2".

      Default: "auto"

   opt (JincResizeMT only) -
      Sets which cpu optimizations to use.
      -1: Auto-detect without AVX-512.
      0: Use C++ code.
      1: Use SSE4.1 code.
      2: Use AVX2 code.
      3: Use AVX-512 code.

      Default: -1.

   initial_capacity (JincResizeMT only) -
      Initial memory allocation size.
      Lower size forces more further memory reallocating that leads to initial slower startup but
        avoids excessive memory allocation.
      Must be greater than 0.

      Default: Max(target_width * target_height, src_width * src_height).

   initial_factor (JincResizeMT only) -
      The initial factor used for the first memory reallocation.
      After the first memory reallocation the factor starts to lower for the next reallocations.
      "initial_factor=1" ensures that the next memory allocation is the minimal possible.
      Must be equal to or greater than 1.0.

      Default: 1.5.

  range -
      This parameter specify the range the output video data has to comply with.
      Limited range is 16-235 for Y, 16-240 for U/V. Full range is 0-255 for all planes.
      Alpha channel is not affected by this paramter, it's always full range.
      Values are adjusted according bit depth of course. This parameter has no effect
      for float datas.
      0 : Automatic mode. If video is YUV mode is limited range, if video is RGB mode is
          full range, if video is greyscale (Y/Y8) mode is Y limited range.
      1 : Force full range whatever the video is.
      2 : Force limited Y range for greyscale video (Y/Y8), limited range for YUV video,
          no effect for RGB video.
      3 : Force limited U/V range for greyscale video (Y/Y8), limited range for YUV video,
          no effect for RGB video.
      4 : Force special camera range (16-255) for greyscale video (Y/Y8) and YUV video,
          no effect for RGB video.

      Default: 1

   logicalCores -
      If threads is set to 0, it will specify if the number of threads will be the number
      of logical CPU (true) or the number of physical cores (false). If your processor doesn't
      have hyper-threading or threads<>0, this parameter has no effect.

      Default: true (bool)

   MaxPhysCore -
      If true, the threads repartition will use the maximum of physical cores possible. If your
      processor doesn't have hyper-threading or the SetAffinity parameter is set to false,
      this parameter has no effect.

      Default: true (bool)

   SetAffinity -
      If this parameter is set to true, the pool of threads will set each thread to a specific core,
      according MaxPhysCore parameter. If set to false, it's leaved to the OS.
      If prefecth>number of physical cores, it's automaticaly set to false.

      Default: false (bool)

  sleep -
      If this parameter is set to true, once the filter has finished one frame, the threads of the
      threadpool will be suspended (instead of still running but waiting an event), and resume when
      the next frame will be processed. If set to false, the threads of the threadpool are always
      running and waiting for a start event even between frames.

      Default: false (bool)

  prefetch - (added negative trim feature)
      This parameter will allow to create more than one threadpool, to avoid mutual resources acces
      if "prefetch" is used in the avs script.
      0 : Will set automaticaly to the prefetch value use in the script. Well... that's what i wanted
          to do, but for now it's not possible for me to get this information when i need it, so, for
          now, 0 will result in 1. For now, if you're using "prefetch" in your script, put the same
          value on this parameter.

      Default: 0

  ThreadLevel -
      This parameter will set the priority level of the threads created for the processing (internal
      multithreading). No effect if threads=1.
      1 : Idle level.
      2 : Lowest level.
      3 : Below level.
      4 : Normal level.
      5 : Above level.
      6 : Highest level.
      7 : Time critical level (WARNING !!! use this level at your own risk)

      Default : 6

The logicalCores, MaxPhysCore, SetAffinity and sleep are parameters to specify how the pool of thread
will be created and handled, allowing if necessary each people to tune according his configuration.

Notes:
AVX-512 benchmark were worse than AVX-2 on my config, this is why for now it's not on auto-detect.
