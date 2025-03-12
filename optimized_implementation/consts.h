#ifndef CONSTS_H
#define CONSTS_H

#define _8XQ            0
#define _8XQINV         8
#define _8XWmont        16
#define _8XWmont_QINV   24
#define _8XF1           32
#define _8XF1_QINV      40
#define _8XF2           48
#define _8XF2_QINV      56
#define _ZETAS_QINV     64
#define _8idx           64+3328+3328+3327+3327

#if NIMS_TRI_NTT_MODE == 1
#define _ZETAS          64+2176 // 2176 = 128 + 2048
#define _INVZETAS_QINV  64+2176+2176
#define _INVZETAS       64+2176+2176+2176

#elif NIMS_TRI_NTT_MODE == 3
#define _ZETAS          64+512
#define _INVZETAS_QINV  64+512+1200
#define _INVZETAS       64+512+1200+512

#elif NIMS_TRI_NTT_MODE == 5
#define _ZETAS          64+3328
#define _INVZETAS_QINV  64+3328+3328
#define _INVZETAS       64+3328+3328+3327

#elif NIMS_TRI_NTT_MODE == 55
#define _ZETAS          64+4352
#define _INVZETAS_QINV  64+4352+4352
#define _INVZETAS       64+4352+4352+4096+192+56+8
#endif

#ifndef __ASSEMBLER__

#include "align.h"
typedef ALIGNED_INT32(64+N+N+N+N+N+N+N+N) qdata_t;

extern const qdata_t qdata;

#endif
#endif
