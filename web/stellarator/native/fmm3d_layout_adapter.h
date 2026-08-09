#ifndef BIESOLVER_FMM3D_LAYOUT_ADAPTER_H
#define BIESOLVER_FMM3D_LAYOUT_ADAPTER_H

#include "lfortran_array.h"

#include <stdint.h>

void biesolver_lfmm3d_t_cd_p_rowmajor(
    double eps, int64_t nsource, struct r64 *source,
    struct r64 *charge, struct r64 *dipvec,
    int64_t ntarg, struct r64 *targ, struct r64 *pottarg,
    int64_t *ier, int64_t *adapter_status, double *elapsed_seconds);

#endif
