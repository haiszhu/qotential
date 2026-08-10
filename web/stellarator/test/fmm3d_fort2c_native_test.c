#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void lfmm3d_t_cd_p_(double *eps, int64_t *nsource, double *source,
                     double *charge, double *dipvec, int64_t *ntarg,
                     double *targ, double *pottarg, int64_t *ier);

static uint64_t random_state = UINT64_C(0x4d595df4d0f33173);

static double random_unit(void)
{
    random_state = random_state*UINT64_C(6364136223846793005) + 1;
    return (double)(random_state >> 11)*0x1.0p-53;
}

int main(void)
{
    const double inv4pi = 0.07957747154594766788444188168626;
    const double threshold2 = 1.0e-30;
    int64_t nsource = 10000, ntarg = nsource, ier = -1;
    double eps = 1.0e-12;
    double *source = calloc((size_t)3*nsource, sizeof(double));
    double *charge = calloc((size_t)nsource, sizeof(double));
    double *dipvec = calloc((size_t)3*nsource, sizeof(double));
    double *fmm = calloc((size_t)ntarg, sizeof(double));
    double *direct = calloc((size_t)ntarg, sizeof(double));
    double numerator = 0.0, denominator = 0.0, max_abs = 0.0;
    double max_abs_direct = 0.0;

    if (!source || !charge || !dipvec || !fmm || !direct) return 2;
    for (int64_t j = 0; j < nsource; ++j) {
        source[3*j] = random_unit();
        source[3*j+1] = random_unit();
        source[3*j+2] = random_unit();
        charge[j] = random_unit()-0.5;
        dipvec[3*j] = 0.1*(random_unit()-0.5);
        dipvec[3*j+1] = 0.1*(random_unit()-0.5);
        dipvec[3*j+2] = 0.1*(random_unit()-0.5);
    }

    lfmm3d_t_cd_p_(&eps, &nsource, source, charge, dipvec,
                    &ntarg, source, fmm, &ier);
    if (ier != 0) {
        fprintf(stderr, "fort2c FMM returned ier=%lld\n", (long long)ier);
        return 3;
    }

    for (int64_t i = 0; i < ntarg; ++i) {
        double value = 0.0;
        for (int64_t j = 0; j < nsource; ++j) {
            const double dx = source[3*i]-source[3*j];
            const double dy = source[3*i+1]-source[3*j+1];
            const double dz = source[3*i+2]-source[3*j+2];
            const double r2 = dx*dx+dy*dy+dz*dz;
            if (r2 <= threshold2) continue;
            const double rinv = 1.0/sqrt(r2);
            const double cd = inv4pi*rinv;
            const double dot = dx*dipvec[3*j] +
                               dy*dipvec[3*j+1] + dz*dipvec[3*j+2];
            value += charge[j]*cd + dot*cd/r2;
        }
        direct[i] = value;
        if (fabs(value) > max_abs_direct) max_abs_direct = fabs(value);
        if (!isfinite(fmm[i])) {
            fprintf(stderr, "non-finite fort2c FMM output at %lld\n",
                    (long long)i);
            return 4;
        }
        const double delta = fmm[i]-direct[i];
        numerator += delta*delta;
        denominator += direct[i]*direct[i];
        if (fabs(delta) > max_abs) max_abs = fabs(delta);
    }

    const double relative_l2 = sqrt(numerator/denominator);
    const double normalized_max = max_abs/fmax(1.0, max_abs_direct);
    printf("FMM3D_FORT2C_NATIVE_CLANG n=%lld eps=%.1e rel_l2=%.17e "
           "max_abs=%.17e normalized_max=%.17e\n",
           (long long)nsource, eps, relative_l2, max_abs, normalized_max);
    free(direct);
    free(fmm);
    free(dipvec);
    free(charge);
    free(source);
    return normalized_max <= 2.0e-11 ? 0 : 5;
}
