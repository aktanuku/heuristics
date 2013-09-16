
#ifndef SMOOTHPEAK_H
#define SMOOTHPEAK_H
#include"init_fp.h"

#define PI 3.14159265358979323846264

extern void smoothPeak(float aPeak, float aRadius, float aDelta, int aNumBins, fingerprint* F,int row, float aWidth);
extern double gaussianIntegral(double x);

#endif

