#ifndef SPRINGS_HINCLUDED
#define SPRINGS_HINCLUDED

#define MAX_NUM_SPRINGS_PER_PARTICLE 32 /* must not be < nSmooth */

typedef struct {
  int iOrder;
  float fZeroStrainLength;
  float fYoungsModulus;
  float fStressLimit;
  } SPRING;

typedef struct {
  double dMeanYoungsModulus;
  double dMeanStrainLimit;
  double dYoungsStdDev;
  double dStrainStdDev;
  double dMaxStrainLimit;
  double dLinkageLength;
  double dZeroStrainLength;
  double dDamp;
  double dPorosityFactor;
  int bReadSpringsData;
  } SPRING_PARAMS;

#define SPRHEAD_SIZE SSHEAD_SIZE /* same struct */

#define SPRDATA_SIZE (8*(MAX_NUM_SPRINGS_PER_PARTICLE))

#define SPRINGS_STRESS_NONE 0
#define SPRINGS_STRESS_MOVE_TENSILE 1
#define SPRINGS_STRESS_MOVE_SHEAR 2
#define SPRINGS_STRESS_FORCE_TENSILE 3
#define SPRINGS_STRESS_FORCE_SHEAR 4
#define SPRINGS_STRESS_FORCE_TWIST 5
#define SPRINGS_STRESS_FORCE_RADIAL 6

#define SPRINGS_STRESS_TEST 0

#if ((SPRINGS_STRESS_TEST) >= 1 && (SPRINGS_STRESS_TEST) <= 2)
#define SPRINGS_STRESS_MOVE
#endif /* !((SPRINGS_STRESS_TEST) >= 1 && (SPRINGS_STRESS_TEST) <= 2) */

#define REFORM_SPRINGS_NONE 0
#define REFORM_SPRINGS_ALLPARTICLES 1
#define REFORM_SPRINGS_NOWALLPARTICLES 2 /* exclude particles stuck to walls (iColor < 0) */
#define REFORM_SPRINGS_NOLOOSEPARTICLES 3 /* do not allow non-wall particles to stick to other non-wall particles */

#define REFORM_SPRINGS 0

#if (REFORM_SPRINGS != 0)
#define ZERO_STRAIN_LENGTH (2.)
#endif /* REFORM_SPRINGS */

/* #define YOUNGS_INCREASE */ /* HACK to allow Young's modulus to increase in time (only functional when REFORM_SPRINGS is defined) */

#ifdef YOUNGS_INCREASE
#define YOUNGS_INCREASE_PER_STEP (0.0000002) /* at each step, Y will increase by dMeanYoungsModulus multiplied by this value */
#endif /* YOUNGS_INCREASE */

/* #define SPRINGS_COLOR_GREYSCALE */ /* alternate coloring of springs */

#endif /* !SPRINGS_HINCLUDED */
