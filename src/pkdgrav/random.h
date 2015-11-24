#ifndef RANDOM_H
#define RANDOM_H

double randUniform(void);
void randSeedGenerator(int iProcID);
double randGaussian(void);
double randPoisson(double dMean);

#endif /* RANDOM_H */
