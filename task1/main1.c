/* FKA121 - Computational Physics - Exercise 3 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.141592653589
#define nbr_iterations 1000000
#define nbr_dim 3

double calculate_wave_function(double[nbr_dim], double[nbr_dim], double);
double calculate_probability(double[nbr_dim], double[nbr_dim], double);
double local_energy (double[nbr_dim], double[nbr_dim], double);
double pow(double x, double y);
double calculate_angle(double[nbr_dim], double[nbr_dim]);


/* Main Program */
int main(){

  size_t i,j;
  double sum_tmp;
  double rand_nbr;
  double m1[nbr_dim] = {0};
  double m2[nbr_dim] = {0};
  double n1[nbr_dim] = {0};
  double n2[nbr_dim] = {0};
  double p_m, p_n;
  double I_value;
  double variance;
  double error;
  double delta = 1;
  double alpha = 0.1;

  size_t nbr_switching_state = 0;

  const gsl_rng_type *T; /* static info about rngs */
  gsl_rng *q; /* rng instance */
  gsl_rng_env_setup(); /* setup the rngs */
  T = gsl_rng_default; /* specify default rng */
  q = gsl_rng_alloc(T); /* allocate default rng */
  gsl_rng_set(q, time(NULL)); /* Initialize rng */


  double *energy = malloc(sizeof(double) * nbr_iterations);
  double *electron1_distance = malloc(sizeof(double) * nbr_iterations);
  double *electron2_distance = malloc(sizeof(double) * nbr_iterations);
  double *theta = malloc(sizeof(double) * nbr_iterations);

  //initialize random coordinates
  for(j = 0; j < nbr_dim; j++){
    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    m1[j] = rand_nbr;

    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    m2[j] = rand_nbr;
  }

  p_m = calculate_probability(m1, m2, alpha);
  sum_tmp = 0;
  // Main loop
  for(i = 0; i < nbr_iterations; i++){

    for(j = 0; j < nbr_dim; j++){
      rand_nbr = gsl_rng_uniform(q); /*generate random number 0-1 (repeatable)*/
      n1[j] = m1[j] + delta*(rand_nbr - 0.5);
      rand_nbr = gsl_rng_uniform(q); /*generate random number 0-1 (repeatable)*/
      n2[j] = m2[j] + delta*(rand_nbr - 0.5);
    }

    p_m = calculate_probability(m1,m2, alpha);
    p_n = calculate_probability(n1,n2, alpha);

    rand_nbr = gsl_rng_uniform(q);
    if((p_n / p_m) > rand_nbr){
      m1[0] = n1[0];
      m1[1] = n1[1];
      m1[2] = n1[2];

      m2[0] = n2[0];
      m2[1] = n2[1];
      m2[2] = n2[2];


      energy[i] = local_energy(m1, m2, alpha);
      nbr_switching_state++;
    } else {
      energy[i] = local_energy(m1, m2, alpha);
    }
    //save electron distance from origo
    electron1_distance[i] = sqrt(m1[0]*m1[0] + m1[1]*m1[1] + m1[2]*m1[2]);
    electron2_distance[i] = sqrt(m2[0]*m2[0] + m2[1]*m2[1] + m2[2]*m2[2]);

    //save angle
    theta[i] = calculate_angle(m1,m2);

  }// End main loop

  printf("Precent switched states:%f\n", nbr_switching_state / (double) nbr_iterations);
  sum_tmp = 0;
  for(i = 0; i < nbr_iterations; i++){
    sum_tmp += energy[i];
  }
  I_value = sum_tmp / nbr_iterations;


  printf("Avg energy =%f \nNumber of iterations=%d\n", I_value, nbr_iterations);
  // Calculate variance
  sum_tmp = 0;
  for(i = 0; i < nbr_iterations; i++){
    sum_tmp += energy[i] * energy[i];
  }
  variance = sum_tmp / nbr_iterations - I_value * I_value;
  error = sqrt(variance / nbr_iterations);

  FILE *fp;

  //Write to file, electron distance from origo data
  fp = fopen("electron_dist.dat","w");
  for (i = 0; i < nbr_iterations; i++){
    fprintf(fp, "%e \t %e", electron1_distance[i], electron2_distance[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);

  //Write to file, theta distribution data
  fp = fopen("theta_dist.dat","w");
  for (i = 0; i < nbr_iterations; i++){
    fprintf(fp, "%e", theta[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);

  // Deallocate rng
  gsl_rng_free (q);

}//End MAIN


double calculate_wave_function(double r1[nbr_dim], double r2[nbr_dim], double alpha){
  double norm_r1, norm_r2, norm_r1r2, value;

  norm_r1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
  norm_r2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
  norm_r1r2 =sqrt((r1[0]-r2[0])*(r1[0]-r2[0]) +
                 (r1[1]-r2[1])*(r1[1]-r2[1]) + (r1[2]-r2[2])*(r1[2]-r2[2]));

  value = exp(-2*norm_r1) * exp(-2*norm_r2) *
    exp(norm_r1r2 / (2 * (1+alpha*norm_r1r2)));
  return value;
}


double calculate_probability(double r1[nbr_dim], double r2[nbr_dim], double alpha){
  double value;
  value = calculate_wave_function(r1, r2, alpha);
  return value*value;
}


double local_energy (double m1[nbr_dim], double m2[nbr_dim], double alpha){
  double E_l;
  double length_m1;
  double length_m2;
  double cross_mult;
  double m1_squared, m2_squared;
  double length_m12;
  double a_r12;

  length_m12 = sqrt((m1[0]-m2[0])*(m1[0]-m2[0]) + (m1[1]-m2[1])*(m1[1]-m2[1]) + (m1[2]-m2[2])*(m1[2]-m2[2]));
  length_m1 = sqrt((m1[0]*m1[0])+(m1[1]*m1[1])+(m1[2]*m1[2]));
  length_m2 = sqrt((m2[0]*m2[0])+(m2[1]*m2[1])+(m2[2]*m2[2]));

  cross_mult = m1[0]*m2[0]+m1[1]*m2[1]+m1[2]*m2[2];
  m1_squared = ((m1[0]*m1[0])+(m1[1]*m1[1])+(m1[2]*m1[2]));
  m2_squared = ((m2[0]*m2[0])+(m2[1]*m2[1])+(m2[2]*m2[2]));
  a_r12 = (1.0+alpha*length_m12);

  E_l = -4.0 + (m1_squared / length_m1 - cross_mult / length_m1 - cross_mult / length_m2 + m2_squared / length_m2) / (length_m12*pow(a_r12,2)) - 1.0/ (length_m12*pow(a_r12,3)) - 0.25/ (pow(a_r12, 4))+ 1.0 / length_m12;

  return E_l;
}


double calculate_angle(double r1[nbr_dim], double r2[nbr_dim]){
  double norm_r1, norm_r2, prod, theta;

  norm_r1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
  norm_r2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);

  prod = (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]);

  theta = acos(prod / (norm_r1 * norm_r2));

  return theta;
}
