/* FKA121 - Computational Physics - Exercise 3 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.141592653589
#define nbr_points 10000000
#define nbr_dim 3

double calculate_weight(double coords[nbr_dim]);
double calculate_function(double p[nbr_dim]);
double calculate_function_g(double p[nbr_dim]);


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
  double delta = 2.5;
  size_t nbr_switching_state = 0;

  const gsl_rng_type *T; /* static info about rngs */
  gsl_rng *q; /* rng instance */
  gsl_rng_env_setup(); /* setup the rngs */
  T = gsl_rng_default; /* specify default rng */
  q = gsl_rng_alloc(T); /* allocate default rng */
  gsl_rng_set(q, time(NULL)); /* Initialize rng */


  double *function_val = malloc(sizeof(double) * nbr_points);

  //initialize random coordinates
  for(j = 0; j < nbr_dim; j++){
    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    n1[j] = rand_nbr;

    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    n2[j] = rand_nbr;
  }

  p_m = calculate_weight(m);

  sum_tmp = 0;
  // Main loop
  for(i = 0; i < nbr_points; i++){
    for(j = 0; j < nbr_dim; j++){
      rand_nbr = gsl_rng_uniform(q); /*generate random number 0-1 (repeatable)*/
      n1[j] = m1[j] + delta*(rand_nbr - 0.5);
      rand_nbr = gsl_rng_uniform(q); /*generate random number 0-1 (repeatable)*/
      n2[j] = m2[j] + delta*(rand_nbr - 0.5);
    }

    p_m = calculate_weight(m);
    p_n = calculate_weight(n);

    rand_nbr = gsl_rng_uniform(q);
    if((p_n / p_m) > rand_nbr){
      m1[0] = n1[0];
      m1[1] = n1[1];
      m1[2] = n1[2];

      m2[0] = n2[0];
      m2[1] = n2[1];
      m2[2] = n2[2];
     

      function_val[i] = calculate_function_g(n);
      nbr_switching_state++;
    } else {
      function_val[i] = calculate_function_g(m);
    }
  }// End main loop

  sum_tmp = 0;
  for(i = 0; i < nbr_points; i++){
    sum_tmp += function_val[i];
  }
  I_value = sum_tmp / nbr_points;

  // Calculate variance
  sum_tmp = 0;
  for(i = 0; i < nbr_points; i++){
    sum_tmp += function_val[i] * function_val[i];
  }
  variance = sum_tmp / nbr_points - I_value * I_value;
  error = sqrt(variance / nbr_points);

  // Deallocate rng
  gsl_rng_free (q);

  printf("Integral: %e, with error(+-): %e\n", I_value, error);
  printf("nbr_switching_state: %lf\n", (double) nbr_switching_state / nbr_points);

}

double calculate_function(double p[nbr_dim]){
  double tmp1, tmp2, fvalue;
  tmp1 = p[0]*p[0] + p[0]*p[0]*p[1]*p[1] + p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
  tmp2 = exp(-(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
  fvalue = pow((double) PI, -3.0/2.0) * tmp1 * tmp2;
  return fvalue;
}

double calculate_function_g(double p[nbr_dim]){
  double fvalue;
  fvalue = p[0]*p[0] + p[0]*p[0]*p[1]*p[1] + p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
  return fvalue;
}

double calculate_weight(double coords[3]){
  double prob;
  prob = pow((double) PI, -3.0/2.0) *
    exp(-(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]));
    return prob;
}

double local_energy (double m1[nbr_dim], double m2[nbr_dim], double alpha){
  double E_l;

  

}



// // Setup random number generator
// double rand;
// const gsl_rng_type *T; /* static info about rngs */
// gsl_rng *q; /* rng instance */
// gsl_rng_env_setup(); /* setup the rngs */
// T = gsl_rng_default; /* specify default rng */
// q = gsl_rng_alloc(T); /* allocate default rng */
// gsl_rng_set(q,time(NULL)); /* Initialize rng */
// rand = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
// // Deallocate rng
// gsl_rng_free (q);
