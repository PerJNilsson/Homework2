/* FKA121 - Computational Physics - Exercise 3 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.141592653589
#define nbr_iterations 1E7
#define nbr_dim 3

double calculate_wave_function(double[nbr_dim], double[nbr_dim], double);
double calculate_probability(double[nbr_dim], double[nbr_dim], double);
double local_energy (double[nbr_dim], double[nbr_dim], double);
double pow(double x, double y);
double calculate_angle(double[nbr_dim], double[nbr_dim]);

int auto_corr_fun(double * N, int nbr_of_lines, int k);
double block_average(double* data, int nbr_of_lines, int B);
void center_data(double* N, int nbr_of_lines);


/* Main Program */
int main(){

  size_t i,j;
  double sum_tmp, sum_squared_tmp;
  double rand_nbr;
  double m1[nbr_dim] = {0};
  double m2[nbr_dim] = {0};
  double n1[nbr_dim] = {0};
  double n2[nbr_dim] = {0};
  double p_m, p_n;
  double I_value;
  double variance;
  double error;
  double delta = 0.975;
  double alpha = 0.1434;
  int equi_phase = 2000;

  size_t nbr_switching_state = 0;

  const gsl_rng_type *T; /* static info about rngs */
  gsl_rng *q; /* rng instance */
  gsl_rng_env_setup(); /* setup the rngs */
  T = gsl_rng_default; /* specify default rng */
  q = gsl_rng_alloc(T); /* allocate default rng */
  gsl_rng_set(q, time(NULL)); /* Initialize rng */


  double *energy = malloc(sizeof(double) * nbr_iterations);
  double *theta = malloc(sizeof(double) * nbr_iterations);

  int nbr_average_runs = 100;
  double avg_energy_all_runs = 0;
  double variance_all_runs = 0;
  double save_I_value[nbr_average_runs];

  for (int kx=0; kx<nbr_average_runs; kx ++){


  //initialize random coordinates
  for(j = 0; j < nbr_dim; j++){
    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    m1[j] = 5;

    rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
    rand_nbr -= 0.5;
    m2[j] = 5;
  }

  p_m = calculate_probability(m1, m2, alpha);
  sum_tmp = 0;
  
  // Equilibrium phase
  for (int ix=0; ix<equi_phase; ix++){
    for(j = 0; j < nbr_dim; j++){
      rand_nbr = gsl_rng_uniform(q); //generate random number 0-1 (repeatable)
      n1[j] = m1[j] + delta*(rand_nbr - 0.5);
      rand_nbr = gsl_rng_uniform(q); //generate random number 0-1 (repeatable)
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
    }
  }
  
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

  }// End main loop

  sum_tmp = 0;
  sum_squared_tmp = 0;
  for(i = 0; i < nbr_iterations; i++){
    sum_tmp += energy[i];
  }
  I_value = sum_tmp / nbr_iterations;
  save_I_value[kx] = I_value;
  printf("%d\n", kx);
  } // End the overall average loop

  double sigma_all2; 
  // Calculate variance
  double tmp_sum_var = 0;
  double tmp_sum_avg = 0;
  for (i=0; i<nbr_average_runs; i++){
    tmp_sum_var += save_I_value[i]*save_I_value[i] / (double)nbr_average_runs;
    tmp_sum_avg += save_I_value[i] / (double) nbr_average_runs;
  }
  sigma_all2 = sqrt(tmp_sum_var - tmp_sum_avg*tmp_sum_avg);

  printf("only I: %f\n", sigma_all2);
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

  E_l = -4.0 + (m1_squared / length_m1 - cross_mult / length_m1 - cross_mult / length_m2 + m2_squared / length_m2) / (length_m12*pow(a_r12,2)) - 1.0/ (length_m12*pow(a_r12,3)) - 0.25/(pow(a_r12, 4))+ 1.0 / length_m12;

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


int auto_corr_fun(double* data, int nbr_of_lines, int k){
  center_data(data, nbr_of_lines);
  double mean_ik[k];
  double phi[k];
  int good_guess=0;
  double s = exp(-2.0);
  for (int i=0; i<k; i++){
    mean_ik[i]=0;
    phi[i]=0;
  }
  double mean_squared = 0;
  FILE * valuessave;
  for (int i = 0; i<k; i++) {
    for (int j=0; j<nbr_of_lines-i; j++){
      mean_ik[i] += data[j]*data[i+j]/(double)(nbr_of_lines-i);
    }
  }
  for (int i=0; i<nbr_of_lines; i++){
    mean_squared += data[i]*data[i]/(double) nbr_of_lines;
  }
  valuessave = fopen("values.dat", "w");
  for (int i=0; i<k; i++){
    phi[i] = (mean_ik[i]) / (mean_squared);
    fprintf(valuessave, "%d \t %f\n", i, phi[i]);
    if(sqrt((phi[good_guess]-s)*(phi[good_guess]-s))> sqrt((phi[i]-s)*(phi[i]-s))){
      good_guess = i;
    }
  }
  fclose(valuessave);
  return good_guess;
}

double block_average(double* data, int nbr_of_lines, int B){
  center_data(data, nbr_of_lines);
  double good_guess_B;
  int j;
  j = nbr_of_lines / B;
  double F[j];
  for (int i=0; i<j; i++){
    F[i] = 0;
  }

  for (int i=0; i<j; i++){
    for (int k=0; k<B; k++){
      F[i] += data[k+i*B]/(double)B;
    }
  }
  double variance_F = 0;
  double mean_F = 0;
  for (int i=0; i<j; i++){
    variance_F += F[i]*F[i]/(double)j;
    mean_F += F[i]/(double)j;
  }
  variance_F = variance_F - mean_F*mean_F;

  double variance_f = 0;
  for (int i=0; i<nbr_of_lines; i++){
    variance_f += data[i]*data[i];
  }
  variance_f = variance_f / (double)nbr_of_lines;
  good_guess_B = (double)B * variance_F / variance_f;
  return good_guess_B;
}


void center_data(double*data, int nbr_of_lines){
  double main_mean=0;
  for (int i=0; i<nbr_of_lines; i++){
    main_mean += data[i];
  }
  main_mean = main_mean / (double)nbr_of_lines;
  for (int i=0; i<nbr_of_lines; i++){
    data[i]=data[i]-main_mean;
  }
}
