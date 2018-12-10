/* FKA121 - Computational Physics - Exercise 3 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.141592653589
#define nbr_iterations 100000
#define nbr_dim 3

double calculate_wave_function(double[nbr_dim], double[nbr_dim], double);
double calculate_probability(double[nbr_dim], double[nbr_dim], double);
double local_energy (double[nbr_dim], double[nbr_dim], double);
double pow(double x, double y);
double calculate_angle(double[nbr_dim], double[nbr_dim]);

int auto_corr_fun(double * N, int nbr_of_lines, int k);
double block_average(double* data, int nbr_of_lines, int B);
void center_data(double* N, int nbr_of_lines);


double gradient_E_alpha(double * energy, double m1[], double m2[], double alpha, int length_E, double * ln_wavefun);
void get_ln_wavefun(double * ln_wavefun, double m1[], double m2[], double alpha, int iteration);
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
  double delta = 0.975;
  double alpha = 0.1;
  int equi_phase = 1500;
  size_t different_alphas= 1000;
  size_t nbr_switching_state = 0;
  size_t different_beta = 6;
  double gamma_alpha;
  double beta_alpha; // 0.5 < beta <= 1.0
  double grad_E_a;
  size_t iterations_per_alpha = 1;
  double sum_avg_alpha = 0;

  double save_beta[different_beta];

  double * pointer_save = malloc(sizeof(double)*different_alphas*different_beta);
  double ** matrix_save = malloc(sizeof(double)*different_alphas);
  for(i=0, j=0; i<different_alphas; i++, j+=different_beta){
    matrix_save[i] = pointer_save +j;
  }

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
  double *alpha_energy_values = malloc(sizeof(double)*different_alphas);

  double * alpha_values = malloc(sizeof(double)*different_alphas);
  double * ln_wavefun = malloc(sizeof(double)*nbr_iterations);

  alpha = 0.25;
  save_beta[0] = 0.51; save_beta[1] = 0.60; save_beta[2] =0.70; save_beta[3] = 0.80;save_beta[4] = 0.90; save_beta[5] = 1.0;
  alpha_values[0] = alpha;
  for (i=0;i<different_beta;i++){
    matrix_save[0][i] = save_beta[i];
  }

  for (size_t kx=0; kx<different_beta; kx++) {
    beta_alpha = save_beta[kx];
    for (size_t jx=1; jx<different_alphas; jx++){
      nbr_switching_state = 0;
      //initialize random coordinates
      for(j = 1; j < nbr_dim; j++){
        rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
        rand_nbr -= 0.5;
        m1[j] = 50;

        rand_nbr = gsl_rng_uniform(q); /* generate random number 0-1 (repeatable) */
        rand_nbr -= 0.5;
        m2[j] = 50;
      }

      p_m = calculate_probability(m1, m2, alpha);
      sum_tmp = 0;

      // Equilibrium phase
      for (int ix=0; ix<equi_phase; ix++){
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
          get_ln_wavefun(ln_wavefun, m1, m2, alpha, i);
          nbr_switching_state++;
        } else {
          energy[i] = local_energy(m1, m2, alpha);
          get_ln_wavefun(ln_wavefun, m1, m2, alpha, i);
        }

      }// End main loop;
      gamma_alpha = pow(jx,-1.0*beta_alpha);
      grad_E_a = gradient_E_alpha(energy, m1, m2, alpha, nbr_iterations, ln_wavefun);
      alpha = alpha_values[jx-1] - gamma_alpha*grad_E_a;
      alpha_values[jx] = alpha;
      matrix_save[jx][kx] = alpha;
      sum_tmp = 0;
      for(i = 0; i < nbr_iterations; i++){
        sum_tmp += energy[i];
      }
      I_value = sum_tmp / nbr_iterations;



      alpha_energy_values[jx] += I_value / iterations_per_alpha;
      //printf("Avg energy =%f \n", I_value);
      //printf("Precent switched states:%f\n", nbr_switching_state / (double) nbr_iterations);
  }//end different alphas, jx
  printf("alpha=%f\n", alpha);
  sum_avg_alpha +=alpha / (double)different_beta;
  save_beta[kx] = beta_alpha;
  } // end different betas
  printf("Final alpha value=%f\n", sum_avg_alpha);
  // Saving results to disk


  FILE * save_beta_alpha;
  save_beta_alpha = fopen("beta_alpha.dat", "w");
  for (i=0; i<different_alphas; i++) {
    for (j=0; j<different_beta; j++){
      fprintf(save_beta_alpha, "%f \t", matrix_save[i][j]);
    }
    fprintf(save_beta_alpha, "\n");
  }
  fclose(save_beta_alpha);

  FILE* alpha_energy;
  alpha_energy = fopen("alpha_energies.dat", "w");
  for (i=0; i<different_alphas; i++) {
    fprintf(alpha_energy, "%f \t %f\n", alpha_values[i], alpha_energy_values[i]);
    }
  fclose(alpha_energy);

  // Calculate variance
  sum_tmp = 0;
  for(i = 0; i < nbr_iterations; i++){
    sum_tmp += energy[i] * energy[i];
  }
  variance = sum_tmp / nbr_iterations - I_value * I_value;
  error = sqrt(variance / nbr_iterations);


  // Checking statistical inefficiency
  int k= 20;
  int max_block = nbr_iterations/1000;
  double s_block[max_block];

  for (int B = 1; B<max_block+1; B++){
    s_block[B-1] = block_average(energy, nbr_iterations, B);
  }
  int s_acf = auto_corr_fun(energy, nbr_iterations, k);

  double s_block_avg_avg = 0;
  int start_avg = max_block/2;
  for (i=start_avg; i<max_block; i++){
    s_block_avg_avg += s_block[i] / (double)(max_block-start_avg);
  }

  printf("s_corr = %d\ns_block =%f\n", s_acf,s_block_avg_avg);


  // Deallocate rng
  gsl_rng_free (q);

  FILE *fp;

  //Write to file, energy_data_eq
  fp = fopen("energy_data_eq.dat","w");
  for (i = 0; i < nbr_iterations; i++){
    fprintf(fp, "%ld \t %e", i, energy[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);

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

double gradient_E_alpha(double * energy, double m1[], double m2[], double alpha, int length_E, double* ln_wavefun){
  double grad_E_alpha;

  double length_m12 = sqrt((m1[0]-m2[0])*(m1[0]-m2[0]) + (m1[1]-m2[1])*(m1[1]-m2[1]) + (m1[2]-m2[2])*(m1[2]-m2[2]));
  double length_m1 = sqrt((m1[0]*m1[0])+(m1[1]*m1[1])+(m1[2]*m1[2]));
  double length_m2 = sqrt((m2[0]*m2[0])+(m2[1]*m2[1])+(m2[2]*m2[2]));

  double cross_mult = m1[0]*m2[0]+m1[1]*m2[1]+m1[2]*m2[2];
  double m1_squared = ((m1[0]*m1[0])+(m1[1]*m1[1])+(m1[2]*m1[2]));
  double m2_squared = ((m2[0]*m2[0])+(m2[1]*m2[1])+(m2[2]*m2[2]));
  double a_r12 = (1.0+alpha*length_m12);


  double average_1 = 0;
  double average_2a = 0;
  double average_2b = 0;

  for (int i=0; i<length_E; i++){
    average_1 += energy[i]*ln_wavefun[i] / (double) length_E;
    average_2a += energy[i] / (double) length_E;
    average_2b += ln_wavefun[i] / (double) length_E;
  }

  grad_E_alpha = 2*(average_1-average_2a*average_2b);
  return grad_E_alpha;
}

void get_ln_wavefun(double * ln_wavefun, double m1[], double m2[], double alpha, int iteration){

  double length_m12 = sqrt((m1[0]-m2[0])*(m1[0]-m2[0]) + (m1[1]-m2[1])*(m1[1]-m2[1]) + (m1[2]-m2[2])*(m1[2]-m2[2]));

  ln_wavefun[iteration] = -length_m12*length_m12 / pow(1.0+alpha*length_m12, 2) / 2;

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

