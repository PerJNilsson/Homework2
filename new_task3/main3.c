#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.141592653589
#define nbr_iterations 500000
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
  double sum_tmp, sum_sq_tmp;
  double rand_nbr;
  double m1[nbr_dim] = {0};
  double m2[nbr_dim] = {0};
  double n1[nbr_dim] = {0};
  double n2[nbr_dim] = {0};
  double p_m, p_n;
  double I_value;
  double variance;
  double error;
  double tmp_var;
  double delta = 0.975;
  double alpha = 0.1;
  int equi_phase = 1500;
  size_t different_alphas = 40;
  size_t nbr_switching_state = 0;
  size_t iterations_per_alpha = 1000;
  double s = 10.5;

  const gsl_rng_type *T; /* static info about rngs */
  gsl_rng *q; /* rng instance */
  gsl_rng_env_setup(); /* setup the rngs */
  T = gsl_rng_default; /* specify default rng */
  q = gsl_rng_alloc(T); /* allocate default rng */
  gsl_rng_set(q, time(NULL)); /* Initialize rng */


  double * energy = malloc(sizeof(double) * nbr_iterations);
  double * energy_squared = malloc(sizeof(double) * nbr_iterations);
  double * alpha_energy_values = malloc(sizeof(double)*different_alphas);
  double * alpha_energy_std = malloc(sizeof(double)*different_alphas);
  double * alpha_values = malloc(sizeof(double)*different_alphas);
  double * statistical_s = malloc(sizeof(double)*different_alphas);
  double * same_alpha_energy = malloc(sizeof(double)*iterations_per_alpha);
  double * same_alpha_sigma = malloc(sizeof(double)*different_alphas);
  
  

  for (size_t jx=0; jx<different_alphas; jx++){
    alpha = 0.05 + jx*0.20 / (different_alphas-1.0);
    alpha_values[jx] = alpha;
    printf("alpha=%f\n", alpha);
    for (size_t kx =0; kx<iterations_per_alpha; kx++){
      nbr_switching_state = 0;
      
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
          nbr_switching_state++;
        } else {
          energy[i] = local_energy(m1, m2, alpha);
        }
       
      }// End main loop
      sum_sq_tmp = 0;
      sum_tmp = 0;
      for(i = 0; i < nbr_iterations; i++){
        sum_tmp += energy[i];
        sum_sq_tmp += energy[i]*energy[i] /nbr_iterations;
      }
      I_value = sum_tmp / nbr_iterations;

      // Using the statistical inefficiency to calculate the variance
      tmp_var = (s/nbr_iterations)* (sum_sq_tmp - I_value*I_value);
      
      same_alpha_energy[kx] = I_value;
      alpha_energy_values[jx] += I_value / iterations_per_alpha;
      alpha_energy_std[jx] += tmp_var / (iterations_per_alpha*iterations_per_alpha);//I_value*I_value / iterations_per_alpha;
    } // end iterations per alpha, kx

    // Using the average for each alpha to calculate the variance
    sum_sq_tmp = 0;
    sum_tmp = 0;
    for (int i=0; i<iterations_per_alpha; i++){
            sum_tmp += same_alpha_energy[i] / iterations_per_alpha;
            sum_sq_tmp += same_alpha_energy[i]*same_alpha_energy[i] / iterations_per_alpha;
    }
    same_alpha_sigma[jx] = sqrt(sum_sq_tmp - sum_tmp*sum_tmp)/sqrt(iterations_per_alpha);
  } // end different alphas, jx
  
  FILE* alpha_energy;
  alpha_energy = fopen("alpha_energies.dat", "w");
  for (i=0; i<different_alphas; i++) {
          alpha_energy_std[i] = sqrt(alpha_energy_std[i]);//sqrt((alpha_energy_std[i] - alpha_energy_values[i]*alpha_energy_values[i]));
    fprintf(alpha_energy, "%f \t %f \t %f\n", alpha_values[i], alpha_energy_values[i], alpha_energy_std[i]);
    }
  fclose(alpha_energy);

  FILE* alpha_energy2;
  alpha_energy2 = fopen("alpha_energies2.dat", "w");
  for (i=0; i<different_alphas; i++) {
          fprintf(alpha_energy2, "%f \t %f \t %f\n", alpha_values[i], alpha_energy_values[i], same_alpha_sigma[i]);
    }
  fclose(alpha_energy2);

  // Deallocate rng
  gsl_rng_free (q);
}//End MAIN

// calculated the wavefunction given from problem paper
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

// Calculate the probability of finding e
double calculate_probability(double r1[nbr_dim], double r2[nbr_dim], double alpha){
  double value;
  value = calculate_wave_function(r1, r2, alpha);
  return value*value;
}

// Calculate the local energy given from problem
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

// AFC
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

