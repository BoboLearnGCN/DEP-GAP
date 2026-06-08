#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
# include <chrono> // time
#include <unistd.h> // getopt
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h>  //mkdir
#include <algorithm>    // sort
#include <vector>  //

using namespace std;

double kernel_dnorm(double x, double mean, double sigma_sq) {
    return -pow(x - mean, 2) / 2 / sigma_sq;
}

double my_max(const double *_seq, int _length) {
    double temp_max = _seq[0];
    for (int i = 1; i < _length; ++i) {
        if (temp_max < _seq[i]) {
            temp_max = _seq[i];
        }
    }
    return temp_max;
}

int which_my_max(const double *_seq, int _length) {
    double temp_max = _seq[0];
    int index = 0;
    for (int i = 1; i < _length; ++i) {
        if (temp_max < _seq[i]) {
            temp_max = _seq[i];
            index = i;
        }
    }
    return index;
}

int rand_cate(const double *_prop, omprng _rng) {
    int res = 0;
    double u = _rng.runif();
    while (u > _prop[res]) {
        u = u - _prop[res];
        res++;
    }
    return res;
}

int rand_cate_uni(int _K, omprng _rng) {
    int res = 0;
    auto *prob = new double[_K];
    for (int k = 0; k < _K; ++k) {
        prob[k] = 1.0 / _K;
    }

    double u = _rng.runif();
    while (u > prob[res]) {
        u = u - prob[res];
        res++;
    }
    delete[] prob;
    return res;
}

void _update_mu(int _D, int *_D_vec, int _K,//dimension
                double *_m_mu, double _tau_mu_sq,//prior
                double *_Y_D, // observed
                int **_S, // observed + latent
                double *_gamma, double *_sigma_sq, //parameter
                omprng _rng,
                double *_mu) {

    int cell_index = 0, k;
    auto *post_mu = new double[_K](), *post_sigma_sq = new double[_K]();

    for (int d = 0; d < _D; d++) {
        for (int i = 0; i < _D_vec[d]; i++) {
            k = _S[d][i];
//            if (k > 0) {
//                //numerator
//                post_mu[k] += (_Y_D[cell_index] - _alpha - _gamma[d]) / _sigma_sq[d];
//                //denominator
//                post_sigma_sq[k] += 1 / _sigma_sq[d];
//            }

            //numerator
            post_mu[k] += (_Y_D[cell_index] - _gamma[d]) / _sigma_sq[d];
            //denominator
            post_sigma_sq[k] += 1 / _sigma_sq[d];
            cell_index++;
        }
    }

    // add prior information
    for (k = 0; k < _K; ++k) {
        post_sigma_sq[k] = 1 / (post_sigma_sq[k] + 1 / _tau_mu_sq);
        post_mu[k] = (post_mu[k] + _m_mu[k] / _tau_mu_sq) * post_sigma_sq[k];
        _mu[k] = _rng.rnorm(post_mu[k], sqrt(post_sigma_sq[k]));
    }

    delete[] post_mu;
    delete[] post_sigma_sq;
}


void _update_gamma(int _D, int *_D_vec, // dimension
                   double _tau_gamma_sq, // prior
                   double *_Y_D, //observed
                   int **_S, //observed + latent
                   double *_mu, double *_sigma_sq, //parameters
                   omprng _rng,
                   double *_gamma) {

    int cell_index = _D_vec[0];
    int k;
    auto *post_mu = new double[_D](),
            *post_sigma_sq = new double[_D]();

    for (int d = 1; d < _D; ++d) {
        for (int i = 0; i < _D_vec[d]; ++i) {
            k = _S[d][i];
            // numerator
            post_mu[d] += (_Y_D[cell_index] - _mu[k]) / _sigma_sq[d];

            cell_index++;
        }
        // denominator
        post_sigma_sq[d] += _D_vec[d] / _sigma_sq[d];
    }

    // add prior information
    for (int b = 1; b < _D; b++) {
        post_sigma_sq[b] = 1 / (post_sigma_sq[b] + 1 / _tau_gamma_sq);
        post_mu[b] = post_mu[b] * post_sigma_sq[b];
        // sampling
        _gamma[b] = _rng.rnorm(post_mu[b], sqrt(post_sigma_sq[b]));
    }

    delete[] post_mu;
    delete[] post_sigma_sq;
}


void _update_sigma_sq(int _D, int *_D_vec, // dimension
                      double _a, double _b, // prior
                      double *_Y_D, //observed
                      int **_S, // observed + latent
                      double *_mu, double *_gamma, // parameters
                      omprng _rng,
                      double *_sigma_sq) {

    int cell_index = 0, k;
    auto *_a_post = new double[_D](), *_b_post = new double[_D]();

    for (int d = 0; d < _D; ++d) {
        _a_post[d] = _D_vec[d] / 2.0 + _a;
        _b_post[d] = _b;
        for (int i = 0; i < _D_vec[d]; ++i) {
            k = _S[d][i];
            // para_2
            _b_post[d] += pow(_Y_D[cell_index] - _mu[k] - _gamma[d], 2) / 2;

            cell_index++;
        }
        _sigma_sq[d] = 1.0 / _rng.rgamma(_a_post[d], 1.0 / (_b_post[d]));
    }

    delete[] _a_post;
    delete[] _b_post;

}

//Gibbs
void _update_S(int _d, int _n_d, int _K, int _G, // dimension
               double **_Y_D, // observed
               double **_mu, double **_gamma, double **_sigma_sq, double **_psi, //parameters
               int **_neighbor_list,
               omprng _rng,
               int *_S) {

    int k_neigh;
    auto *prob = new double[_K]();
    int index = _d * _n_d;

    for (int i = 0; i < _n_d; ++i) {
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = 0.0;
            // feature information
            for (int g = 0; g < _G; ++g) {
                prob[k_t] += kernel_dnorm(_Y_D[g][index],
                                          _mu[g][k_t] + _gamma[g][_d], _sigma_sq[g][_d]);
            }

            // potts penalty
            for (int l = 0; l < 4; ++l) {
                if (_neighbor_list[i][l] >= 0) {
                    k_neigh = _S[_neighbor_list[i][l]];
                    prob[k_t] += _psi[k_t][k_neigh];
                }
            }
        }
        double prob_max = my_max(prob, _K);
        double prob_norm = 0;
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = exp(prob[k_t] - prob_max);
            prob_norm += prob[k_t];
        }

        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = prob[k_t] / prob_norm;
        }

        _S[i] = rand_cate(prob, _rng);
        index++;
    }

    delete[] prob;
}


// Gibbs
void _proposal_S(int _n_d, int _K, // dimension
                 double **_psi, //parameters
                 int **_neighbor_list,
                 omprng _rng,
                 int *_S) {
    int k_neigh;
    auto *prob = new double[_K]();

    for (int i = 0; i < _n_d; ++i) {
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = 0.0;
        }

        for (int k_t = 0; k_t < _K; ++k_t) {
            // potts penalty
            for (int l = 0; l < 4; ++l) {
                if (_neighbor_list[i][l] >= 0) {
                    k_neigh = _S[_neighbor_list[i][l]];
                    prob[k_t] += _psi[k_t][k_neigh];
                }
            }
        }

        double prob_max = my_max(prob, _K);
        double prob_norm = 0;
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = exp(prob[k_t] - prob_max);
            prob_norm += prob[k_t];
        }

        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = prob[k_t] / prob_norm;
        }

        _S[i] = rand_cate(prob, _rng);
    }

    delete[] prob;
}

void _update_psi(int _D, int *_n_vec_D, int _K,// dimension
                 int **_S, // latent
                 double _eta, double _tau_psi_sq, double _tau_0_sq, // prior
                 int ***_neighbor_list,
                 omprng _rng,
                 double **_psi) {
    auto **psi_candidate = new double *[_K];
    //initialize
    for (int k = 0; k < _K; ++k) {
        psi_candidate[k] = new double[_K]();
    }

    int **S_candidate = new int *[_D];
    for (int d = 0; d < _D; ++d) {
        S_candidate[d] = new int[_n_vec_D[d]]();
    }

    int k_candidate, k_t, nei_candidate, nei_t;
    double log_r;

    for (int i = 0; i < _K - 1; ++i) {
        for (int j = i + 1; j < _K; ++j) {

            for (int _i = 0; _i < _K - 1; ++_i) {
                for (int _j = _i + 1; _j < _K; ++_j) {
                    psi_candidate[_i][_j] = psi_candidate[_j][_i] = _psi[_i][_j];
                }
            }
            psi_candidate[i][j] = psi_candidate[j][i] = _rng.rnorm(_psi[i][j], sqrt(_tau_0_sq));

            for (int d = 0; d < _D; ++d) {
                for (int i_s = 0; i_s < _n_vec_D[d]; ++i_s) {
                    S_candidate[d][i_s] = _S[d][i_s];
                }
                _proposal_S(_n_vec_D[d], _K, _psi,
                            _neighbor_list[d], _rng, S_candidate[d]);
            }

            // prior
            log_r = kernel_dnorm(psi_candidate[i][j], _eta, _tau_psi_sq)
                    - kernel_dnorm(_psi[i][j], _eta, _tau_psi_sq);

            // potts
            for (int d = 0; d < _D; ++d) {
                for (int _i = 0; _i < _n_vec_D[d]; ++_i) {
                    for (int l = 0; l < 4; ++l) {
                        if (_neighbor_list[d][_i][l] >= 0) {
                            k_candidate = S_candidate[d][_i];
                            nei_candidate = S_candidate[d][_neighbor_list[d][_i][l]];
                            k_t = _S[d][_i];
                            nei_t = _S[d][_neighbor_list[d][_i][l]];

                            log_r += psi_candidate[k_t][nei_t] + _psi[k_candidate][nei_candidate] -
                                     psi_candidate[k_candidate][nei_candidate] - _psi[k_t][nei_t];
                        }
                    }
                }
            }

            if (log_r > log(_rng.runif())) {
                _psi[i][j] = _psi[j][i] = psi_candidate[i][j];
            }
        }
    }

    for (int k = 0; k < _K; ++k) {
        delete[] psi_candidate[k];
    }
    delete[] psi_candidate;

    for (int d = 0; d < _D; ++d) {
        delete[] S_candidate[d];
    }
    delete[] S_candidate;

}

int _location_S(int _dim_x, int _dim_y, int _x, int _y) {
    if (_x == -1) { return -1; }
    if (_x == _dim_x) { return -1; }
    if (_y == -1) { return -1; }
    if (_y == _dim_y) { return -1; }
    else { return ((_y) * _dim_x + _x); }
}

void _neighbor_S(int _dim_x, int _dim_y, int _x, int _y, int *_neighbor) {
    _neighbor[0] = _location_S(_dim_x, _dim_y, _x - 1, _y);
    _neighbor[1] = _location_S(_dim_x, _dim_y, _x, _y - 1);
    _neighbor[2] = _location_S(_dim_x, _dim_y, _x + 1, _y);
    _neighbor[3] = _location_S(_dim_x, _dim_y, _x, _y + 1);
}


int main(int argc, char **argv) {
    ////////////////////////
    // 1. Load count data //
    ////////////////////////
    // allocate memory for the observation

    int G, K, B, N = 0, dim_X, dim_Y, core_number;
    int n_iter, n_record;
    int seed;

    string output_list, data_name, label_name;

    int opt;
    while ((opt = getopt(argc, argv, ":g:k:b:n:x:y:c::d:l:o::t::s::")) != -1) {
        switch (opt) {
            case 'g':
                G = atoi(optarg);
                cout << "The number of features is " << G << "." << endl;
                break;
            case 'k':
                K = atoi(optarg);
                cout << "The number of labels is " << K << "." << endl;
                break;
            case 'b':
                B = atoi(optarg);
                cout << "The number of batches is " << B << "." << endl;
                break;
            case 'x':
                dim_X = atoi(optarg);
                cout << "The dimension 1 is " << dim_X << "." << endl;
                break;
            case 'y':
                dim_Y = atoi(optarg);
                cout << "The dimension 2 is " << dim_Y << "." << endl;
                break;
            case 'c':
                if (optarg) {
                    core_number = atoi(optarg);
                } else {
                    core_number = 1;
                }
                cout << "The core number is " << core_number << "." << endl;
                omp_set_num_threads(core_number);
                cout << "Actually, there are total " << omp_get_max_threads() << "." << endl;
                break;
            case 'd':
                data_name.assign(optarg);
                cout << "Data file name is " << data_name.c_str() << endl;
                break;
            case 'l':
                label_name.assign(optarg);
                cout << "Label file name is " << label_name.c_str() << endl;
                break;
            case 'o':
                if (optarg) {
                    output_list.assign(optarg);
                } else {
                    output_list.assign("./");
                }
//                output_list.assign(optarg);
                cout << "Output directory is " << output_list.c_str() << endl;
                break;
            case 't':
                if (optarg) {
                    n_iter = atoi(optarg);
                } else {
                    n_iter = 3000;
                }
//                n_iter = atoi(optarg);
                cout << "The iteration time is " << n_iter << "." << endl;
                n_record = n_iter / 2;
                cout << "The burn-in time is " << n_record << "." << endl;
                break;
            case 's':
                if (optarg) {
                    seed = atoi(optarg);
                } else {
                    seed = 666;
                }
//                seed = atoi(optarg);
                cout << "set seed as " << seed << "." << endl;
                break;
            case '?':
                cout << "Unknown option " << optarg << endl;
                break;
            case ':':
                cout << "Missing option for " << optarg << endl;
                break;
            default: /* '?' */
                cout << "Error Usage: " << argv[0]
                     << " [-g GeneNumber] [-k CellTypeNumber] [-b BatchNumber] [-x DimX] [-y DimY] [-c CoreNumber] [-d DataName] [-o OutputDir] [-t Iteration] [-s Seed]"
                     << endl;
                exit(EXIT_FAILURE);

        }
    }

    int *n_vec = new int[B];

    for (int b = 0; b < B; ++b) {
        n_vec[b] = dim_X * dim_Y;
        N += n_vec[b];
        cout << "The sample size of " << b + 1 << " batch is " << n_vec[b] << "." << endl;
    }
    cout << "The total sample size is " << N << endl;
    cout << endl;

    cout << "The iteration is " << n_iter << "." << endl;
    cout << "The burn-in is " << n_record << "." << endl;
    cout << "set seed as " << seed << "." << endl;
    cout << endl;

    // loading labels
//    ifstream label_data;
//    label_data.open("Y_labels.txt");
//    if (!label_data) {
//        cout << "Unable to open the file: Y_labels.txt" << endl;
//        exit(1); // terminate with error
//    }

    int **S_t = new int *[B];
    for (int b = 0; b < B; ++b) {
        S_t[b] = new int[n_vec[b]]();
//        for (int j = 0; j < n_vec[b]; ++j) {
//            label_data >> S_t[b][j];
//        }
    }
//    label_data.close();
//    cout << "Finish label assignment." << endl;




    // loading observed data (feature)
    ifstream obs_data;
    obs_data.open(data_name);
    if (!obs_data) {
        cout << "Unable to open the data!" << endl;
        exit(1); // terminate with error
    }
    auto **Y = new double *[G];
    for (int g = 0; g < G; ++g) {
        Y[g] = new double[N]();
        for (int i = 0; i < N; ++i) {
            obs_data >> Y[g][i];
        }
    }

    obs_data.close();
    cout << "Load all feature data successfully!" << endl;

    // set-up spatial information
    int ***neighbor_list = new int **[B];
    int index_S;
    for (int d = 0; d < B; ++d) {
        neighbor_list[d] = new int *[n_vec[d]];
        index_S = 0;
        for (int y = 0; y < dim_Y; ++y) {
            for (int x = 0; x < dim_X; ++x) {
                neighbor_list[d][index_S] = new int[4]();
                _neighbor_S(dim_X, dim_Y, x, y, neighbor_list[d][index_S]);

                index_S++;
            }
        }
    }


    ///////////////////////////////////////////////////////
    // 2. Set hyper-parameters and Initialize parameters //
    ///////////////////////////////////////////////////////
    cout << "Set initial values." << endl;
    omprng MCMC_Rng;
    cout << "Initializing rng." << endl;
    MCMC_Rng.fixedSeed(seed);
    MCMC_Rng.setNumThreads(core_number);

    //hyper-parameters
    double tau_mu_sq = 5, tau_gamma_sq = 5,
            a_s = 1, b_s = 0.01, eta = 0, tau_psi_sq = 0.01, tau_0_sq = 0.1;
    auto **m_mu = new double *[G];

    ///////////////////////////////////////////////////////
    // 3. Data Generation //
    ///////////////////////////////////////////////////////
    ofstream output_temp;

    auto start_MCMC = chrono::system_clock::now();
    auto end_MCMC = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds_MCMC = end_MCMC - start_MCMC;

    string out_file;


    // initialization: data-dependent initial values

    auto **mu_t = new double *[G]();
    auto **gamma_t = new double *[G]();
    auto **sigma_sq_t = new double *[G]();
    auto **psi_t = new double *[K]();

    for (int g = 0; g < G; ++g) {
        mu_t[g] = new double[K]();
        m_mu[g] = new double[K]();

        gamma_t[g] = new double[B]();
        sigma_sq_t[g] = new double[B]();
    }

    for (int k = 0; k < K; ++k) {
        psi_t[k] = new double[K]();
    }

    obs_data.open(label_name);
    if (!obs_data) {
        cout << "Unable to open the label!" << endl;
        exit(1); // terminate with error
    }

    for (int b = 0; b < B; ++b) {
        for (int i = 0; i < n_vec[b]; ++i) {
//            S_t[b][i] = rand_cate_uni(K, MCMC_Rng);
            obs_data >> S_t[b][i];
        }
    }
    obs_data.close();
    cout << "Set initial labels." << endl;

    auto ***raw_sums = new double **[G];
    int ***count_k = new int **[G];
    int cell_k;
    int index;
    for (int g = 0; g < G; ++g) {
        raw_sums[g] = new double *[B];
        count_k[g] = new int *[B];
        index = 0;
        for (int b = 0; b < B; ++b) {
            raw_sums[g][b] = new double[K]();
            count_k[g][b] = new int[K]();

            for (int j = 0; j < n_vec[b]; ++j) {
                cell_k = S_t[b][j];
                raw_sums[g][b][cell_k] += Y[g][index];
                count_k[g][b][cell_k]++;

                index++;
            }
        }
    }

    for (int g = 0; g < G; ++g) {
        for (int k = 0; k < K; ++k) {
            mu_t[g][k] = m_mu[g][k] = raw_sums[g][0][k] / (count_k[g][0][k] + 1);
        }

        for (int b = 1; b < B; ++b) {
            gamma_t[g][b] = (raw_sums[g][b][0]) / n_vec[b] - mu_t[g][0];
        }
    }


#pragma omp parallel for
    for (int g = 0; g < G; ++g) {
        _update_sigma_sq(B, n_vec,
                         a_s, b_s,
                         Y[g], S_t,
                         mu_t[g], gamma_t[g],
                         MCMC_Rng,
                         sigma_sq_t[g]);
    }

    for (int i = 0; i < K - 1; ++i) {
        for (int j = i + 1; j < K; ++j) {
            psi_t[i][j] = psi_t[j][i] = -1.0;
        }
    }

    //////////////////////
    // 4. MCMC sampling //
    //////////////////////
    // allocate memory to store MCMC samplings
    int n_output = n_iter / 10;
    double ***mu_record = new double **[n_output];
    double ***gamma_record = new double **[n_output];
    double ***sigma_sq_record = new double **[n_output];
    int ***S_record = new int **[n_output];
    double ***psi_record = new double **[n_output];

    for (int i = 0; i < n_output; ++i) {
        mu_record[i] = new double *[G];
        gamma_record[i] = new double *[G];
        sigma_sq_record[i] = new double *[G];
        for (int g = 0; g < G; ++g) {
            mu_record[i][g] = new double[K]();
            gamma_record[i][g] = new double[B]();
            sigma_sq_record[i][g] = new double[B]();
        }

        S_record[i] = new int *[B];
        for (int b = 0; b < B; ++b) {
            S_record[i][b] = new int[n_vec[b]]();
        }

        psi_record[i] = new double *[K];
        for (int k = 0; k < K; ++k) {
            psi_record[i][k] = new double[K]();
        }
    }

    start_MCMC = chrono::system_clock::now();
    cout << "Start MCMC sampling." << endl;


    int t_iter = 0;
    for (int out = 0; out < 10; ++out) {


        for (int iter = 0; iter < n_output; ++iter) {
#pragma omp parallel for
            for (int g = 0; g < G; ++g) {

                //update mu
                _update_mu(B, n_vec, K,
                           m_mu[g], tau_mu_sq,
                           Y[g], S_t,
                           gamma_t[g], sigma_sq_t[g],
                           MCMC_Rng,
                           mu_t[g]);

                //update gamma
                _update_gamma(B, n_vec,
                              tau_gamma_sq,
                              Y[g], S_t,
                              mu_t[g], sigma_sq_t[g],
                              MCMC_Rng,
                              gamma_t[g]);

                //update sigma_sq
                _update_sigma_sq(B, n_vec,
                                 a_s, b_s,
                                 Y[g], S_t,
                                 mu_t[g], gamma_t[g],
                                 MCMC_Rng,
                                 sigma_sq_t[g]);

            }

            // update S
#pragma omp parallel for
            for (int b = 0; b < B; ++b) {

                _update_S(b, n_vec[b], K, G,
                          Y,
                          mu_t, gamma_t, sigma_sq_t, psi_t,
                          neighbor_list[b],
                          MCMC_Rng,
                          S_t[b]);

/*            cout << "update S!" << endl;*/

            }

            //record samplings
            for (int g = 0; g < G; ++g) {
                for (int k = 0; k < K; ++k) {
                    mu_record[iter][g][k] = mu_t[g][k];
                }

                for (int b = 0; b < B; ++b) {
                    gamma_record[iter][g][b] = gamma_t[g][b];
                    sigma_sq_record[iter][g][b] = sigma_sq_t[g][b];
                }
            }


            for (int b = 0; b < B; ++b) {
                for (int i = 0; i < n_vec[b]; ++i) {
                    S_record[iter][b][i] = S_t[b][i];
                }
            }

            // update psi
            _update_psi(B, n_vec, K,
                        S_t,
                        eta, tau_psi_sq, tau_0_sq,
                        neighbor_list,
                        MCMC_Rng,
                        psi_t);

            //record samplings of psi
            for (int j = 0; j < K - 1; ++j) {
                for (int k = j + 1; k < K; ++k) {
                    psi_record[iter][j][k] = psi_record[iter][k][j] = psi_t[j][k];
                }
            }

            t_iter++;
            if (t_iter % 100 == 0) {
                cout << "Iteration " << t_iter << endl;
            }

        }

        //////////////////////////////////////
        // 5. output the posterior sampling //
        //////////////////////////////////////


        cout << "output samplings now" << endl;

        out_file = output_list + "/mu_" + to_string(out) + ".txt";
        output_temp.open(out_file.c_str(), ios::out | ios::app);
        for (int g = 0; g < G; ++g) {
            for (int k = 0; k < K; ++k) {
                for (int iter = 0; iter < n_output; ++iter) {
                    output_temp << mu_record[iter][g][k];
                    output_temp << " ";
                }
                output_temp << endl;
            }
        }
        output_temp.close();
        cout << "mu finished!" << endl;

        out_file = output_list + "/gamma_" + to_string(out) + ".txt";
        output_temp.open(out_file.c_str(), ios::out | ios::app);
        for (int g = 0; g < G; ++g) {
            for (int b = 0; b < B; ++b) {
                for (int iter = 0; iter < n_output; ++iter) {
                    output_temp << gamma_record[iter][g][b];
                    output_temp << " ";
                }
                output_temp << endl;
            }
        }
        output_temp.close();
        cout << "gamma finished!" << endl;

        out_file = output_list + "/sigma_square_" + to_string(out) + ".txt";
        output_temp.open(out_file.c_str(), ios::out | ios::app);
        for (int g = 0; g < G; ++g) {
            for (int b = 0; b < B; ++b) {
                for (int iter = 0; iter < n_output; ++iter) {
                    output_temp << sigma_sq_record[iter][g][b];
                    output_temp << " ";
                }
                output_temp << endl;
            }
        }
        output_temp.close();
        cout << "sigma square finished!" << endl;


        out_file = output_list + "/S_" + to_string(out) + ".txt";
        output_temp.open(out_file.c_str(), ios::out | ios::app);
        auto *temp_count = new double[K]();

        for (int b = 0; b < B; ++b) {
            for (int j = 0; j < n_vec[b]; ++j) {
                for (int iter = 0; iter < n_output; ++iter) {
                    output_temp << S_record[iter][b][j];
                    output_temp << " ";
                }
                output_temp << endl;
            }
        }
        output_temp.close();
        cout << "S finished!" << endl;

        out_file = output_list + "/psi_" + to_string(out) + ".txt";
        output_temp.open(out_file.c_str(), ios::out | ios::app);
        for (int k = 0; k < K - 1; ++k) {
            for (int j = k + 1; j < K; ++j) {
                for (int iter = 0; iter < n_output; ++iter) {
                    output_temp << psi_record[iter][k][j];
                    output_temp << " ";
                }
                output_temp << endl;
            }
        }
        output_temp.close();
        cout << "psi finished!" << endl;


//        //free the memory
//
//        for (int iter = 0; iter < n_output; ++iter) {
//            for (int g = 0; g < G; ++g) {
//                delete[] mu_record[iter][g];
//                delete[] gamma_record[iter][g];
//                delete[] sigma_sq_record[iter][g];
//            }
//            delete[] mu_record[iter];
//            delete[] gamma_record[iter];
//            delete[] sigma_sq_record[iter];
//
//            for (int b = 0; b < B; ++b) {
//                delete[] S_record[iter][b];
//            }
//            delete[] S_record[iter];
//
//            for (int k = 0; k < K; ++k) {
//                delete[] psi_record[k];
//            }
//            delete[] psi_record[iter];
//        }
//        delete[] mu_record;
//        delete[] gamma_record;
//        delete[] sigma_sq_record;
//        delete[] S_record;
//        delete[] psi_record;

    }


    end_MCMC = chrono::system_clock::now();
    elapsed_seconds_MCMC = end_MCMC - start_MCMC;
    cout << "Time of " << n_iter << " iterations of MCMC sampling is: "
         << elapsed_seconds_MCMC.count() / 60.0 << "min" << endl;


    for (int g = 0; g < G; ++g) {
        delete[] mu_t[g];
        delete[] gamma_t[g];
        delete[] sigma_sq_t[g];
    }
    delete[] mu_t;
    delete[] gamma_t;
    delete[] sigma_sq_t;

    for (int b = 0; b < B; ++b) {
        delete[] S_t[b];
    }
    delete[] S_t;
    for (int k = 0; k < K; ++k) {
        delete[] psi_t[k];
    }
    delete[] psi_t;

    cout << "MCMC finished!" << endl;
    cout << "" << endl;


    return 0;
}

