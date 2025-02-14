using namespace std;
#include <cstdio>
#include <vector>
// the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>

// for gamma function
#include <cmath>

vector<vector<int>> state_count(const vector<vector<uint8_t>> &data,
                                vector<int> &children, vector<int> &parents,
                                vector<int> &n_states, int &parallel) {
  // if parallel == 0 then use single thread, if parallel == 1 then use multi
  // thread(CPU), if parallel == 2 then use GPU (unimplemented)
  if (children.size() == 1) {
    int node_x = children.at(0);
    if (parents.empty()) {
      int x = n_states.at(node_x);
      vector<vector<int>> counts(x, vector<int>(1, 0));
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          counts.at(data.at(i).at(node_x)).at(0) += 1;
        }
      } else if (parallel == 1) {
#pragma omp parallel
        {
          vector<int> temp(x, 0);
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            temp.at(data.at(i).at(node_x)) += 1;
          }
          for (int j = 0; j < x; j++) {
#pragma omp atomic
            counts.at(j).at(0) += temp.at(j);
          }
        }
      }
      return counts;
    } else {
      // return the state counts of X, Y | Z shape: countmap[state of
      // child][state of parents]
      int x = n_states.at(node_x);
      int y = 1;
      for (int i = 0; i < (int)parents.size(); i++) {
        y = y * n_states.at(parents.at(i));
      }
      vector<vector<int>> counts(x, vector<int>(y, 0));
      // count the number of each state
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          int yy = 0;
          for (int j = 0; j < (int)parents.size(); j++) {
            if (j == 0) {
              yy = data.at(i).at(parents.at(j));
            } else {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
          }
          counts.at(data.at(i).at(node_x)).at(yy) += 1;
        }
      } else if (parallel == 1) {
        int yy;
#pragma omp parallel private(yy)
        {
          vector<vector<int>> temp(x, vector<int>(y, 0));
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            yy = data.at(i).at(parents.at(0));
            for (int j = 1; j < (int)parents.size(); j++) {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
            temp.at(data.at(i).at(node_x)).at(yy) += 1;
          }
          for (int j = 0; j < x; j++) {
            for (int k = 0; k < y; k++) {
#pragma omp atomic
              counts.at(j).at(k) += temp.at(j).at(k);
            }
          }
        }
      }
      return counts;
    }
  } else {  // if children.size() == 2
    if (parents.empty()) {
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      vector<vector<int>> counts(len, vector<int>(1, 0));
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          counts
              .at(data.at(i).at(children.at(0)) * yc +
                  data.at(i).at(children.at(1)))
              .at(0) += 1;
        }
      } else if (parallel == 1) {
#pragma omp parallel
        {
          vector<int> temp(len, 0);
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            temp.at(data.at(i).at(children.at(0)) * yc +
                    data.at(i).at(children.at(1))) += 1;
          }
          for (int j = 0; j < len; j++) {
#pragma omp atomic
            counts.at(j).at(0) += temp.at(j);
          }
        }
      }

      return counts;
    } else {
      // return the state counts of X, Y | Z shape: countmap[state of
      // child][state of parents]
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      int y = 1;
      for (int i = 0; i < (int)parents.size(); i++) {
        y = y * n_states.at(parents.at(i));
      }
      vector<vector<int>> counts(len, vector<int>(y, 0));
      // count the number of each state
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          int yy = 0;
          for (int j = 0; j < (int)parents.size(); j++) {
            if (j == 0) {
              yy = data.at(i).at(parents.at(j));
            } else {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
          }
          counts
              .at(data.at(i).at(children.at(0)) * yc +
                  data.at(i).at(children.at(1)))
              .at(yy) += 1;
        }
      } else if (parallel == 1) {
        int yy;
#pragma omp parallel private(yy)
        {
          vector<vector<int>> temp(len, vector<int>(y, 0));
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            yy = data.at(i).at(parents.at(0));
            for (int j = 1; j < (int)parents.size(); j++) {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
            temp.at(data.at(i).at(children.at(0)) * yc +
                    data.at(i).at(children.at(1)))
                .at(yy) += 1;
          }
          for (int j = 0; j < len; j++) {
            for (int k = 0; k < y; k++) {
#pragma omp atomic
              counts.at(j).at(k) += temp.at(j).at(k);
            }
          }
        }
      }

      return counts;
    }
  }
}

float natori_independent_score(const vector<vector<uint8_t>> &data, int &node_x,
                               int &node_y, vector<int> &parents,
                               vector<int> &n_states, int &parallel) {
  // return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5;  // hyperparameter
  vector<int> node_x_vec(1, node_x);
  vector<int> node_y_vec(1, node_y);
  if (parents.empty()) {
    // no parents
    vector<vector<int>> count;
    count = state_count(data, node_x_vec, parents, n_states, parallel);
    int r = count.size();  // number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
    vector<vector<int>> count2;
    count2 = state_count(data, node_y_vec, parents, n_states, parallel);
    r = count2.size();  // number of states of node_x
    n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count2.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count2.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  } else {
    // have parents
    vector<vector<int>> count;
    count = state_count(data, node_x_vec, parents, n_states, parallel);
    int q = count.at(0).size();  // number of states of parents
    int r = count.size();        // number of states of node_x
    vector<float> n_ij(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
    vector<vector<int>> count2;
    count2 = state_count(data, node_y_vec, parents, n_states, parallel);
    q = count2.at(0).size();  // number of states of parents
    r = count2.size();        // number of states of node_x
    vector<float> n_ij2(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij2.at(j) += count2.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score += lgamma(alpha + count2.at(k).at(j)) - lgamma(alpha);
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij2.at(j));
    }
  }
  // calculate the score
  return score;
}

float natori_dependent_score(const vector<vector<uint8_t>> &data, int &node_x,
                             int &node_y, vector<int> &parents,
                             vector<int> &n_states, int &parallel) {
  // return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5;

  vector<int> children{node_x, node_y};
  if (parents.empty()) {
    // no parents
    vector<vector<int>> count;
    count = state_count(data, children, parents, n_states, parallel);
    int r = count.size();  // number of states of node_x]
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  } else {
    // have parents
    vector<vector<int>> count;
    count = state_count(data, children, parents, n_states, parallel);
    int q = count.at(0).size();  // number of states of parents
    int r = count.size();        // number of states of node_x
    vector<float> n_ij(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
  }
  return score;
}

bool ci_test_by_Bayesfactor(const vector<vector<uint8_t>> &data, int &node_x,
                            int &node_y, vector<int> &Z, vector<int> &n_states, int &parallel) {
  // CI test for X _|_ Y | Z
  float independent_score = 0.0;
  float dependent_score = 0.0;
  independent_score += natori_independent_score(data, node_x, node_y, Z,
                                                n_states, parallel);
  dependent_score +=
      natori_dependent_score(data, node_x, node_y, Z, n_states, parallel);

  if (independent_score > dependent_score) {
    // cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<"
    // |"<<independent_score<<">"<<dependent_score<< endl;
    return true;
  } else {
    // cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<"
    // |"<<independent_score<<"<"<<dependent_score<< endl;
    return false;
  }
}

// double qnorm(double u) {
//   static double a[9] = {1.24818987e-4,   -1.075204047e-3, 5.198775019e-3,
//                         -0.019198292004, 0.059054035642,  -0.151968751364,
//                         0.319152932694,  -0.5319230073,   0.797884560593};
//   static double b[15] = {-4.5255659e-5,   1.5252929e-4,    -1.9538132e-5,
//                          -6.76904986e-4,  1.390604284e-3,  -7.9462082e-4,
//                          -2.034254874e-3, 6.549791214e-3,  -0.010557625006,
//                          0.011630447319,  -9.279453341e-3, 5.353579108e-3,
//                          -2.141268741e-3, 5.35310549e-4,   0.999936657524};
//   double w, y, z;
//   int i;

//   if (u == 0.) return 0.5;
//   y = u / 2.;
//   if (y < -6.) return 0.;
//   if (y > 6.) return 1.;
//   if (y < 0.) y = -y;
//   if (y < 1.) {
//     w = y * y;
//     z = a[0];
//     for (i = 1; i < 9; i++) z = z * w + a[i];
//     z *= (y * 2.);
//   } else {
//     y -= 2.;
//     z = b[0];
//     for (i = 1; i < 15; i++) z = z * y + b[i];
//   }

//   if (u < 0.) return (1. - z) / 2.;
//   return (1. + z) / 2.;
// }

// double qchi(double x2, int n) {
//   double w, pw, x, qc;
//   int i, i1;

//   if (n < 1) {
//     fprintf(stderr, "Error : 自由度 < 1 in qchi()!\n");
//     return 0.;
//   }
//   if (x2 <= 0.) return 1.;
//   if (x2 > 400.) return 0.;
//   if (n > 10) {
//     w = 2. / (9. * (double)n);
//     pw = pow(x2 / (double)n, 1. / 3.);
//     return 1. - qnorm((pw - 1. + w) / sqrt(w));
//   }
//   w = exp(-x2 / 2.);
//   if (n == 2) return w;
//   x = sqrt(x2);
//   if (n == 1) return 2. * (1. - qnorm(x));
//   if ((n % 2) == 0) {
//     i1 = 2;
//     qc = w;
//   } else {
//     i1 = 1;
//     qc = 2. * (1. - qnorm(x));
//     w *= (0.797884560750774 / x);
//   }
//   for (i = i1; i <= n - 2; i += 2) {
//     w *= (x2 / (double)i);
//     qc += w;
//   }
//   return qc;
// }

// bool ci_test_by_Chi_squared_test(const
//   vector<int> &children, vector<int> &parents,
//   vector<int> &n_states, int &parallel = 1) {
//   //CI test for X _|_ Y | Z using Chi-squared test; Test statistics G = 2 *
//   sum_x,y,z(n_xyz * log(n_xyz / (n_xz n_yz / n))) ~ chi^2_{d.f= (|X| - 1) *
//   (|Y| - 1) * |Z|} (from
//   http://www.ai.lab.uec.ac.jp/wp-content/uploads/2019/04/41f06ccbd0ac30d15c7728117770b105.pdf)
//   float threthold = 0.05;
//   vector<vector<int>>count;
//   count = state_count(data, children, parents, n_states, parallel);
//   int r_x = n_states.at(children.at(0));
//   int r_y = n_states.at(children.at(1));
//   int r_z = n_states.at(node_z);
//   int n = data.size();
//   vector<vector<double>> n_xz(r_x, vector<double>(r_z, 0));
//   vector<vector<double>> n_yz(r_y, vector<double>(r_z, 0));
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         n_xz.at(i).at(k) += count.at(i * r_y + j).at(k);
//         n_yz.at(j).at(k) += count.at(i * r_y + j).at(k);
//         //cout << count.at(i * r_y + j).at(k)<< ", ";
//       }
//       //cout << endl;
//     }
//   }
//   vector<double> sum_x(r_z, 0);
//   vector<double> sum_y(r_z, 0);
//   vector<double> sum_z(r_z, 0);
//   for (int i = 0; i < r_x; i++) {
//     for (int k = 0; k < r_z; k++) {
//       sum_x.at(k) += n_xz.at(i).at(k);
//     }
//   }
//   for (int j = 0; j < r_y; j++) {
//     for (int k = 0; k < r_z; k++) {
//       sum_y.at(k) += n_yz.at(j).at(k);
//     }
//   }
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         sum_z.at(k) += count.at(i * r_y + j).at(k);
//       }
//     }
//   }
//   //cout<< "n_xz"<<endl;
//   for (int i = 0; i < r_x; i++) {
//     for (int k = 0; k < r_z; k++) {
//       //cout << n_xz.at(i).at(k)/ sum_x.at(k)<< ", ";
//     }
//     //cout << endl;
//   }
//   //cout<< "n_yz"<<endl;
//   for (int j = 0; j < r_y; j++) {
//     for (int k = 0; k < r_z; k++) {
//       //cout << n_yz.at(j).at(k)/ sum_y.at(k)<< ", ";
//     }
//     //cout << endl;
//   }
//   double G = 0.0;
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         //cout << "ijk:" << (double)count.at(i * r_y + j).at(k)
//         /(double)n << "xzyz:" <<((double)n_xz.at(i).at(k) *
//         (double)n_yz.at(j).at(k) / (sum_x.at(k) * sum_y.at(k)))
//         *sum_z.at(k)/(double)n << "ijk/xzyz:"<< ((double)count.at(i * r_y +
//         j).at(k) / (double)n )/ (((double)n_xz.at(i).at(k) *
//         (double)n_yz.at(j).at(k) / (sum_x.at(k) *
//         sum_y.at(k)))*sum_z.at(k)/(double)n)<< endl;
//         // G = G + 2 * (double)count.at(i * r_y + j).at(k) *
//         log(((double)count.at(i * r_y + j).at(k) / (double)n )/
//         (((double)n_xz.at(i).at(k) * (double)n_yz.at(j).at(k) / (sum_x.at(k)
//         * sum_y.at(k)))*sum_z.at(k)/(double)n)); G = G + 2 *
//         (double)count.at(i * r_y + j).at(k) * log((double)count.at(i * r_y +
//         j).at(k) * (double)n / ((double)n_xz.at(i).at(k) *
//         (double)n_yz.at(j).at(k)));
//       }
//     }
//   }
//   //cout<< "G:" <<G<< endl;
//   int dim = (r_x - 1) * (r_y - 1) * r_z;
//   double p_value = qchi(G, dim);
//   //cout<< "p_value:" <<p_value<< endl;
//   bool flag = false;
//   if (p_value <= threthold || 1 - p_value <= threthold){
//     flag = true;
//   }
//   if (flag){
//     return true;//dependent, reject null hypothesis, find V-structure
//     //cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<< node_z <<
//     endl;
//   }
//   else{
//     //cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<< node_z
//     << endl; return false;
//   }
// }

// bool ci_test_by_CMI(const vector<vector<uint8_t>>
// &data, int &node_x, int &node_y, int &node_z, vector<int> &n_states) {
//   //CI test for X _|_ Y | Z using CMI; Test statistics cmi(x,y|z) = sum_x,y,z
//   p(x,y,z) * log(p(x,y|z)/p(x|z)p(y|z)) double threthold = 0.003; vector<int>
//   children(2); children.at(0) = node_x; children.at(1) = node_y; vector<int>
//   parents(1, node_z); vector<vector<int>>count; count = state_count(data,
//   children, parents, n_states, parallel); int r_x = n_states.at(node_x); int
//   r_y = n_states.at(node_y); int r_z = n_states.at(node_z); int n =
//   data.size(); vector<vector<double>> n_xz(r_x, vector<double>(r_z, 0));
//   vector<vector<double>> n_yz(r_y, vector<double>(r_z, 0));
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         n_xz.at(i).at(k) += count.at(i * r_y + j).at(k);
//         n_yz.at(j).at(k) += count.at(i * r_y + j).at(k);
//         //cout << count.at(i * r_y + j).at(k)<< ", ";
//       }
//       //cout << endl;
//     }
//   }
//   vector<double> sum_x(r_z, 0);
//   vector<double> sum_y(r_z, 0);
//   vector<double> sum_z(r_z, 0);
//   for (int i = 0; i < r_x; i++) {
//     for (int k = 0; k < r_z; k++) {
//       sum_x.at(k) += n_xz.at(i).at(k);
//     }
//   }
//   for (int j = 0; j < r_y; j++) {
//     for (int k = 0; k < r_z; k++) {
//       sum_y.at(k) += n_yz.at(j).at(k);
//     }
//   }
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         sum_z.at(k) += count.at(i * r_y + j).at(k);
//       }
//     }
//   }

//   vector<vector<double>> p_xyz(r_x * r_y, vector<double>(r_z, 0.0));
//   vector<vector<double>> p_xy_z(r_x * r_y, vector<double>(r_z, 0.0));
//   vector<vector<double>> p_xz(r_x, vector<double>(r_z, 0.0));
//   vector<vector<double>> p_yz(r_y, vector<double>(r_z, 0.0));
//   //cout << "p_xy_z" << endl;
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         p_xyz.at(i * r_y + j).at(k) = (double)count.at(i * r_y + j).at(k) /
//         (double)n; p_xy_z.at(i * r_y + j).at(k) = (double)count.at(i * r_y +
//         j).at(k) / (double)sum_z.at(k);
//         //cout << p_xy_z.at(i * r_y + j).at(k) << ", ";
//       }
//       //cout << endl;
//     }
//   }
//   //cout << "p_xz" << endl;
//   for (int i = 0; i < r_x; i++) {
//     for (int k = 0; k < r_z; k++) {
//       p_xz.at(i).at(k) = (double)n_xz.at(i).at(k) / (double)sum_x.at(k);
//       //cout << p_xz.at(i).at(k) << ", ";
//     }
//     //cout << endl;
//   }
//   //cout << "p_yz" << endl;
//   for (int j = 0; j < r_y; j++) {
//     for (int k = 0; k < r_z; k++) {
//       p_yz.at(j).at(k) = (double)n_yz.at(j).at(k) / (double)sum_y.at(k);
//       //cout << p_yz.at(j).at(k) << ", ";
//     }
//     //cout << endl;
//   }
//   //cout << "log" << endl;
//   double cmi = 0.0;
//   for (int i = 0; i < r_x; i++) {
//     for (int j = 0; j < r_y; j++) {
//       for (int k = 0; k < r_z; k++) {
//         double temp = 0;
//         temp = p_xyz.at(i * r_y + j).at(k) * log(p_xy_z.at(i * r_y + j).at(k)
//         / (p_xz.at(i).at(k) * p_yz.at(j).at(k))); cmi = cmi + temp;
//         //cout <<"(" << p_xy_z.at(i * r_y + j).at(k) / (p_xz.at(i).at(k)*
//         p_yz.at(j).at(k)) <<"," << log(p_xy_z.at(i * r_y + j).at(k)
//         /(p_xz.at(i).at(k) * p_yz.at(j).at(k)))<<")";
//         //cout << temp << ", ";
//       }
//       //cout << endl;
//     }
//   }

//   //cout << "CMI: " << cmi << endl;
//   if(cmi < threthold){
//     //cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<<node_z<<endl;
//     return false;
//   }else{
//     //cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<<node_z<<endl;
//     return true;
//   }
// }