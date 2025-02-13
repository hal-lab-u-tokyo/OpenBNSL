#pragma once

vector<vector<int>> state_count(const vector<vector<uint8_t>> &data,
                                vector<int> &children, vector<int> &parents,
                                vector<int> &n_states, int &parallel);

float natori_independent_score(const vector<vector<uint8_t>> &data, int &node_x,
                               int &node_y, vector<int> &parents,
                               vector<int> &n_states, float &ESS,
                               int &parallel);

float natori_dependent_score(const vector<vector<uint8_t>> &data, int &node_x,
                             int &node_y, vector<int> &parents,
                             vector<int> &n_states, float &ESS, int &parallel);

bool ci_test_by_Bayesfactor(const vector<vector<uint8_t>> &data, int &node_x,
                            int &node_y, vector<int> &Z, vector<int> &n_states,
                            float &ESS, int &parallel);