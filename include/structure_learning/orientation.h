#include "base/PDAG2.h"

/*
    orient edges in a PDAG to a maximally oriented graph.
    orient rules are based on rule 1~3 from Meek,C.:Causal Inference and
   Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in
   Artificial Inteligence (UAl-95),p.403-410 (195)
*/
void orientation(PDAG &G, const vector<int> &sepsets);
void orientation(int level, PDAG &G, const vector<int> &sepsets);