#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include "base/DAG.h"
namespace py = pybind11;

// 制約ベース構造学習アルゴリズムのためのCPDAGクラス
// DAGはCPDAGの特別な例であり，DAGの一般化といえるので，CPDAGのみでよいものと思われる
// 一旦，無向辺は二つの有効辺，すなわち x -> y かつ y -> x として表現する
// 制約ベースでは，完全無向グラフから始めて，辺を削除する過程で徐々に疎になっていく
// ここで，CPDAGのデータ構造は隣接行列，隣接リスト，辺リスト等いろいろ考えられる
// (1)隣接行列: [i, j] = 1 ならば i -> j である．密なグラフ向け
// (2)隣接リスト: {i: [a_1, a_2, ...]} ならば i -> a_1, i -> a_2, ... である．疎なグラフ向け
// (3)非隣接リスト: {i: [a_1, a_2, ...]} ならば i -> a_1, i -> a_2, ... でない．密なグラフ向け
// 正直どれがいいかわからないが，完全無向グラフを使うので一旦(1)の隣接行列で実装する
// のちのちアルゴリズムごとに変更するのもアリ
// 適宜，制約ベース構造学習で必要な関数を生やす

// pgmpyとのインターフェースについて
// pgmpyではDAG, CPDAGクラスはnetworkxのDiGraphを継承して楽してる
// networkxのDiGraphのコンストラクタは，隣接リストを引数に取るっぽい
// なので，隣接リストを返す関数get_adj_list()が必要
// こんな感じで使えるはず: nx.DiGraph(dag.get_adj_list())

CPDAGwithAdjMat::CPDAGwithAdjMat(size_t n){
    this->n = n;
    this->adj_mat.resize((n*n + 63) / 64, 1);
}

bool CPDAGwithAdjMat::has_edge(int from, int to) {
    int idx = from * this->n + to;
    return (this->adj_mat[idx / 64] >> (idx % 64)) & 1;
}

void CPDAGwithAdjMat::add_edge(int from, int to) {
    int idx = from * this->n + to;
    if (this->adj_mat[idx / 64] >> (idx % 64) & 1) throw std::runtime_error("Edge already exists");
    this->adj_mat[idx / 64] |= 1 << (idx % 64);
}

void CPDAGwithAdjMat::remove_edge(int from, int to) {
    int idx = from * this->n + to;
    if (!(this->adj_mat[idx / 64] >> (idx % 64) & 1)) throw std::runtime_error("Edge does not exist");
    this->adj_mat[idx / 64] &= ~(1 << (idx % 64));
}

std::unordered_map<int, std::vector<int>> CPDAGwithAdjMat::get_adj_list() {
    std::unordered_map<int, std::vector<int>> _adj_list;
    for (int i = 0; i < this->n; i++) {
        std::vector<int> adj;
        for (int j = 0; j < this->n; j++) {
            if (this->has_edge(i, j)) {
                adj.push_back(j);
            }
        }
        _adj_list[i] = adj;
    }
}