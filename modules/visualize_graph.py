import matplotlib.pyplot as plt
import networkx as nx

# 有向グラフの描画
def display_graph_info(G):
    # # ノードのリストを表示
    # print("Nodes in the graph:")
    # print(G.nodes())

    # # エッジのリストを表示
    # print("\nEdges in the graph:")
    # print(G.edges())

    # 各ノードの隣接ノードを表示
    if nx.is_directed(G):
        print("\nNeighbors of each node:")
        for node in G.nodes():
            print(f"{node}: {list(G.successors(node))}")
    if nx.is_directed(G):
        Gc = nx.DiGraph()
    else:
        Gc = nx.Graph()
    Gc.add_nodes_from(G.nodes())
    Gc.add_edges_from(G.edges())
    pos = nx.spring_layout(Gc)  # ノードの配置を計算
    #print("pos:", pos)
    # nx.draw(Gc, pos, with_labels=True, node_size=3000, node_color='lightblue', font_size=12, font_weight='bold', arrowstyle='-|>', arrowsize=15)
    # plt.title("Directed Graph")
    # plt.show()
