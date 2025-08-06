"""
Functions for identifying the causal effect vector in a causal graph.
"""

import random
import pickle
import os
import warnings
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import multiprocess as mp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

warnings.filterwarnings("ignore", category=FutureWarning, message="adjacency_matrix will return a scipy.sparse array instead of a matrix in Networkx 3.0.")


def find_rv(g_dir, g_bid, v):
    """
    Find the removable ancestors of a node v
    Input:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.
    - v: Vertex for which rv needs to be calculated.

    Output:
    - rv: List of vertices that form rv based on the algorithm.
    """
    # Find the ancestors of vertex v in the directed graph
    av = nx.ancestors(g_dir, v)
    rv = []

    # Create a set of neighbors of v in the bidirectional graph
    sv = set(g_bid.neighbors(v))
    sv.add(v)  # Include v itself in the set

    # Iterate over the ancestors in av
    for _, u in enumerate(av):
        # Create a set of neighbors of u in the bidirectional graph
        su = set(g_bid.neighbors(u))
        su.add(u)  # Include u itself in the set

        # Check if Su is not a subset of Sv, and if not, add u to rv
        if not su.issubset(sv):
            rv.append(u)

    return rv



def augment_graph(adj_matrix, capacities):
    """
    Augment the graph to accomodate max flow with nodes capacity
    Input:
    - adj_matrix: The adjacency matrix of the original graph.
    - capacities: A list of capacities for each vertex in the original graph.

    Output:
    - Augmented adjacency matrix for solving the max-flow problem.
    """
    num_vertices = len(adj_matrix)
    augmented_matrix = np.zeros((2 * num_vertices, 2 * num_vertices), dtype=int)

    # Iterate over the original adjacency matrix
    for i in range(num_vertices):
        for j in range(num_vertices):
            if adj_matrix[i, j] > 0:
                augmented_matrix[i, j + num_vertices] = adj_matrix[i, j]  # v_in to v_out

    for i in range(num_vertices):
        augmented_matrix[i + num_vertices, i] = capacities[i]  # v_out to v_i with capacity

    return augmented_matrix



def check_criterion(g_dir, g_bid, v, q):
    """
    Check identifiability of the causal effect vector from q to v
    Input:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.
    - v: Vertex for which the criterion is checked.
    - q: A set of vertices to be checked against the parent set of v.

    Output:
    - True if the criterion is satisfied, False otherwise.
    """
    p = len(g_dir.nodes())

    # Find the parent set of vertex v in the directed graph
    pav = set(g_dir.predecessors(v))

    # Check if q is not a subset of pav, and if not, add it to the parent set
    if not q.issubset(pav):
        q = q & pav

    if len(q) == 0:
        return True

    # Find the ancestors of v in the directed graph
    av = nx.ancestors(g_dir, v)

    # Find the set rv using the find_rv function
    rv = find_rv(g_dir, g_bid, v)

    # Create the adjacency matrix of g_dir
    adj = nx.adjacency_matrix(g_dir)

    # Calculate the max-flow values for q and pav using the maxflow_q function
    val_q = maxflow_q(adj, av, rv, v, q, p)
    val_pa = maxflow_q(adj, av, rv, v, pav.difference(q),p)
    # Calculate the difference in max-flow values and check against the size of q
    val = val_q - val_pa
    return val == len(q)

def maxflow_q(adj, av, rv, v, q, p):
    """
    Calculate max flow based on given parameters
    Input:
    - adj: The adjacency matrix of the augmented graph.
    - av: Set of ancestors of vertex v in the directed graph.
    - rv: List of vertices forming rv for vertex v.
    - v: The source vertex.
    - q: The set of vertices forming the sink set.
    - pav: The parent set of vertex v in the directed graph.

    Output:
    - Maximum flow value for the given parameters.
    """
    newpar = np.zeros(p)
    newpar[list(q)] = 1
    adj[:, v] = newpar
    adj = adj.todense()

    schil = np.zeros(p)
    schil[rv] = 1

    a = np.insert(adj, 0, schil, 0)
    a = np.insert(a, 0, np.zeros(p + 1), 1)

    av = [int(x) + 1 for x in av]
    rv = [int(x) + 1 for x in rv]
    q = [int(x) + 1 for x in q]
    v = v + 1

    mf_vert = np.insert(av, 0, [0], 0)
    mf_vert = np.insert(mf_vert, len(av), [v], 0)

    a = a * (len(q) + 1)
    a = a[mf_vert,:][:,mf_vert]

    mf_cap = np.ones(len(mf_vert) + 1)
    mf_cap[0] = (len(q) + 1)
    mf_cap[len(mf_cap) - 1] = (len(q) + 1)
    aug = augment_graph(a, mf_cap)

    graph = csr_matrix(aug)
    val = maximum_flow(graph, 0,  aug.shape[0]-1).flow_value

    return val

def wholematrix_criterion(g_dir,g_bid):
    """
    Checks the identifiability criterion for the entire matrix of the causal graph.

    Args:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.

    Returns:
    - True if the entire matrix satisfies the identifiability criterion, False otherwise.
    """
    for v in range(len(g_dir.nodes())):
        if not check_criterion(g_dir,g_bid,v,set(g_dir.predecessors(v))):
            return False
    return True


def rand_graphs_sample(n = 25, m = 50):
    """
    Generates a random ADMG.

    Args:
    - n: Number of nodes in the graph. Default is 25.
    - m: Number of edges in the graph. Default is 50.

    Returns:
    - Directed and bidirectional graphs representing the original graph.
    """
    m_d = random.randint(0,min(m, n*(n-1)/2))
    g_dir = nx.gnm_random_graph(n, m_d, directed=False)
    a_dir = nx.adjacency_matrix(g_dir)
    for i in range(n):
        for j in range(i):
            a_dir[i,j] = 0
    g_dir = nx.DiGraph(a_dir.todense())

    m_b = m-m_d
    g_bid = nx.gnm_random_graph(n, m_b, directed=False)
    return g_dir, g_bid

def rand_graphs_exp(n = 25, ran = 100, rep = 100):
    """
    Simulates random graphs and calculates the proportion of identifiable ADMGs.

    Args:
    - n: Number of nodes in the graph. Default is 25.
    - ran: Number of iterations. At each iterations the number of edges is the graph
    increases of a factor of n.
    - rep: Number of repetitions per iteration. Default is 5000.

    Returns:
    - Array containing the proportion of identifiable ADMGs for each iteration.
    """
    per = np.zeros(ran)
    m_max = n*(n-1)
    for i, k in enumerate(np.linspace(0.01, 0.99, ran)):
        m = int(m_max*k)
        t = 0
        for _ in range(rep):
            g_dir, g_bid = rand_graphs_sample(n, m)
            t += wholematrix_criterion(g_dir,g_bid)
        per[i] = t/rep
    return per

def rand_graphs_exp_par(n, ran = 100, rep = 625):
    """
    Simulates random graphs in parallel and calculates tthe proportion of identifiable ADMGs.

    Args:
    - n: Number of nodes in the graph.
    - ran: Number of iterations. At each iterations the number of edges is the graph
    increases of a factor of n.
    - rep: Number of repetitions per iteration. Default is 625.

    Returns:
    - Array containing the count of identifiable ADMGs for each iteration.
    """
    counts = np.zeros(ran)
    m_max = n*(n-1)
    for i, m in enumerate(np.linspace(0.01, 0.99, ran)*m_max):
        t = 0
        for _ in range(rep):
            g_dir, g_bid = rand_graphs_sample(n, int(m))
            t += wholematrix_criterion(g_dir,g_bid)
        counts[i] = t
    return counts

def parallel_rand_graphs_exp(n):
    """
    Executes random graph simulations in parallel.

    Args:
    - n: Number of nodes in the graph.

    Returns:
    - Sum of results from parallel simulations.
    """
    pool = mp.Pool(8)

    results = pool.map(rand_graphs_exp_par, [n]*8)

    pool.close()

    return sum(results)

def save_object(obj):
    """
    Saves the object to a pickle file.

    Args:
    - obj: Object to be saved.
    """
    try:
        with open(obj.head+".pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)

def load_object(filename):
    """
    Loads the object from a pickle file.

    Args:
    - filename: Name of the pickle file.

    Returns:
    - Loaded object.
    """
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported):", ex)


class Countdata:
    """
     Represents data containing counts of identifiable ADMGs.

     Attributes:
     - n: Number of nodes in the graph.
     - rep: Number of repetitions per iteration.
     - count: Array containing the counts of identifiable ADMGs for each iteration.
     - edges: Array containing the number of edges for each iteration.
     - head: Head of the data.
     """
    def __init__(self, n, reps, counts, ran):
        self.n = n
        self.rep = reps
        self.count = counts
        self.edges = [n*(m+1) for m in range(ran)]
        self.head = f"Data/id_count-n_ {n} -medg_ {n*ran}-rep_{reps}"

def prod_plot_mult_plots(files = None,
                         files_dir = None,
                         label_fontsize=17,
                         subtitle_fontsize=25,
                         legend_fontsize=18):
    """
    Produces multiple bar plots showing the proportion of identifiable and non-identifiable ADMGs.

    Args:
    - files (list): A list of Countdata objects containing data for plotting. If None, files_dir must be provided.
    - files_dir (str): Directory path containing Countdata objects. Default is None.
    Either the file list or the directory must be provided.
    - label_fontsize: Font size for the axis labels. Default is 17.
    - subtitle_fontsize: Font size for the subplot titles. Default is 25.
    - legend_fontsize: Font size for the legend. Default is 18.
    """

    if files is None:
        if files_dir is None:
            raise ValueError('Specify location for the files, or enter file list.')
        file_names = [f for f in os.listdir(files_dir) if not os.path.isfile(f) and f[0] == "i"]
        files = []
        for file_name in file_names:
            files = [(load_object(files_dir + file_name))] + files

    fig, axs = plt.subplots(1, len(files), figsize=(12, 8), sharey=True)

    for i, f in enumerate(files):
        counts = np.concatenate(([0], f.rep - f.count, [f.rep]))
        n = f.n
        x_axis = np.concatenate(([0], np.linspace(0.01, 0.99, len(f.count)), [1]))

        ax = axs[i]
        # Plot results
        ax.plot(x_axis, counts, color='gray', linewidth=2)
        ax.fill_between(x_axis, counts, label = 'Non-Identifiable', color='black', alpha=0.4)
        ax.fill_between(x_axis, counts, f.rep, label = 'Identifiable', color='gray', alpha=0.3)

        ax.set_title(f'p = {n}', fontsize=subtitle_fontsize)
        ax.set_xlabel('Density', fontsize=label_fontsize)
        ax.set_facecolor('white')  # Set background color

        ax.tick_params(axis='both', which='major', labelsize=label_fontsize)
        ax.spines['top'].set_color('black')
        ax.spines['right'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.xaxis.set_tick_params(color='black')
        ax.yaxis.set_tick_params(color='black')

    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()  # Reverse the handles list
    labels.reverse()   # Reverse the labels list

    fig.legend(handles,
               labels,
               loc='lower center',
               bbox_to_anchor=(0.5, -0.05),
               fontsize=legend_fontsize,
               ncol=2)  # Display legend items in 2 columns

    plt.subplots_adjust(wspace=0.05)

    plt.show()