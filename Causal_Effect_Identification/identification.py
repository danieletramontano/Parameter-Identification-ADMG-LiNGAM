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


def find_Rv(g_dir, g_bid, v):
    """
    Find the removable ancestors of a node v
    Input:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.
    - v: Vertex for which Rv needs to be calculated.

    Output:
    - Rv: List of vertices that form Rv based on the algorithm.
    """
    # Find the ancestors of vertex v in the directed graph
    Av = nx.ancestors(g_dir, v)
    Rv = []
    
    # Create a set of neighbors of v in the bidirectional graph
    Sv = {n for n in g_bid.neighbors(v)}
    Sv.add(v)  # Include v itself in the set
    
    # Iterate over the ancestors in Av
    for i, u in enumerate(Av):
        # Create a set of neighbors of u in the bidirectional graph
        Su = {n for n in g_bid.neighbors(u)}
        Su.add(u)  # Include u itself in the set
        
        # Check if Su is not a subset of Sv, and if not, add u to Rv
        if not Su.issubset(Sv):
            Rv.append(u)

    return Rv



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



def check_criterion(g_dir, g_bid, v, Q):
    """
    Check identifiability of the causal effect vector from Q to v
    Input:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.
    - v: Vertex for which the criterion is checked.
    - Q: A set of vertices to be checked against the parent set of v.

    Output:
    - True if the criterion is satisfied, False otherwise.
    """
    
    p = len(g_dir.nodes())
    
    # Find the parent set of vertex v in the directed graph
    pav = {n for n in g_dir.predecessors(v)}
    
    # Check if Q is not a subset of pav, and if not, raise a ValueError
    if not Q.issubset(pav):
        Q = Q & pav
    
    if len(Q)==0:
        return True
    
    # Find the ancestors of v in the directed graph
    Av = nx.ancestors(g_dir, v)
    
    # Find the set Rv using the find_Rv function
    Rv = find_Rv(g_dir, g_bid, v)
    
    # Create the adjacency matrix of g_dir
    Adj = nx.adjacency_matrix(g_dir)
    
    # Calculate the max-flow values for Q and pav using the maxflow_Q function
    val_Q = maxflow_Q(Adj, Av, Rv, v, Q, pav,p)
    val_pa = maxflow_Q(Adj, Av, Rv, v, pav.difference(Q), pav,p)
    
    # Calculate the difference in max-flow values and check against the size of Q
    val = val_Q - val_pa
    if val == len(Q):
        return True
    else:
        return False

def maxflow_Q(Adj, Av, Rv, v, Q, pav,p):
    """
    Calculate max flow based on given parameters
    Input:
    - Adj: The adjacency matrix of the augmented graph.
    - Av: Set of ancestors of vertex v in the directed graph.
    - Rv: List of vertices forming Rv for vertex v.
    - v: The source vertex.
    - Q: The set of vertices forming the sink set.
    - pav: The parent set of vertex v in the directed graph.

    Output:
    - Maximum flow value for the given parameters.
    """
    newpar = np.zeros(p)
    newpar[list(Q)] = 1
    Adj[:, v] = newpar
    Adj = Adj.todense()
    
    schil = np.zeros(p)
    schil[Rv] = 1
    
    A = np.insert(Adj, 0, schil, 0)
    A = np.insert(A, 0, np.zeros(p + 1), 1)

    Av = [int(x) + 1 for x in Av]
    Rv = [int(x) + 1 for x in Rv]
    Q = [int(x) + 1 for x in Q]
    v = v + 1

    mf_vert = np.insert(Av, 0, [0], 0)
    mf_vert = np.insert(mf_vert, len(Av), [v], 0)
    
    A = A * (len(Q) + 1)
    A = A[mf_vert,:][:,mf_vert]
    
    mf_cap = np.ones(len(mf_vert) + 1)
    mf_cap[0] = (len(Q) + 1)
    mf_cap[len(mf_cap) - 1] = (len(Q) + 1)
    Aug = augment_graph(A, mf_cap)

    graph = csr_matrix(Aug)
    val = maximum_flow(graph, 0,  Aug.shape[0]-1).flow_value
    
    return val

def wholematrix_crit(g_dir,g_bid):
    """
    Checks the identifiability criterion for the entire matrix of the causal graph.
    
    Args:
    - g_dir: Directed graph (networkx.DiGraph) representing the original graph.
    - g_bid: Bidirectional graph (networkx.Graph) representing the original graph.
    
    Returns:
    - True if the entire matrix satisfies the identifiability criterion, False otherwise.
    """
    for v in range(len(g_dir.nodes())):
        if check_criterion(g_dir,g_bid,v,{n for n in g_dir.predecessors(v)}) == False:
            return False
    return True


def ranmixdigraph(n = 25, m = 50):
    """
    Generates a random ADMG.
    
    Args:
    - n: Number of nodes in the graph. Default is 25.
    - m: Number of edges in the graph. Default is 50.
    
    Returns:
    - Directed and bidirectional graphs representing the original graph.
    """
    m_d = random.randint(0,m)
    g_dir = nx.gnm_random_graph(n, m_d, directed=False)
    Adir = nx.adjacency_matrix(g_dir)
    for i in range(n):
        for j in range(i):
            Adir[i,j] = 0
    g_dir=nx.DiGraph(Adir.todense())
    
    m_b = m-m_d
    g_bid = nx.gnm_random_graph(n, m_b, directed=False)
    return g_dir, g_bid

def rangraphsimul(n = 25, ran = 10, rep = 5000):
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
    per=np.zeros(ran)
    for k in range(ran):
        m = n*(k+1)
        t = 0
        for i in range(rep):
            g_dir, g_bid = ranmixdigraph(n, m)
            if wholematrix_crit(g_dir,g_bid)==True:
                t = t+1
        per[k]=t/rep
    return per

def rangraphsimul_par(n, ran = 10, rep = 625):
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
    counts=np.zeros(ran)
    for k in range(ran):
        m = n*(k+1)
        t = 0
        for i in range(rep):
            g_dir, g_bid = ranmixdigraph(n, m)
            if wholematrix_crit(g_dir,g_bid)==True:
                t = t+1
        counts[k]=t
    return counts

def parallel_ranggraph(n):
    """
    Executes random graph simulations in parallel.
    
    Args:
    - n: Number of nodes in the graph.
    
    Returns:
    - Sum of results from parallel simulations.
    """
    pool = mp.Pool(8)
    nlist = np.floor(np.ones(8)*n)
    
    results = pool.map(rangraphsimul_par, [int(n) for n in nlist])
    
    # Step 3: Don't forget to close
    pool.close()    
    
    return sum(results)
    
def prod_plot(counts=np.zeros(10),n=25, ran=10, reps=8000, w=22):
    """
    Produces a bar plot showing the proportion of identifiable and non-identifiable ADMGs.

    Args:
    - counts: Array containing the counts of identifiable ADMGs for each iteration.
    - n: Number of nodes in the graph. Default is 25.
    - ran: Number of iterations. Default is 10.
    - reps: Number of repetitions per iteration. Default is 8000.
    - w: Width of the bars. Default is 22.
    """
    type_counts = {
        'Non-Idenitifiable': reps-counts,
        'Identifiable': counts,
        }
    nedges = [n*(m+1) for m in range(ran)]

    width = w # the width of the bars: can also be len(x) sequence


    fig, ax = plt.subplots()
    bottom = np.zeros(10)

    for t,t_count in type_counts.items():
        p = ax.bar(nedges, t_count, width, label=t, bottom=bottom)
        bottom += t_count

        ax.bar_label(p, label_type='center')

    # ax.set_title('Proportion of Identifiable ADMGs with {} nodes'.format(n))
    ax.set_xlabel('Number of Edges')
    ax.legend()

    plt.show()

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


class count_data:
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
        self.head = "Data/id_count-n_ {} -medg_ {}-rep_{}".format(n, n*ran, reps)
        
class minim_data:
    """
    Represents minimal data.
    
    Attributes:
    - df: Dataframe containing the data.
    - head: Head of the data.
    """
    def __init__(self, data, head):
        self.df = data
        self.head = head
    
def prod_plot_mult_plots(files = None,
                         files_dir = None,
                         title_fontsize=20,
                         label_fontsize=13, 
                         subtitle_fontsize=20,
                         legend_fontsize=18,
                         inner_fontsize=11):
    """
    Produces multiple bar plots showing the proportion of identifiable and non-identifiable ADMGs.
    
    Args:
    - files (list): A list of count_data objects containing data for plotting. If None, files_dir must be provided.
    - files_dir (str): Directory path containing count_data objects. Default is None. 
    Either the file list or the directory must be provided.
    - title_fontsize: Font size for the main title. Default is 16.
    - label_fontsize: Font size for the axis labels. Default is 13.
    - subtitle_fontsize: Font size for the subplot titles. Default is 20.
    - legend_fontsize: Font size for the legend. Default is 13.
    - inner_fontsize: Font size for the bar labels. Default is 11.
    """
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), sharey=True)  
    ran = 10
    
    if files is None:
        if files_dir is None:
            raise ValueError('Specify location for the files, or enter file list.')
        else:
            file_names = [f for f in os.listdir(files_dir) if not os.path.isfile(f) and f[0] == "i"]
            files = []
            for file_name in file_names:
                files = [(load_object(files_dir + file_name))] + files                

    for i, f in enumerate(files):
        counts = f.count
        n = f.n
        reps = f.rep
        if n == 25:
            w = 22
        else:
            w = 44
        type_counts = {
            'Non-Identifiable': reps - counts,
            'Identifiable': counts,
        }
        nedges = [n * (m + 1) for m in range(ran)]

        width = w

        ax = axes[i]
        
        bottom = np.zeros(10)

        for t, t_count in type_counts.items():
            color = 'black' if t == 'Non-Identifiable' else 'grey'
            alpha =  0.4 if t == 'Non-Identifiable' else 0.3
            p = ax.bar(nedges, t_count, width, label=t, bottom=bottom, color=color, edgecolor='black', alpha = alpha)

            bottom += t_count

            # ax.bar_label(p, label_type='center', fontsize=inner_fontsize, color='black')

        ax.set_title('p = {}'.format(n), fontsize=subtitle_fontsize)
        ax.set_xlabel('Number of Edges', fontsize=label_fontsize)
        ax.set_facecolor('white')  # Set background color
        
        # Adjusting tick label sizes
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