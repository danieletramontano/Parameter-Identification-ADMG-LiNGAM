load("model_identification.m2")

-- Testing model equivalence of two graphs
n = 3;

edg1 = {(0,1),(1,2)};
b1 = {(1,2)};

-- equivalent models
edg2 = {(0,1),(0, 2),(1,2)};
b2 = {(1,2)};

me = model_equivalence(n, edg1, edg2, b1, b2)

-- non-equivalent models
edg3 = {(0,1),(0, 2),(1,2)};
b3 = {(0, 2), (1,2)};

me = model_equivalence(n, edg1, edg3, b1, b3)
