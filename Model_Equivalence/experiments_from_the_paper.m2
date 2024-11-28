load("model_identification.m2")

-- IV equivalence class
n = 3;

dlist3 = daglist(3);
blist31 = bidlist(3,1);

edg1 = {(0,1),(1,2)};
b1 = {(1,2)};

eqclass = findeqsearch(n,dlist3,blist31,edg1,b1);

for e in eqclass do(
     print(e);
     );

