--check equivalence of two ADMGs
checkequivalencerat = {n,edg1,edg2,b1,b2} -> {
     R = frac(QQ[l_(0,0)..l_(n-1,n-1),b_(0,0)..b_(n-1,n-1),m_(0,0)..m_(n-1,n-1)]);
    
    -- adjacency matrix of the bidiricted part of the first graph
    Bad1 =  mutableIdentity(R,n);
    for e in b1 do(
	Bad1_e = Bad1_(e_1,e_0) = 1;
	);
    
    Bad2 =  mutableIdentity(R,n);
    for e in b2 do(
	Bad2_e = Bad2_(e_1,e_0) = 1;
	);
    
    sum1 = 0;
    sum2 = 0;
    
    for i from 0 to n-1 do (
	for j from 0 to i-1 do(
	    if Bad1_(i,j) != 0 then sum1 = sum1+1;
	    if Bad2_(i,j) != 0 then sum2 = sum2+1;
	    );
	);
    --if(sum1 != sum2) then  return(ideal(1),ideal(1));
    
    
    
    -- adjacency matrix directed part 1st graph
    L1 = mutableIdentity(R,n);
    for e in edg1 do (
	L1_e = -l_e;      
	);
    
     -- adjacency matrix directed part 2nd graph 
    L2 = mutableIdentity(R,n);
    for e in edg2 do (
	L2_e = -m_e;      
	);
    
    -- A matrix for 1\in 2
    res1 = amatrixrat(n, L1, transpose L2);
    A1 = res1_1;
    
    -- A matrix for 2\in 1
    res2 = amatrixrat(n, L2, transpose L1);
    A2 = res2_1;
    
     
    -- Candidate bidiracted graphs
    --B1cand = {};
    B1cand = new MutableList;
    for i from 0 to n-1 do(
	for j from 0 to i-1 do(
	    vcheck = false;
	    for k1 from 0 to n-1 do(
		for k2 from 0 to n-1 do(
		    if(vcheck==false and Bad2_(k1,k2) != 0 and A2_(i,k1) != 0 and A2_(j,k2) != 0) then (
		    	--B1cand = append(B1cand,(i,j));
			B1cand##B1cand = (i, j);
		    	vcheck = true;
		    	)
		    );
		);
	    );
	);
        	
        -- B2cand = {};
	B2cand = new MutableList;	
	
    	for i from 0 to n-1 do(
	    for j from 0 to i-1 do(
	    	vcheck = false;
	    	for k1 from 0 to n-1 do(
		    for k2 from 0 to n-1 do(
		    	if(vcheck==false and Bad1_(k1,k2) != 0 and A1_(i,k1) != 0 and A1_(j,k2) != 0) then (
		    	    -- B2cand = append(B2cand,(i,j));
			    B2cand##B2cand = (i, j);
		    	    vcheck = true;
		    	    )
		    	);
		    );
	    	);
	    );

    -- Check candidate bidirected sets
    Cad1 = mutableIdentity(R,n);
    for e in B1cand do(
	Cad1_e = Cad1_(e_1,e_0) = 1;
	);

    Cad2 = mutableIdentity(R,n);
    for e in B2cand do(
	Cad2_e = Cad2_(e_1,e_0) = 1;
	);    
    
    
    D1 = matrix Cad1 - matrix Bad1;
    D2 = matrix Cad2 - matrix Bad2;

    
    -- eqlist1 = {};
    eqlist1 = new MutableList;
    Q1 = QQ[m_(0,0)..m_(n-1,n-1)];
    S1 = Q1[l_(0,0)..l_(n-1,n-1)];        
    A2 = substitute(A2, S1);
    
    
    for i from 0 to n-1 do(
	for j from 0 to n-1 do(
	    A2 = substitute(A2, m_(i,j) => randomrat(1000000));
	    );	
	);
    S1 = QQ[l_(0,0)..l_(n-1,n-1)];        
    A2 = substitute(A2, S1);    
    
    -- elist1 = {};
    -- eqlist2 = {};
    elist1 = new MutableList;
    eqlist2 = new MutableList;
    Q2 = QQ[l_(0,0)..l_(n-1,n-1)];
    S2 = Q2[m_(0,0)..m_(n-1,n-1)];
    A1 = substitute(A1,S2);
    
    for i from 0 to n-1 do(
	for j from 0 to n-1 do(
	    A1 = substitute(A1, l_(i,j) => randomrat(1000000));
	    );	
	);
    S2 = QQ[m_(0,0)..m_(n-1,n-1)];        
    A1 = substitute(A1, S2);    
    
    -- elist2 = {};
    elist2 = new MutableList;
    
    
    for i from 0 to n-1 do(
	for j from 0 to i do(
	    if D1_(i,j) == 1 then (
		for k1 from 0 to n-1 do(
		    for k2 from 0 to n-1 do(
			if(Bad2_(k1,k2) != 0 and A2_(i,k1) != 0 and A2_(j,k2) != 0) then(
    			    -- eqlist1 = append(eqlist1,A2_(i,k1)*A2_(j,k2));
			    eqlist1##eqlist1 = A2_(i,k1)*A2_(j,k2);
			    -- elist1 =  append(elist1, {{(i,k1),(j,k2)}});			   		    
			    elist1##elist1 = {{(i,k1),(j,k2)}};
			    ) 
			);
		    );
		);
	    
	    if D2_(i,j) == 1 then (
		for k1 from 0 to n-1 do(
		    for k2 from 0 to n-1 do(
			if(Bad1_(k1,k2) != 0 and A1_(i,k1) != 0 and A1_(j,k2) != 0) then(
    	    	    	    -- eqlist2 = append(eqlist2,A1_(i,k1)*A1_(j,k2));
			    eqlist2##eqlist2 = A1_(i,k1)*A1_(j,k2);
			    -- elist2 =  append(elist2, {{(i,k1),(j,k2)}});
			    elist2##elist2 = {{(i,k1),(j,k2)}};
			    ) 
			);
		    );
		);
	    );
	);
    	    	    	    
    	I1 = ideal(toList eqlist1);
	--dim I1
	

	I2 = ideal(toList eqlist2);
	--dim I2
	
	return(I1,I2);
   };

amatrixrat = {n, L, L2} -> {
    -- R = QQ[l_(0,0)..l_(n-1,n-1),b_(0,0)..b_(n-1,n-1)];
    L = matrix(L);
    -- B matrix 
    S = inverse(L)*det(L);
    
    S = transpose S;
    
    B = mutableIdentity(R,n);
    for i from 0 to n-1 do(
	for j from 0 to n-1 do (
	    if S_(i,j)!=0 then B_(i,j)=b_(i,j);
	    );
	);
    
    -- A matrix
    A = (matrix L2)*(matrix S);
    return((matrix L2)*(matrix B), A);
    };

amatrix = {n, L, L2} -> {
    -- R = QQ[l_(0,0)..l_(n-1,n-1),b_(0,0)..b_(n-1,n-1)];
    L = matrix(L);
    -- B matrix 
    S = matrix(mutableIdentity(R,n));
    
    L= L-S;
    
    
    D = matrix(mutableIdentity(R,n));
    
    for i from 1 to n-1 do(
	D = D*L;
	S = S+D;
	);  
    S = transpose S;
    
    B = mutableIdentity(R,n);
    for i from 0 to n-1 do(
	for j from 0 to n-1 do (
	    if S_(i,j)!=0 then B_(i,j)=b_(i,j);
	    );
	);
    
    -- A matrix
    A = (matrix L2)*S;
    return(B, A);
    };



-- generate random rational
randomrat = {M} -> {
    
    num = random(-M, M);
    den = random(-M, M);
    
    if den != 0 then return(num/den);
    
    return(num);
    
    };


-- generate random dag
randag = {n,enum} -> {
    if enum > binomial(n,2) then return("there are too many edges");
    if enum == 0 then return({});
    -- edgelist = {};
    edgelist = new MutableList;
        
    cord = random toList (0..n-1);
    
    for i from 0 to n-2 do(
	for j from i+1 to n-1 do(
	    -- edgelist = append(edgelist, (cord_(i),cord_(j)));
	    edgelist##edgelist = (cord_(i),cord_(j));
	    );
    	);
--    print(edgelist);
    edgelist = (random edgelist);
    edgelist = edgelist_{0..enum-1};
    return(sort edgelist);  
    };




-- generate random bidirected part

ranbid = {n,bnum} -> {
    if bnum > binomial(n,2)  then return("there are too many edges");
    edgelist = {};

    cord = random toList (0..n-1);

    for i from 0 to n-2 do(
	for j from i+1 to n-1 do(
	    edgelist = append(edgelist, (cord_(i),cord_(j)));
	    );
    	);
    
    edgelist = random edgelist;
    edgelist = edgelist_{0..bnum-1}; 
    ordlist = {};
    for edge in edgelist do(
	if edge_1 < edge_0 then ordlist = append(ordlist,(edge_1,edge_0))	
	else ordlist = append(ordlist,edge);
	);
    return(sort ordlist);         
    };
    
--search equivalent graph by sampling
findeqsamp = {n,bnum,cnum,edg1,b1,eqclass} -> {
    for enum from 0 to binomial(n,2) do(
    	for i from 1 to cnum do (
	    edg2 = randag(n,enum);
	    b2 = ranbid(n,bnum);    	  
	    isequal = false;
	    for graph in eqclass do(
		if(edg2 == graph_0 and b2 == graph_1) then isequal = true;				
		);
	    if(isequal == false) then (
	    	ideals =  checkequivalence(n,edg1,edg2,b1,b2);
	    	if (dim ideals_0 != -1 and dim ideals_1 != -1 ) then eqclass=append(eqclass,{edg2,b2});
	    	)
	    );
    	);
    return(eqclass);
    };

-- decide model equivalence of two graphs

model_equivalence = {n, edg1, edg2, b1, b2} -> {
    ideals =  checkequivalencerat(n,edg1,edg2,b1,b2);
    if (dim ideals_0 != -1 and dim ideals_1 != -1 ) then return true;
    return false
    };


--search equivalent graph by exaustive search
findeqsearch = {n,dlist,blist,edg1,b1} -> {
    eqclass = {};
    for edg2 in dlist do(
    	for b2 in blist do (
	    ideals =  checkequivalencerat(n,edg1,edg2,b1,b2);
	    if (dim ideals_0 != -1 and dim ideals_1 != -1 ) then eqclass=append(eqclass,{edg2,b2});
	    );
    	);
    return(eqclass);
    };

 --- list all the dags with n nodes 
 daglist = {n} -> {
     orders = permutations toList (0..n-1);
     dlist = {{}};
     
     for o in orders do(
	 edgelist = {};
	 for i from 0 to n-2 do (
	     for j from i+1 to n-1 do(
		 edgelist = append(edgelist,(o_i,o_j));
		 );
	     );	 
	 edgeperm = permutations edgelist;
    	 for el in edgeperm do(
	     for e from 1 to binomial(n,2) do(
		 dlist = append(dlist, sort el_{0..e-1});
		 dlist = unique dlist;
		 );	     
	     );
	 );
     return(dlist);
     };


-- list all bidirected graphs with n nodes and k edges
bidlist = {n,k} -> {     
     blist = {};
     edgelist = {};
     for i from 0 to n-2 do (
     	 for j from i+1 to n-1 do(
	     edgelist = append(edgelist,(i,j));
	     );
     	 );	 
     edgeperm = permutations edgelist;
     for el in edgeperm do(
	 blist = append(blist, sort el_{0..k-1});
	 blist = unique blist;
	 );
     return(blist);
    }


