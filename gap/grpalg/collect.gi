#############################################################################
##
#W    collect.gi                 From the ModIsomExt package           Leo Margolis
#W                                                                     Tobias Moede
##
## functions to work with tables (in the sense of ModIsom) which correspond to quotients of the augmentation ideal I of a group algebra of a finite p-group over a field of characteristic p modulo a power of I. In particular, to generate the table for this quotient

# Input: A list
# Output: The nonzero positions of the list
BindGlobal("PosNonzero", function(list)
    return PositionsProperty(list, x -> not IsZero(x));
end);


# Input: Linear combination of words. A word is a list of pairs of natural numbers. We think e.g. g1^2*g3 = [[1,2],[3,1]]
# Output: Same linear combination (in the math sense) where each word 
#         appears at most once
BindGlobal("ReducedSumOfWords", function( s )
local todo, ret, i, j, sc;

	sc := StructuralCopy(s);
	todo := [1..Length(sc)];
	ret := [];

	for i in todo do
		if i = false then
			continue;
		fi;

		for j in [i+1..Length(sc)] do
			if sc[i][2] = sc[j][2] then
				sc[i][1] := sc[i][1] + s[j][1];
				todo[j] := false;
			fi;
		od;

		ret[Length(ret)+1] := i;
	od;

	return Filtered(sc{ret}, x -> x[1] <> 0*x[1]);
end);

# If the same Jennings basis elements appear one after the other they get multiplied, i.e. (g-1)(h-1)(h-1)(k-1) -> (g-1)(h-1)^2(k-1)
# Reduces the word w
# Input: Jennings-word
# Output: Same word where two eqaul factors don't appear one after the other
BindGlobal("ReducedWord", function( w )
local i, wc;

	wc := StructuralCopy(w);
	
	for i in [1..Length(wc)-1] do
		if wc[i][1] = w[i+1][1] then
			wc[i][2] := wc[i][2] + wc[i+1][2];
			wc[i+1][1] := false;
		fi;
	od;

	return Filtered(wc, x -> x[1] <> false);
end);

# A word is sorted if it corresponds to a Jennings basis element, i.e. the basic factors (g-1) appear in the same order as in the given basis. No exponent can exceed p-1
# Check whether a reduced word w is a sorted word wrt. a prime p
# Input: Jennings-word and characteristic of underlying group algebra
# Output: Boolean whether the word is a basis element
BindGlobal("IsSortedWord", function( w, p )
local i;

	# Check exponents in [0,..,p-1]
	for i in [1..Length(w)] do
		if w[i][2] > p-1 then
			return false;
		fi;
	od;

	# Check if gen indices ascending
	for i in [1..Length(w)-1] do
		if w[i][1] >= w[i+1][1] then
			return false;
		fi;
	od;

	return true;
end);


# Expresses the commutators of pairs of Jennings pcgs and their p-th powers as products in the 
# Jennings pcgs
# Input: record containing pcgs of a Jennings-series as returned by PcgsJenningsSeries
# Output: Record with lists of exponents corresponding to the pcgs
BindGlobal("ExpsOfComPow", function(s)
local coms, pows, p, combs, c, g, h, i;

    coms := [ ];
    pows := [ ];
    p := PrimePGroup(Group(s.pcgs));

    combs := Combinations([1..Length(s.pcgs)], 2);

    for c in combs do
        if c[2] > c[1] then
		    g := s.pcgs[c[1]];
		    h := s.pcgs[c[2]];

            Add(coms, [c, ExponentsOfPcElement(s.pcgs, h^-1*g^-1*h*g)]);
	    fi;
    od;

    for i in [1..Length(s.pcgs)] do
	    g := s.pcgs[i];
	    Add(pows, ExponentsOfPcElement(s.pcgs, g^p));
    od;

    return rec( coms := coms, pows := pows);
end);


# Gets all the relevant Jennings basis information of the modular group algebra A
# Input: Modular group algebra
# Output: Record with Jennings basis data
BindGlobal("PrecompWeightedBasisOfRad", function( A )
local G, p, n, js, jb, jw, ww, wb, wc, df, i, h, ecp;

    # set up
    G := UnderlyingMagma(A);
    p := PrimePGroup(G);
    n := Length(Pcgs(G));

    # get a special basis of G
    js := PcgsJenningsSeries(G);
    ecp := ExpsOfComPow(js);
    js.coms := ecp.coms;
    js.pows := ecp.pows;    

    jb := js.pcgs;
    jw := js.weights;
    df := List( jb, x -> (x-One(A)) );

    # set up for basis with weights
    wc := Tuples( [ 0 .. p-1 ], n ); RemoveSet( wc, 0*[1 .. n] );
    wc := List(wc, Reversed);
    ww := [];
    wb := [];

    # determine weights for a basis of A
    for i in [1..Length(wc)] do
        ww[i] := Sum( [1..n], x -> wc[i][x]*jw[x] );
    od;

    # sort and return
    h := Sortex(ww);
    return rec( jen := js,
		 	df := df,
 		        weights := ww,
                        exps := Permuted(wc, h),
		        poswone := Positions(List(Permuted(wc, h), Sum),1),
		        rnk := RankByWeights(ww) );
end);


# Compute weight of given word.
# Input: Jennings word
# Output: Its weight
BindGlobal("WordWeight", function(T, w)
local weights, g, b;

    weights := T.pre.jen.weights;
    g := 0;

    for b in w do
        g := g + weights[b[1]]*b[2];
    od;

    return g;
end);


# Translates word to the corresponding exponent. E.g. if |G| = p^4: [[2,2],[3,1]] -> [0,2,1,0]
# Input: Underlying table and reduced word
# Output: Exponent corresponding to the word, a list of integers length n for |G| = p^n  
BindGlobal("SortedWordToExp", function(T, w)
local l, L, b;

    l := Length(T.pre.exps[1]);
    L := ListWithIdenticalEntries(l, 0);

    for b in w do
        L[b[1]] := b[2];
    od;

    return L;
end);


# Translation reversing the preceeding function
BindGlobal("ExpToWord", function(exp)
local w, sw, pos, p;
  
    sw := 0;
    w := [ ];
    pos := PosNonzero(exp);
    
    for p in pos do
        sw := sw + 1;
        w[sw] := [p, exp[p]];
    od;

    return w;
end);


# Translation of linear combination of words to a vector over underlying field with respect to Jennings basis
# Input: Table, Linear combination of words
# Output: Vector of length Size(Jennings basis) = |G| - 1
BindGlobal("LinComSortWordToVec", function(T, lw)
local F, v, exp, pos, s;

    F := T.fld;
    v := ListWithIdenticalEntries(T.dim, Zero(F));

    for s in lw do
        exp := SortedWordToExp(T, s[2]);
        pos := Position(T.pre.exps, exp);
        v[pos] := v[pos] + s[1];
    od;

    return v;
end);


# Sort linear combination of words into reduced (done) and not reduced (todo)
# Input: Linear combination of words, characteristic of field
# Output: Record with two lists
BindGlobal("TodoDoneLinComWords", function(lw, p)
local todo, done, sd, st, s;

    todo := [ ];
    st := 0;
    done := [ ];
    sd := 0;

    for s in lw do
        if IsSortedWord(s[2], p) then
            sd := sd + 1;
            done[sd] := s;
        else
            st := st + 1;
            todo[st] := s;
        fi;
    od;

    return rec(todo := todo, done := done );
end);


# Distributively multiply a word on the left, a word on the right and a linear combination 
# of words with coeffcients in the middle
# Input: Table, word which is factor on left, linear combination in the middle, word factor on right, 
#        a factor from the field with which the result is multiplied, the maximal weight of the words 
#        we want
# Output: Linear combination of words. All summands have weight <= wmax
BindGlobal("MultDistWords", function(T, lw, rep, rw, fac, wmax) 
local res, sres, slw, t, respart, srespart, s;

    res := [ ];
    sres := 0;
    slw := Size(lw);

    # Each summand in rep will contribute one word in the linear comb
    for t in rep do 
        sres := sres + 1;

        # Get factor of the summand
        respart := [ t[1]*fac, StructuralCopy(lw) ]; 
        srespart := slw;
  
        # Write the word from the middle behind lw
        for s in t[2] do 
            srespart := srespart + 1;
            respart[2][srespart] := s;
        od;
 
        # Add the word on the right
        for s in rw do  
            srespart := srespart + 1;
            respart[2][srespart] := s;
        od;
  
        # Check if resulting word has weight small enough
        if WordWeight(T, respart[2]) <= wmax then 
            respart[2] := ReducedWord(respart[2]);
            sres := sres + 1;
            res[sres] := respart;
        fi;  
    od; 

    return res;
end);


# Transforms a linear combination of words into an expression in the Jennings basis
# Input: Table, linear comb of words
# Output: Same linear comb expressed in Jennings basis
BindGlobal("TransLinComToBas", function(T, lcw) # Eigabe muss reduziert sein in der momentanen Implementierung
local p, te, todo, done, sd, wmax, todotemp, stdt, s, lw, rep, rw, slw, srw, con, i, j, mult, d, t;

    p := Characteristic(T.fld);
    te := TodoDoneLinComWords(lcw, p);
    todo := te.todo;
    done := te.done;
    sd := Size(done);
    wmax := T.pre.weights[T.dim];

    while todo <> [ ] do
        # The non-reduced word involved in the product
        todotemp := [ ]; 
        stdt := 0;

        for s in todo do
            # left factor word
            lw := [ ]; 
            # replaced word
            rep := [ ]; 
            # right factor word
            rw := [ ]; 
            # sizes of left and right words
            slw := 0; 
            srw := 0;

            con := false;

            for i in [1..Size(s[2])] do
                if con then 
                    continue;
                fi;
              
                # First found a p-th power
                if s[2][i][2] >= p then 
                    # Put rest of power to left side
                    if s[2][i][2] > p then 
                        slw := slw + 1;
                        lw[slw] := [s[2][i][1], s[2][i][2]-p];
                    fi;
        
                    for j in [i+1..Size(s[2])] do
                        srw := srw + 1;
                        rw[srw] := [s[2][j][1], s[2][j][2]];
                    od;

                    # Part to be replaced
                    rep := T.powwords[s[2][i][1]]; 
                    mult := MultDistWords(T, lw, rep, rw, s[1], wmax);
                    te := TodoDoneLinComWords(mult, p);
 
                    for d in te.done do
                        sd := sd + 1;
                        done[sd] := d;
                    od;

                    for t in te.todo do
                        stdt := stdt + 1;
                        todotemp[stdt] := t;
                    od;

                    # Skip rest of the word
                    con := true; 
                fi; 

                if con then 
                    continue;
                fi;

                # if (g_i - 1)(g_j - 1) with i > j is found
                if i < Size(s[2]) and s[2][i][1] > s[2][i+1][1] then 
                    if s[2][i][2] > 1 then
                        slw := slw + 1;
                        lw[slw] := [s[2][i][1], s[2][i][2]-1];
                    fi;
        
                    if s[2][i+1][2] > 1 then 
                        srw := srw + 1;
                        rw[srw] := [s[2][i+1][1], s[2][i+1][2]-1];
                    fi;
 
                    for j in [i+2..Size(s[2])] do
                        srw := srw + 1;
                        rw[srw] := [s[2][j][1], s[2][j][2]];
                    od;
        
                    rep := T.commwords[s[2][i][1]][s[2][i+1][1]];
                    mult := MultDistWords(T, lw, rep, rw, s[1], wmax);
                    te := TodoDoneLinComWords(mult, p);

                    for d in te.done do
                        sd := sd + 1;
                        done[sd] := d;
                    od;

                    for t in te.todo do
                        stdt := stdt + 1;
                        todotemp[stdt] := t;
                    od;

                    con := true;
                fi; 

                if con then 
                    continue;
                fi;

                slw := slw + 1;

                # Nothing to be replaced was found at this moment
                lw[slw] := [s[2][i][1], s[2][i][2]]; 
            od;
        od;
        
        todo := todotemp;
    od;

    done := ReducedSumOfWords(done);
    
    return done;
end);


##################################################################################################
# Computes (g_i-1)^p as expression in Jennings basis
# Input: Table and maximal weight we want
# Output: None. T.powwords is filled.
BindGlobal("WordFillTablePowers", function(T, lvl)
local F, p, exp, pos1, i, j, pows, pos, combs, v, c, w, s, m, tup, coefprod, expprod, mm, vec;

    F := T.fld;
    p := Characteristic(T.fld);

    for i in T.pre.poswone do
        exp := T.pre.exps[i];
        pos1 := Position(exp, 1);
         # Fill first the p-th power of (g_i-1) for g_i a Jennings pcgs
        j := Position(T.pre.exps, (p-1)*exp); 
        
        # Check if weight of p-th power exceeds level
        if p*T.pre.weights[i] <= lvl then 
            pows := StructuralCopy(T.pre.jen.pows[pos1]);
            pos := PosNonzero(pows);
            pows := pows{pos};
            combs := Combinations([1..Length(pos)]);
            Remove(combs, 1); 
            v := [ ];
    
            # Sum over all non-trivial subsets of Jennings series pcgs involved in expression of g_i^p
            for c in combs do 
                w := [ ];
                s := Length(c);
                m := Maximum(pows{c});
      
                # For each subset sum over all tuples of exponents not exceeding the exponents involved in g_i^p
                for tup in IteratorOfTuples([1..m], s) do 
                    # Check that exponents are not too big
                    if tup <= pows{c} then
                        # Check if the weight of the summand is still small enough
                        if Sum(List([1..s], x -> tup[x]*T.pre.jen.weights[pos[c[x]]])) <= lvl then 
                            coefprod := Product([1..s], x-> Binomial(pows[c[x]], tup[x] ));
                            expprod := ListWithIdenticalEntries(Length(exp), 0);
           
                            for mm in [1..s] do
                                expprod[pos[c[mm]]] := tup[mm];
                            od;

                            vec := ExpToWord(expprod);
                            Add(w, [coefprod*One(F), vec]);
                        fi;
                    fi;
                od; 
      
                Append(v, w);
            od;

            T.powwords[pos1] := v;

        # Weight too large, so product is trivial
        else 
            T.powwords[pos1] := [ ]; 
        fi;
    od;
end);


# Computes (g_i-1)(g_j-1) as expression in Jennings basis
# Input: Table and maximal weight we want
# Output: None. T.commwords is filled
BindGlobal("WordFillTableComm", function(T, lvl)
local F, p, exp, pos1, i, j, pows, pos, combs, v, c, w, s, m, tup, coefprod, expprod, mm, vec, expj, pos2, poscom;

    F := T.fld;
    p := Characteristic(F);

    for i in T.pre.poswone do
        exp := T.pre.exps[i];
        pos1 := Position(exp, 1);

        # Now fill product of (g_i-1)(g_j-1)
        for j in T.pre.poswone do
            expj := T.pre.exps[j];
            pos2 := Position(expj, 1);

            # Check if weight is small enough
            if T.pre.weights[i] + T.pre.weights[j] <= lvl then 
                # Element is in the Jennings basis
                if pos2 > pos1 then 
                    T.commwords[pos1][pos2] := [[One(F), ExpToWord(exp+expj)]];
                fi;

                # Element is in the Jennings basis
                if pos1 = pos2 and p <> 2 then 
                    T.commwords[pos1][pos2] := [[One(F), ExpToWord(exp+expj)]];
                fi;

                # Element is not in the Jennings basis
                if pos2 < pos1 then 
                    poscom := Position(List(T.pre.jen.coms, x -> x[1]), [pos2,pos1]);
                    pows := StructuralCopy(T.pre.jen.coms[poscom][2]);
                    pows[pos1] := 1;
                    pows[pos2] := 1;
                    pos := PosNonzero(pows);
                    pows := pows{pos};
                    combs := Combinations([1..Length(pos)]);
                    v := [ ];

                    for c in combs do
                        # Cases for empty product or just the base elements correpsonding to (g_i-1), (g_j-1) 
	                    if not c in [ [ ], [1], [2] ] then 
                            w := [ ];
                            s := Length(c);
                            m := Maximum(pows{c});
                           
                            for tup in IteratorOfTuples([1..m], s) do
                                if tup <= pows{c} then
                                    if Sum(List([1..s], x -> tup[x]*T.pre.jen.weights[pos[c[x]]])) <= lvl then
         	                        coefprod := Product([1..s], x-> Binomial(pows[c[x]], tup[x] ));
		                        expprod := ListWithIdenticalEntries(Length(exp), 0);
		                                
                                        for mm in [1..s] do
                                            expprod[pos[c[mm]]] := tup[mm];
                                        od;
	          
                                        vec := ExpToWord(expprod);
                                	Add(w, [coefprod*One(F), vec]);
                                    fi;
                                fi;
                            od;
   
	                    Append(v, w);
                        fi;
                    od;
      
                    T.commwords[pos1][pos2] := StructuralCopy(v);
                fi;
    
            else
                T.commwords[pos1][pos2] := [ ];
            fi;
        od;
    od;
end);


# Generates a table containing Jennings data
# Input: Modular group algebra
# Output: Table containing Jennings basis data of the algebra
BindGlobal("PreSetTable", function(A)
local pre, tabs;

    pre := PrecompWeightedBasisOfRad(A);
    tabs := rec( dim := 0, # Will contain dimension of I/I^k
                 wgs := [ ], # will contain weights of elements in the Jennings bais of I/I^k 
                 rnk := pre.rnk, # rank of group basis in group algebra
                 wds := [ ], # will contain "definitions" in the sense of ModIsom
                 fld := LeftActingDomain(A), 
                 tab := [ ], # Will contain parts of structure constant table
                 pre := pre, # Jennings basis data
                 commwords := [ ], # Will contain expressions for (g_i-1)(g_j-1) for g_i, g_j pcgs of Jennings series
                 powwords := [ ]); # Will contain (g_i-1)^p

    return tabs;
end);


# Function which fills T.powwords and T.commwords. Also sets T.dim, T.wgs
# Input: Table and maximal weight we want. Typically output from PreSet
# Output: None. T.powwords, T.commwords, T.dim, T.wgs are set
BindGlobal("WordFillTable", function(T, lvl)
local F, p, s1, l, i, pos1, expw1, w1, sortweights1, sortexpw1, f1, posf1, s2, post2, expw2, w2, sortweights2, sortexpw2, f2, posf2, dep, mm, vec, posvec, fac1, posnonzero, v;

    F := T.fld;
    p := Characteristic(F);
    l := Position(T.pre.weights, lvl+1);
    
    if l = fail then
        l := p^Length(T.pre.exps[1]) - 1;
    else
        l := l -1;
    fi;

    T.dim := l;
    T.wgs := T.pre.weights{[1..l]};

    # For technical reasons so that entry is bound
    for i in [1..Length(T.pre.exps[1])] do 
        T.commwords[i] := [ ];
    od;

    WordFillTablePowers(T, lvl);
    WordFillTableComm(T, lvl);
end);


# Fill structure constant table such that products of form (g_i-1)*b are written for 
# all g_i pcgs of Jennings series and b any word in Jennings basis.
# For the other products of Jennings basis elements definitions are written
# Input: Table. Typically output of WordFillTable
# Output: None. T.tab is partly filled and T.wds are written
BindGlobal("TableOfRadByCollection", function( T )
    local n, F, N, i, j, w1, w2, len1, v, d, w, l, k;

    # set up
    n := T.dim;
    F := T.fld;
    N := ListWithIdenticalEntries(n, Zero(F));

    T.tab := [];
    
   # Loop over factors of form (g_i-1)
    for i in T.pre.poswone do 
        T.tab[i] := [];

        # Loop over the whole Jennings basis
        for j in [1..n] do 
            if T.pre.weights[i]+T.pre.weights[j] <= T.pre.weights[n] then 
                # Get word for left factor
                w1 := ExpToWord(T.pre.exps[i]); 
                # Get word for right factor
                w2 := ExpToWord(T.pre.exps[j]); 
                len1 := Length(w1);

	        for v in w2 do
                    len1 := len1 + 1;
	            w1[len1] := v;
                od;
        
                # Reduce factor word
                w1 := ReducedWord(w1); 

                # Express factor word as linear combination in Jennings basis and translate 
                # to vector of underlying field
	        T.tab[i][j] := LinComSortWordToVec(T, TransLinComToBas(T, [[ One(F), w1 ]])); 
 
            # If weight is too large, then the product will be trivial
            else
                T.tab[i][j] := N; 
            fi;
        od;
    od;

    # Generate the definitions. E.g. [0,2,1,3] = [0,1,0,0]*[0,1,1,3] is saved. 
    for i in [1..n] do 
        if not i in T.pre.poswone then
            d := PositionNonZero(T.pre.exps[i]);
            v := 0*T.pre.exps[i]; v[d] := 1; 
            l := Position(T.pre.exps, v);
            w := T.pre.exps[i]-v; 
            k := Position(T.pre.exps, w);
            T.wds[i] := [l,k];
        fi;
    od;
end);


#############################################################################
##
#O TableOfRadQuotient( A, n )
##
##
BindGlobal("TableOfRadQuotient", function( A, n )
local T;

    T := PreSetTable(A);
    WordFillTable(T, n);
    TableOfRadByCollection(T);

    return T;
end );

################################3
## new function, more user friendly
# Set up the table for calculations in I/I^(n+1)
#
## Input: finite p-group G, n so that the table for the quotient I/I^(n+1) is computed (I is the augmentation ideal) and optionally f so that it is over GF(p^f)
## Output: the table, in the sense of ModIsomExt corresponding to I/I^(n+1) of the group algebra GF(p^f)[G]
BindGlobal("ModIsomTable", function(G, L...)
local p, T, n, f, k, kG;
    p := PrimePGroup(G);
    n := L[1];
    if Length(L) = 2 then
        f := L[2];
    else
        f := 1;
    fi;
    k := GF(p^f);
    kG := GroupRing(k, G);
    T := PreSetTable(kG);
    WordFillTable(T, n);
    TableOfRadByCollection(T);
    return T;
end);
