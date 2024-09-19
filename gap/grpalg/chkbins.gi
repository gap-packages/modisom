#############################################################################
##
#W    checkbin.gi               From the ModIsomExt package          Leo Margolis
#W                                                                   Tobias Moede
##  small changes from ModIsomExt, now possible to work over bigger fields than just GF(p)


#############################################################################
##
#F RefineBins( bins, vals )
##
## These functions reamin unchanged from earlier versions
##
BindGlobal( "RefineBinByVals", function( bin, vals )
    return List(Set(vals{bin}), x -> Filtered(bin, y -> vals[y]=x));
end );
    
BindGlobal( "RefineBinsByVals", function( bins, vals )
    local refs;
    refs := Concatenation(List(bins, x -> RefineBinByVals(x, vals )));
    return Filtered(refs, x -> Length(x)>1);
end );

#############################################################################
##
#O MIPBinSplit( p, n, k, start, step, list )
##
##
BindGlobal( "MIPBinSplit", function(p, n, k, start, step, L...)
local grps, algs, tabs, bins, l, m, vals, todo, i, j, d, splits, stime, f, list;

    list := L[1];
    # the size of the field we will work over is p^f. If no size is specified, we work over GF(p)
    if Length(L) = 2 then
        f := L[2];
    else
        f := 1;
    fi;

    if not IsPrimeInt(p) then
        Error("<p> must be a prime integer");
    fi;

    if not (IsPosInt(k) or IsBool(k)) then
        Error("<k> must be a positive integer or a boolean value");
    fi;

    stime := Runtime();
    l := Length(list);
 
    # Check whether we have a list of integers or a list of groups
    if ForAll(list, IsInt) then
        grps := List(list, x -> SmallGroup(p^n, x));
    elif ForAll(list, IsGroup) then
        grps := list;
    else
        Error("<L> must be a list of integers or a list of groups");
    fi;

    algs := List(grps, x -> GroupRing(GF(p^f), x));  # changed here with respect to ModIsomExt
    tabs := [];

    for i in [1..l] do
        tabs[i] := TableOfRadQuotient(algs[i], start);
    od;

    # Init bins
    bins := [[1..l]];
    splits := [];

    Info(InfoModIsom, 1, "Refine bin");

    # Refine by weights
    vals := List( tabs, x -> x.pre.weights );
    bins := RefineBinsByVals( bins, vals );
    Info(InfoModIsom, 1, "  Weights yields bins ", bins);
    todo := Flat(bins);

    if Length(todo) = 0 then 
        return rec( bins := [], splits := [[0, list]], time := Runtime()-stime ); 
    fi;

    # Start up and refine by aut group
    vals := [];
    for i in [1..l] do
        if i in todo then 
            InitCanoForm(tabs[i]); 
            vals[i] := [tabs[i].auto.size, tabs[i].auto.partition];
        fi;
    od;
    bins := RefineBinsByVals( bins, vals );
    Info(InfoModIsom, 1, "  Layer 1 yields bins ", bins);
    todo := Flat(bins);

    if Length(todo) = 0 then 
        return rec( bins := [], splits := splits, time := Runtime()-stime ); 
    fi;

    # Get limit
    m := Maximum( List(tabs, x -> Maximum(x.pre.weights)));
    if not IsBool(k) then m := Minimum(m,k); fi;

    # Loop and refine by canonical forms
    for j in [2..m] do
        if j > start then
	        start := start+step;
	        for i in todo do
	            WordFillTable(tabs[i], start);
	    	    TableOfRadByCollection(tabs[i]);
	        od;
	    fi;	

        vals := [ ];
    
        for i in [1..l] do
            if i in todo then 
                ExtendCanoForm(tabs[i], j);
				# Either use CompleteTable or do as in CompareTables 
		CompleteTable(tabs[i].cano);
                vals[i] := tabs[i].cano.tab;  # this is the only (?) change compared to original function          
            else
	        tabs[i] := "done"; 
	    fi;
        od;
    
        bins := RefineBinsByVals( bins, vals );
        Info(InfoModIsom, 1, "  Layer ", j, " yields bins ", bins);
        
        # Record split information  
        if Length(Flat(bins)) < Length(todo) then    
            Add(splits, [j, list{Difference([1..l], Flat(bins))}]);  
        fi;
        todo := Flat(bins);
    
        if Length(todo) = 0 then 
            return rec( bins := [], splits := splits, time := Runtime()-stime );
        fi;
    od;

    return rec( bins := List(bins, x -> List(x, y -> list[y])), splits := splits, time := Runtime()-stime );
end );

#### more user friendly function.
## Input: list of groups as first argument. Optionally integer f so that the algebras are taken over GF(p^f), otherwise GF(p) is used
BindGlobal("MIPSplitGroupsByAlgebras", function(L...)
local f, list, G, p, n, start, step, listgroups;
   
    if IsInt(L[1]) then
        p := L[1];
        n := L[2];
        if Length(L) = 4 then
            f := L[4];
        else
            f := 1;
        fi;
        list := L[3];
        listgroups := List(list, x -> SmallGroup(p^n, x));
        start := JenningsBoundPairwise(p, n, list); # this will be the smallest value which can potentially separate all the algebras
        step := 1;
        return MIPBinSplit(p, n, false, start, step, list, f);
    elif IsList(L[1]) then  
        if Length(L) = 2 then
            f := L[2];
        else
            f := 1;
        fi;
        list := L[1];
        G := list[1];
        p := PrimePGroup(G);
        n := Log(Size(G), p);
        if IdGroupsAvailable(p^n) then
            start := JenningsBoundPairwise(p, n, list); # this will be the smallest value which can potentially separate all the algebras. But is slow in general
        else
            start := 2;
        fi;
        step := 1;
        return MIPBinSplit(p, n, false, start, step, list, f);
    fi;
end );
