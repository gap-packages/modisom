#####################################################################################
# Replaces the function RefineBins from ModIsom
## Input: (order of groups, list of lists of groups, function to apply, boolean on how much info is printed)
# output: bins are sorted finer by applying func
#InstallGlobalFunction(RefineBins,  function( o, bins, func, prnt ) # this was for a AutDoc experiment
BindGlobal("RefineBins", function( o, bins, func, prnt )
    local i, j, k, res, inv, new;

    for i in [1..Length(bins)] do

        # report
        if prnt then
            Print("  start bin ",i," of ",Length(bins),"\n");
        fi;

        # compute values
        res := [];
        for j in [1..Length(bins[i])] do
            if IsInt(bins[i][j]) then
                res[j] := func(SmallGroup(o, bins[i][j]));
            else
                res[j] := func(bins[i][j]);
            fi;
        od;
        inv := Set(res);

        # split bin
        new := List(inv, x -> []);
        for j in [1..Length(res)] do
            k := Position(inv,res[j]);
            Add( new[k], bins[i][j] );
        od;

        # reset bin
        bins[i] := new;
    od;

    # return result
    if IsInt(bins[1][1]) then
        return SortedList(Concatenation(bins));
    else
        return Concatenation(bins);
    fi;
end) ;


## return GroupId if this is available for this order and abelianization otherwise
## not changed from previous versions
BindGlobal( "GroupInfo", function(G)
    if not IsBool(ID_AVAILABLE(Size(G))) then 
        return IdGroup(G);
    else
        return [Size(G), AbelianInvariants(G)];
    fi;
end );


#############################################################################
##
#F Group Theoretical Invariants. 
## All functions to follow have as input a finite p-group 
##
## information from conjugacy classes: powers and dimension of first Hochschild cohomology
BindGlobal("ConjugacyClassInfo", function(G)
    local p, e, c, s, o, exp, i, cov, r, cc, g, h;

    p := PrimePGroup(G);
    e := Elements(G);
    c := Orbits(G,e);

    # Roggenkamp (actually not him, just dimension of first Hochschild cohomology group)
    s := [Sum(List(c, x -> RankPGroup(Stabilizer(G,x[1]))))];

    # the length
    Add( s, Length(c) );

    # p^n-th powers. Kuelshammer
    repeat
        e := Set( List(e, x -> x^p) );
        o := Orbits(G, e);
        Add( s, Length(o) );
    until Length(e) = 1;

    # Parmenter+Polcino-Milies. Hertweck+Soriano Lemma 2.4. 
    # This can probably be written much more efficiently!
    exp := Log(Exponent(G), p);
    c := Filtered(c, x -> Size(x) > 1); # Filter out center to improve performance
    for i in [1..exp-1] do
        cov := [ ];
        r := 0;
        
        for cc in c do
            g := Representative(cc);
        
            if Order(g) < p^i then # Serves to improve performance. Does it?
                h := g^(p^i);
                if not h^G in cov then
                    if Size(g^G) = Size( h^G ) then
                        r := r + 1;
                        Add(cov, h^G);
                    fi;
                fi;
            fi;
        od;

        Add(s,r);
    od;

    return s;
end);

### information from quotients in the Jennings series aka dimension subgroup series
BindGlobal("JenningsInfo", function(G)
    local s, r, i, a;
    s := JenningsSeries(G);
    r := [ ];
    for i in [1..Length(s)-1] do
        a := [GroupInfo(s[i]/s[i+1])];
        if i <= Length(s)-2 then
            Add(a, GroupInfo(s[i]/s[i+2]));
        fi;
        if i <= (Length(s)-1)/2 then
            Add(a, GroupInfo(s[i]/s[2*i+1]));
        fi;
        Add(r, a);
    od;

    a := [ ];

    # Hertweck's G/M4(G) invariant
    if PrimePGroup(G) <> 2 then 
        # Not just use group info as it is weaker for orders 5^6 and 3^7 (ID_AVAILABLE fails for these orders)
        if Length(s) <= 4 and "IdGroup" in KnownAttributesOfObject(G) then
            Add(a, IdGroup(G)); 
        else
            Add(a, GroupInfo(G/s[4]));
        fi;
    # Not Hertweck, but can help for 2-groups of order 2^9
    else 
        # Not just use group info as it is weaker for orders 2^9 (ID_AVAILABLE fails for these orders)
        if Length(s) <= 3 and "IdGroup" in KnownAttributesOfObject(G) then
            Add(a, IdGroup(G)); 
        else
            Add(a, GroupInfo(G/s[3]));
        fi;
    fi;
    Add(r,a);
    return r;
end);

# for all fields only subsequent quotients are known invariants
BindGlobal("JenningsInfoAllFields", function(G)
    local s, r, i, a;
    s := JenningsSeries(G);
    r := [ ];
    for i in [1..Length(s)-1] do
        a := [GroupInfo(s[i]/s[i+1])];
        Add(r, a);
    od;
    return r;
end);

## Hertweck-Soriano end of page 16
## As this is based on the small group algebra it is included in SandlingInfo below
BindGlobal("HertweckSorianoFrattiniInfo", function(G) 
local p;
  
    p := PrimePGroup(G);
    
    if IsElementaryAbelian(DerivedSubgroup(G)) then
        if Length(JenningsSeries(G)) <= (2*p) then
            return GroupInfo(FrattiniSubgroup(G));
        fi;
    fi;
  
    return false;
end);


## infromation from small group algebra as first done by Sandling
BindGlobal("SandlingInfo", function(G) 
    local s, p, U, r, d, W, V;
    if IsAbelian(G) then # otherwise lower central series has length at least 3 and the following parts work
        return GroupInfo(G);
    fi;
    s := LowerCentralSeries(G);
    p := PrimePGroup(G);
    r := [ ];
    U := Subgroup(G, Concatenation(Pcgs(s[3]), List(Pcgs(s[2]), x -> x^p)));
    if Size(U) = 1 and "IdGroup" in KnownAttributesOfObject(G) then
        Add(r,IdGroup(G));
    else
        Add(r, GroupInfo(G/U));
    fi;

    # Baginski 99 / Margolis+Moede 2022
    if RankPGroup(G) <= 2 and Length(s) >= 4 then 
        V := Subgroup(G, Concatenation(Pcgs(s[4]), List(Pcgs(s[2]), x -> x^p)));
        if Size(V) = 1 and "IdGroup" in KnownAttributesOfObject(G) then
            Add(r,IdGroup(G));
        else
            Add(r, GroupInfo(G/V));
        fi;
    fi;
    # from Hertweck+Soriano 06, p.16
    Add(r, HertweckSorianoFrattiniInfo(G));

    # Brenner+Garcia-Lucas 24, Theorem B
    if p > 2 and Length(s) >= 4 then
        W := Subgroup(G, Concatenation(Pcgs(FrattiniSubgroup(G)), Pcgs(Center(G))) ); #Phi(G)Z(G)
        d := Log(Size(G)/Size(W), p);
        if IsSubgroup(U, Intersection(s[2], Agemo(G, p))) and Size(s[2]/U) = p^Binomial(d, 2) then
            V := Subgroup(G, Concatenation(Pcgs(s[4]), List(Pcgs(s[2]), x -> x^p)));
            if Size(V) = 1 and "IdGroup" in KnownAttributesOfObject(G) then
                Add(r, IdGroup(G));
            else
                Add(r, GroupInfo(G/V));
            fi;
        fi;
    fi;
    return(r);
end);

## function unchanged from earlier version
## calculates number of conjugacy classes of maximal elementary abelian subgroups all all ranks
BindGlobal( "SubgroupsInfo", function(G)
    local lat, cls, max, sub, new, U, p, l, d, i, j, k;

    # all elementary-abelian subgroups
    lat := LatticeByCyclicExtension( G, IsElementaryAbelian );
    cls := ConjugacyClassesSubgroups( lat );

    # get maximal size
    p := PrimePGroup(G);
    l := LogInt(Maximum(List(cls, x->Size(Representative(x)))),p);

    # set up
    max := Filtered( cls, x -> Size(Representative(x)) = p^l );

    for i in Reversed( [1..l-1] ) do
        sub := Filtered( cls, x -> Size(Representative(x)) = p^i );
        new := [];
        for j in [1..Length(sub)] do
            U := Representative(sub[j]);
            k := 1;
            repeat
                d := ForAny( Elements(max[k]), x -> IsSubgroup(x,U) );
                k := k + 1;
            until d = true or k > Length(max);
            if d = false then Add( new, sub[j] ); fi;
        od;
        Append( max, new );
    od;

    d := List(max, x -> RankPGroup(Representative(x)));
    return List([1..Maximum(d)], x -> Length(Filtered(d, y -> y=x)));
end );

###### Some new invariants (new compared to ModIsom versions 1 and 2). First an article of Baginski


## Baginski 99 Corollary 7 and Theorem 9, also Baginski+Zabielski 25 Lemma 7
BindGlobal("BaginskiInfo", function(G) 
local D, N, act, F, r, p;
    D := DerivedSubgroup(G);
    F := FrattiniSubgroup(D);
    p := PrimePGroup(G);
 
    # if-condition only for technical reasons, output is the same
    if Size(D) = Size(F) then
        N := Centralizer(G, D);
    else
        act := AbelianSubfactorAction(G, D, FrattiniSubgroup(D));
        N := Kernel(act[1]);
    fi;
    
    if IsCyclic(G/N) then
        return [GroupInfo(G/FrattiniSubgroup(N)), GroupInfo(N/F)];
    elif Order(G/N) = p then
        return [GroupInfo(G/FrattiniSubgroup(N)), GroupInfo(N/F), JenningsInfoAllFields(N)];
    else
        return false;
    fi;
end);

## Based on Baginski+Caranti 1988 Proposition 1.2
BindGlobal("BaginskiCarantiInfo", function(G) 
    return NilpotencyClassOfGroup(G/FrattiniSubgroup(DerivedSubgroup(G)));
end);

## Based on Theorem 6.11 from Sandling's 1984 survey
BindGlobal("CenterDerivedInfo", function(G) 
local U;
    U := Intersection(Center(G), DerivedSubgroup(G));
   
    return [AbelianInvariants(U), AbelianInvariants(Center(G)/U)];
end);


## Based on Sandling survey Lemma 6.26 for G'
BindGlobal("JenningsDerivedInfo", function(G) 
local s, r, a, i;
  
    s := JenningsSeries(DerivedSubgroup(G));
    r := [ ];

    for i in [1..Length(s)-1] do
        Add(r, GroupInfo(s[i]/s[i+1]));
    od;
  
    return r;
end);


### cases where the nilpotency class is known to be an invariant
# Baginski+Konovalov 05
BindGlobal("NilpotencyClassInfo", function(G) 
local p;
  
    p := PrimePGroup(G);

    # Note that all of these three are invariants anyway by other results
    if Exponent(G) = p or NilpotencyClassOfGroup(G) = 2 or IsCyclic(DerivedSubgroup(G)) then 
        return NilpotencyClassOfGroup(G);
    else
        return false;
    fi;
end);


#### A bit of cohomology. Dimension of first cohomology in included in first step of Jennings series. Dimension of first Hochschild cohoology inlcuded in the conjugacy classes calculations. Here we compute dimension of second cohomology for the trivial module and the dimension of second Hochschild cohomology
BindGlobal("DimensionTwoCohomology", function(G)
local p, M, L, C;

    p := PrimePGroup(G);
    L := List([1..Size(GeneratorsOfGroup(G))], x -> [[Z(p)^0]]);
    M := GModuleByMats(L, GF(p));
    C := TwoCohomology(G,M);

    return Dimension(Image(C.cohom));
end);

BindGlobal("DimensionSecondHochschild", function(G)
local CC, R, Cs, s, C;

    CC := ConjugacyClasses(G);;
    R := List(CC, Representative);;
    Cs := List(R, x -> Centralizer(G, x));;
    s := 0;
    for C in Cs do
        s := s + DimensionTwoCohomology(C);
    od;
    return s;
end);




####################################
# Programs to test if group is covered by theoretical result. 
# Not including results for small groups
BindGlobal("HasMaximalAbelianSubgroup", function(G)
local LM, U;

    LM := MaximalSubgroups(G);
  
    for U in LM do
        if IsAbelian(U) then
            return true;
        fi;
    od;
    
    return false;
end);
 
BindGlobal("HasCyclicSubgroupIOfIndexP2", function(G) 
local LM, U, LUM, V;

    LM := MaximalSubgroups(G);

    for U in LM do
        LUM := MaximalSubgroups(U);
    
        for V in LUM do
            if IsCyclic(V) then 
                return true;
            fi;
        od;
    od; 
    
    return false;
end);

# Auxilary function to apply in next function
# to compute the criterion in Margolis+Stanojkovski 22, Theorem 3.5.
BindGlobal("Theorem35MS22Applies", function(G)
local p, F, UPS, Z2, Gam, LM, U, KG;

    if IsAbelian(G) then
        return GroupInfo(G);
    fi;
    p := PrimePGroup(G);
    F := FrattiniSubgroup(G);
    UPS := UpperCentralSeries(G);
    Z2 := UPS[Size(UPS)-2];
    if Index(G, F) <> p^3 or F <> Z2 then
        return false;
    else
        Gam := Intersection(DerivedSubgroup(G), Center(G)); 
        if Index(DerivedSubgroup(G), Gam) = p^3 then
            KG := Z2;
        else
            LM := MaximalSubgroups(G);
            for U in LM do
                if IsSubgroup(Gam, DerivedSubgroup(U)) then
                    KG := U;
                    break;
                fi;
            od;
        fi;
        if IsSubset(CommutatorSubgroup(KG, DerivedSubgroup(G)), CommutatorSubgroup(CommutatorSubgroup(KG,G),G)) then
            return true;
        else
            return false;
        fi;    
    fi;
end);

# to see if group is metacylic
BindGlobal("IsMetacyclicGroup", function(G)
    local LN, N;
    if IsCyclic(DerivedSubgroup(G)) then
        LN := NormalSubgroups(G);
        for N in LN do 
            if IsCyclic(N) and IsCyclic(G/N) then
                return true;
            fi;
        od;
    fi;
    return false;
end);

### 
## to apply Roehl 90, Theorem 3.1.2
# first the dimensions of Jennings quotients of d-generated free pro-p group. Based on Minac+Rogelstad+Tan 2016, Proposition 3.4
BindGlobal("MinacRogelstadTanWnFunction", function(d, n)
    local res, m;
    res := 0;
    for m in DivisorsInt(n) do
        res := res + (1/n)*MoebiusMu(m)*d^(n/m);
    od;
    return res;
end);

BindGlobal("MinacRogelstadTanCnFunction", function(d, n, p)
    local facs, k, m, res, i;
    facs := FactorsInt(n);
    k := Size(Positions(facs, p));
    m := n/(p^k);
    res := 0;
    for i in [0..k] do
        res := res + MinacRogelstadTanWnFunction(d, (p^i)*m);
    od;
    return res;
end);

BindGlobal("SatisfiesRoehl90", function(G)
    local p, JSS, d, i;
    p := PrimePGroup(G);
    JSS := List(JenningsSeries(G), Size);
    d := Log(Size(G/FrattiniSubgroup(G)), p);
    if Size(JSS) = 1 then # elementary abelian case 
        return true;
    fi;
    for i in [1..Size(JSS)-2] do
        if Log(JSS[i]/JSS[i+1], p) <> MinacRogelstadTanCnFunction(d, i, p) then
            return false;
       fi;
    od; 
    return true;
end);


#### Some theoretical results follow from the other criteria. In particular metacyclic groups (Sandling 96, for prime field) and (elem-ab.)-by-cyclic groups (Baginski 99) are covered.

BindGlobal("IsCoveredByTheory", function(G)
local p, n, D, F, N, act;

    p := PrimePGroup(G);
    n := Log(Size(G),p);

    # 2-groups of maximal class known. Baginski
    if p = 2 and NilpotencyClassOfGroup(G) = n-1 then 
        return true;
    fi;

    # Baginski+Caranti 88
    if p <> 2 and n <= p+1 and HasMaximalAbelianSubgroup(G) and NilpotencyClassOfGroup(G) = n-1 then 
        return true;
    fi;

    # Drensky 89
    if Size(G)/Size(Center(G)) = p^2 then 
        return true;
    fi;

    # Brenner+Garcia-Lucas 24
    if p > 2 and Size(G)/Size(Center(G)) = p^3 then 
        return true;
    fi;
   
 
    # Baginski+Konovalov 05. Rather slow, but can be helpful. Case p<>2 covered by other criteria by Lemma 1 (ArXiv-version) of BK05
    #if p = 2 and HasCyclicSubgroupIOfIndexP2(G) then 
        #return true;
    #fi;

    # as Margolis+Sakurai25, proof of Theorem 3.9, shows the following is covered by other invariants (in particular NormalSubgroupsInfo)
    # Broche+del Rio 20, Theorem 1. The case of odd primes is covered by other invariants
 #   if NilpotencyClassOfGroup(G) = 2 and Size(G/FrattiniSubgroup(G)) = p^2 then 
 #       return true;
 #   fi;


    # Margolis+Stanojkovski 22: Theorem 3.3, Theorem 3.5 
    if p <> 2 and NilpotencyClassOfGroup(G) <= 3 and Exponent(DerivedSubgroup(G)) = p then # note: these conditions imply the class is exactly 3
        if Size(Centralizer(G, DerivedSubgroup(G))) = Size(G)/p and IsAbelian(Centralizer(G, DerivedSubgroup(G))) then 
            return true;
        elif Size(G/FrattiniSubgroup(G)) = p^3 and Size(G)/Size( UpperCentralSeries(G)[Size(UpperCentralSeries(G))-2] ) = p^3 and Theorem35MS22Applies(G) then
            return true;
        fi;
    fi;

    # Margolis+Sakurai+Stanojkovski23, Theorem 5.6
    if p = 2 and IsCyclic(Center(G)) and IsDihedralGroup(G/Center(G)) then
        return true;
    fi;

    # Garcia-Lucas+Margolis24. Odd p covered by other invariants
    if p = 2 and IsCyclic(Center(G)) and NilpotencyClassOfGroup(G) = 2 then
        return true;
    fi;

    # Garcia-Lucas+del Rio 24, Proposition 3.7
    if p <> 2 and IsCyclic(DerivedSubgroup(G)) and IsCyclic( Agemo(G/DerivedSubgroup(G), p, 2) ) then
        return true;
    fi;

    # Baginski+Zabielski 25, Proposition 8
    # if-condition only for technical reasons, output is the same
    if Size(G/FrattiniSubgroup(G)) = p^2 then
        D := DerivedSubgroup(G);
        F := FrattiniSubgroup(D);
        if Size(D) = Size(F) then
            N := Centralizer(G, D);
        else
            act := AbelianSubfactorAction(G, D, FrattiniSubgroup(D));
            N := Kernel(act[1]);
        fi;
        if IsAbelian(N) and Size(G/N) = p then
            return true;
        fi;
    fi; 

    #Roehl 90, Proposition 3.2
    if SatisfiesRoehl90(G) then
        return true;
    fi;

    return false;
end);


# variation of previous function only where MIP is solved for all fields
BindGlobal("IsCoveredByTheoryAllFields", function(G)
local p, n, q;

    p := PrimePGroup(G);
    n := Log(Size(G),p);

    # 2-groups of maximal class known. Baginski
    if p = 2 and NilpotencyClassOfGroup(G) = n-1 then 
        return true;
    fi;

    # Drensky 89
    if Size(G)/Size(Center(G)) = p^2 then 
        return true;
    fi;
 
    # Margolis+Sakurai+Stanojkovski 23, Theorem 5.6
    if p = 2 and IsCyclic(Center(G)) and IsDihedralGroup(G/Center(G)) then
        return true;
    fi;

    # Garcia-Lucas+Margolis 24, case of all fields for p=2
    if p = 2 and IsCyclic(Center(G)) and NilpotencyClassOfGroup(G) = 2 then
        q := AbelianInvariants(G/Center(G));
        if Size(q) <= 2 or q[Size(q)] > q[Size(q)-2] then
            return true;
        fi;
    fi;

    # as Margolis-Sakurai25 shows the following is covered by other invariants
    # Margolis-Sakurai25, Theorem 3.9. The case of odd primes is covered by other invariants
#    if p = 2 and NilpotencyClassOfGroup(G) = 2 and Size(G/FrattiniSubgroup(G)) = 4 then 
#        return true;
#    fi;

   # Margolis+Sakurai 25, Theorem 3.6. For prime field this is covered by other invariants, for all fields those used might not be enough if GroupInfo is not strong enough
   if IsMetacyclicGroup(G) then
       return true; 
   fi;  


    return false;
end);


### Theorem 4.1 from Margolis+Stanojkovski 22
BindGlobal("Theorem41MS22", function(G)
local p, LCS;

    p := PrimePGroup(G);
    if p <> 2 and Size(G)/Size(FrattiniSubgroup(G)) > p^2 and NilpotencyClassOfGroup(G) = 3 and Exponent(LowerCentralSeries(G)[3]) = p then 
        LCS := LowerCentralSeries(G);
        return [AbelianInvariants(LCS[2]), AbelianInvariants(LCS[3])]; # subgroups are certainly abelian
    else
        return false;
    fi;
end);

####
# functions to recognize maximal abelian or elemantary abelian direct factors
BindGlobal("MaximalAbelianDirectFactor", function(G)
local factors, H,e,i,p,F;

    e := Exponent(G/DerivedSubgroup(G));
    p := PrimePGroup(G);
    factors := [ ];
    for i in [1..Log(e,p)] do
        H := Subgroup(G, Union( [GeneratorsOfGroup(Omega(Center(G),p,i)), GeneratorsOfGroup(Agemo(G,p, i)), GeneratorsOfGroup(DerivedSubgroup(G))] )) /Subgroup(G, Union(GeneratorsOfGroup(Agemo(G, p,i)), GeneratorsOfGroup(DerivedSubgroup(G))));
        F := Filtered(AbelianInvariants(H), x-> x=p^i);
        factors := Concatenation(factors, F);
    od;
    return factors;
end);

BindGlobal("MaximalElementaryAbelianDirectFactor", function(G)
local factors, H,e,i,p,F;

    e := Exponent(G/DerivedSubgroup(G));
    p := PrimePGroup(G);
    factors := [ ];
    H := Subgroup(G, Union( [GeneratorsOfGroup(Omega(Center(G),p)), GeneratorsOfGroup(Agemo(G,p)), GeneratorsOfGroup(DerivedSubgroup(G))] )) /Subgroup(G, Union(GeneratorsOfGroup(Agemo(G, p)), GeneratorsOfGroup(DerivedSubgroup(G))));
    F := Filtered(AbelianInvariants(H), x-> x=p);
    factors := Concatenation(factors, F);
    return factors;
end);


#### Invariants from Theorem B of Margolis+Sakurai+Stanojkovski 23, Theorem B, and Garcia-Lucas22, Example 3.7(2)
# The Agemo*-invariant
BindGlobal("ZAstar", function(G, m)
local p, D, A, Z, W, U;

    p := PrimePGroup(G);
    D := DerivedSubgroup(G);
    A := Agemo(G, p, m);
    Z := Center(G);
    W := ClosureGroup(D, A);
    U := Intersection(W, Z);
    return [AbelianInvariants(U), AbelianInvariants(Z/U)]; # group is abelian anyhow, so this is full info
end);

# the data for all values at the same time
BindGlobal("AgemoInvariantAllM", function(G)
local res, p, l, m;

    res := [ ];
    p := PrimePGroup(G);
    l := Log(Exponent(G), p);
    for m in [1..l-1] do
        Add(res, ZAstar(G, m)); 
    od;
    return res;
end);


### The Omega*-invariant
BindGlobal("GammaOstar", function(G, m)
local p, D, A, W;

    p := PrimePGroup(G);
    D := DerivedSubgroup(G);
    A := Omega(Center(G), p, m);
    W := ClosureGroup(D, A);
    return [AbelianInvariants(G/W), AbelianInvariants(W/D)]; # group is abelian anyhow, so this is full info
end);

# the data for all values at the same time
BindGlobal("OmegaInvariantAllM", function(G)
local res, p, l, m;

    res := [ ];
    p := PrimePGroup(G);
    l := Log(Exponent(Center(G)), p); 
    for m in [1..l-1] do
        Add(res, GammaOstar(G, m)); 
    od;
    return res;
end);

# Garcia-Lucas Example 3.7(2) invariant: G/G'Agemo_m(Center(G)) and its dual 
BindGlobal("ModuloAgemoCenter", function(G, m)
local p, D, A, W;

    p := PrimePGroup(G);
    D := DerivedSubgroup(G);
    A := Agemo(Center(G), p, m);
    W := ClosureGroup(D, A);
    return [AbelianInvariants(G/W), AbelianInvariants(W/D)]; # group is abelian anyhow, so this is full info
end);

# data for all values of m
BindGlobal("AgemoCenterInvariantAllM", function(G)
local res, p, l, m;

    res := [ ];
    p := PrimePGroup(G);
    l := Log(Exponent(Center(G)), p);  
    for m in [1..l-1] do
        Add(res, ModuloAgemoCenter(G, m)); 
    od;
    return res;
end);


#### invariants for cyclic derived subgroup and p odd from Garcia-Lucas+del Rio+Stanojkovski 22, Garcia-Lucas+del Rio 24 and Margolis+Sakurai 25

# auxilary program for Corollary E of GLdRS22

BindGlobal("NumberSmallerValues", function(L, m)
local i;

    for i in [1..Size(L)] do
        if L[i] < m then
            return i-1;
        fi;
    od;
    return Size(L);
end);

# appearing in Corollary E of GLdRS22
BindGlobal("TypeInvariants", function(G)
local p, OS, omegas, res, i;
  
    p := PrimePGroup(G);
    OS := List([1..Log(Exponent(G), p)], x -> Agemo(G, p, x));
    omegas := [Size(G)/Size(OS[1])];
    for i in [2..Size(OS)] do
        Add(omegas, Size(OS[i-1])/Size(OS[i]));
    od;
 
    res := [ ];
    for i in [1..omegas[1]] do
        Add(res, NumberSmallerValues(omegas, i));
    od;
    return res;
end);

BindGlobal("CyclicDerivedInfo", function(G)
local p, D, C, JS, a, i, LCS, res1, res2, U;

    if IsAbelian(G) then
       return GroupInfo(G);
    fi;  
    p := PrimePGroup(G);
    D := DerivedSubgroup(G);
    res1 := [ ];
    if IsCyclic(D) then
        if p <> 2 then
            C := Centralizer(G, D);
            JS := JenningsSeries(C);
            a := [ ];
            for i in [1..Size(JS)-1] do
                Add(a, AbelianInvariants(JS[i]/JS[i+1]));
            od; 
            LCS := LowerCentralSeries(G);
            if Length(LCS) = 3 then # just to avoid an error in the following
                U := TrivialSubgroup(G);
            else
                U := Agemo(LCS[3], p);
            fi;
            res1 := [a, Exponent(C), AbelianInvariants(C/D), AbelianInvariants(C), GroupInfo( G/U ), GroupInfo( G/Agemo(D,p,3) ) ]; # from Theorem 4.2 in GLdRS22 and Theorem A in GLdR
            if Size(G/FrattiniSubgroup(G)) = p^2 then
                res2 := [GroupInfo(C), TypeInvariants(G)]; # Corollaries D and E from GLdRS22
                res1 := Concatenation(res1, res2);
            fi;
        fi;
        if RankPGroup(G) = 2 then # Sandling/Margolis+Sakurai
            Add(res1, GroupInfo(G/Agemo(D, p)));
        fi;
    else
        res1 := false;
    fi;
    return res1;
end);

# the results in GLdRS apply to all fields, but those in GLdR do not, so they are removed in the following function
BindGlobal("CyclicDerivedInfoAllFields", function(G)
local p, D, C, JS, a, i, LCS, res1, res2;
 
    if IsAbelian(G) then
       return GroupInfo(G);
    fi;  
    p := PrimePGroup(G);
    D := DerivedSubgroup(G);
    res1 := [ ];
    if IsCyclic(D) then
        if p <> 2 then
            C := Centralizer(G, D);
            JS := JenningsSeries(C);
            a := [ ];
            for i in [1..Size(JS)-1] do
                Add(a, AbelianInvariants(JS[i]/JS[i+1]));
            od; 
            LCS := LowerCentralSeries(G);
            res1 := [a, Exponent(C), AbelianInvariants(C/D), AbelianInvariants(C) ]; # from Theorem 4.2 in GLdRS22 and Theorem A in GLdR
            if Size(G/FrattiniSubgroup(G)) = p^2 then
                res2 := [GroupInfo(C), TypeInvariants(G)]; # Corollaries D and E from GLdRS22
                res1 := Concatenation(res1, res2);
            fi;
        fi;
        if RankPGroup(G) = 2 then # Margolis-Sakurai 25
            Add(res1, GroupInfo(G/Agemo(D, p)));
        fi;
    else
        res1 := false;
    fi;
    return res1;
end);


#############3 the following function implement part of the subgroup lattice which is canonical as described in Garcia-Lucas 24
### written (mostly) by Diego Garcia-Lucas
##The next three functions concern the info described in Lemma 3.2
BindGlobal("OnlyJenningsInfo", function(G)
  return List(JenningsSeries(G), Order);
end);

BindGlobal("CenterNormalInfo", function(G, N)
  local ZI, ZQ;
  ZI := Intersection(Center(G), N);
  ZQ := Subgroup(G, Union(GeneratorsOfGroup(Center(G)), GeneratorsOfGroup(N)))/N;
  return [ AbelianInvariants(ZI), AbelianInvariants(ZQ)];
end);


#L is normal, N contains G'
BindGlobal("LNInfo", function(G, L, N)
  local LN, Q1, Q2;
  LN := Subgroup(G, Union(GeneratorsOfGroup(N), GeneratorsOfGroup(L)));
  Q1 := G/LN;
  Q2 := LN/N;
  return [AbelianInvariants(Q1), AbelianInvariants(Q2)];
end);


# The next three functions are the operations preserving canonicity described in Lemma 3.6.
BindGlobal("OmegaModuloN", function(G, N, p, t)
  local f, Q, R;
  if Size(N)=Size(G) then return G; else
  f := NaturalHomomorphismByNormalSubgroup( G, N );
  Q := Image(f);
  R := Omega(Q, p, t);
  return PreImage(f, R);
  fi;
end);

BindGlobal("AgemoLModuloN", function(G, L, N, p, t)
  local AgemoL;
  AgemoL := Agemo(L, p, t);
  return Subgroup(G, Union(GeneratorsOfGroup(AgemoL), GeneratorsOfGroup(N)));
end);


BindGlobal("OmegaCenterN", function(G, N, p, t)
    local OmegaCenter;
    OmegaCenter := Omega(Center(G), p, t);
    return Subgroup(G, Union(GeneratorsOfGroup(OmegaCenter), GeneratorsOfGroup(N)));
end);

#N is a subgroup containing G'. Then the following function computes the subgroups obtained using the previous operations.
#(Observe that if L=G, t doesn't really matter since the output for different values of t can be obtained iterating of the operation).
BindGlobal("SuccessorsN", function(G, L, N)
  local S, SizeN, successorsN,trivialSubgroup, p;
  p := PrimePGroup(G);
  trivialSubgroup := Group(One(G));
  successorsN := [ ];
  SizeN := Size(N);
  if SizeN = 1 then
    return [ ];
  else
    S := OmegaModuloN(G, N, p, 1);
    if Size(S) <> SizeN then 
      Add(successorsN, ["Omega", S]); 
    fi;
    S := AgemoLModuloN(G, L, N, p, 1);
    if Size(S) <> SizeN then 
      Add(successorsN, ["Agemo", S]);  
    fi;
    S := OmegaCenterN(G, N, p, 1);
    if Size(S) <> SizeN then 
      Add(successorsN, ["OmegaCenter", S]);
    fi;
    return successorsN;
    fi;
end);

#!
BindGlobal("CanonicalNormalSubgroups", function(G)
  local normalSubgroups, normalSubgroupsLeveli, normalSubgroupsLevelip1, e, Gprime, Ab, SN, N, p,  i, SS, SString;
  Gprime := DerivedSubgroup(G);
  normalSubgroups := [ [[ ], Gprime] ];
  p := PrimePGroup(G);
  Ab := G/Gprime;
  e := 3*Log(Exponent(Ab), p);
  normalSubgroupsLeveli := normalSubgroups;
  for i in [1..e] do
    normalSubgroupsLevelip1 := [ ];
    for N in normalSubgroupsLeveli do
        SN := SuccessorsN(G, G, N[2]);
        for SS in SN do
            SString := Concatenation(N[1], [SS[1]]);
            Add(normalSubgroupsLevelip1, [SString, SS[2]]);
        od;
    od;
    normalSubgroupsLeveli := normalSubgroupsLevelip1;
    normalSubgroups := Concatenation(normalSubgroups, normalSubgroupsLevelip1);
  od;
  return normalSubgroups;
end);

#The next function collects the data arising from the Jenning series of the (major part of) canonical normal subgroups.
BindGlobal("NormalSubgroupsInfo", function(G)
  local normalSubgroups,Gprime;
  Gprime := DerivedSubgroup(G);
  return List(CanonicalNormalSubgroups(G), x -> [x[1], OnlyJenningsInfo(x[2]), CenterNormalInfo(G, x[2]), LNInfo(G, x[2], Gprime)]);
end);
#####################
#####################

## now all auxilary functions are collected

################################################################################################################################3
##
#O MIPBinsByGTInternal( p, n, L, coh )
##
## Function applying all the group theoretical invariants to the groups contained in the list of lists L. The groups must be of order p^n. The last input is a boolean deciding on whether the 2-cohomology is computed
BindGlobal( "MIPBinsByGTInternal", function( p, n, L, coh )
local bins, cent, i;

    if not IsPrimeInt(p) then
        Error("<p> must be a prime integer");
    fi;
    
    if not IsList(L) or not (ForAll(Flat(L), IsInt) or ForAll(Flat(L), IsGroup)) then
        Error("<L> must be a list of integers or a list of groups");
    fi;
    bins := [ShallowCopy(L)];

    # refine by abelian invariants
    Info(InfoModIsom, 1, "refine by abelian invariants of group");
    bins := RefineBins( p^n, bins, AbelianInvariants, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # refine by center
    Info(InfoModIsom, 1, "refine by abelian invariants of center (Sehgal/Ward)");
    cent := function(G) return AbelianInvariants(Center(G)); end;
    bins := RefineBins( p^n, bins, cent, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    ### NEW. 
    # refine by center and derived subgroup interplay
    Info(InfoModIsom, 1, "refine by abelian invariants of center with derived subgroup (Sandling)");
    bins := RefineBins( p^n, bins, CenterDerivedInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # refine by lower central series
    # Comment changed. This is really more accurate
    Info(InfoModIsom, 1, "refine by small group algebra (Sandling/Baginski/Brenner+Garcia-Lucas)"); 
    bins := RefineBins( p^n, bins, SandlingInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # CHANGED
    # refine by Jennings series
    Info(InfoModIsom, 1, "refine by Jennings series (Passi+Sehgal/Ritter+Sehgal/Hertweck)");
    bins := RefineBins( p^n, bins, JenningsInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Jennings series of derived group
    Info(InfoModIsom, 1, "refine by Jennings series of derived group (Sandling)");
    bins := RefineBins( p^n, bins, JenningsDerivedInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Baginski critiria
    Info(InfoModIsom, 1, "refine by Baginski criteria");
    bins := RefineBins( p^n, bins, BaginskiInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Baginski-Caranti critiria NilpotencyClass(G/Phi(G'))
    Info(InfoModIsom, 1, "refine by Baginski-Caranti invariant");
    bins := RefineBins( p^n, bins, BaginskiCarantiInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Nilpotency class
    # From Baginski-Konovalov 05
    Info(InfoModIsom, 1, "refine by nilpotency class for some cases (Baginski+Konovalov)"); 
    bins := RefineBins( p^n, bins, NilpotencyClassInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW NEW
    # refine by Omega and Agemo from MargolisSakuraiStanojkovski23 (Theorem B) and Garcial-Lucas22 (Examples 3.7)
    Info(InfoModIsom, 1, "refine by Omega and Agemo (Margolis+Stanojkovski+Sakurai/Garcia-Lucas)");  
    bins := RefineBins( p^n, bins, AgemoInvariantAllM, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    if bins <> [ ] then
        bins := RefineBins( p^n, bins, OmegaInvariantAllM, false );
        bins := Filtered( bins, x -> Length(x)>1 );
    fi;
    if bins <> [ ] then
        bins := RefineBins( p^n, bins, AgemoCenterInvariantAllM, false );
        bins := Filtered( bins, x -> Length(x)>1 );
    fi;
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW NEW
    # refine by Theorem 4.1 of Margolis-Stanojkovski 22
    Info(InfoModIsom, 1, "refine by lower central series (Margolis+Stanojkovski)");
    bins := RefineBins( p^n, bins, Theorem41MS22, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW NEW
    # refine by results on p odd with cyclic derived from GarciaLucas-del Rio-Stanojkovski 22 (Theorem 4.2 and Corollaries D, E) and GarciaLucas-del Rio 24 (Theorem A)
    Info(InfoModIsom, 1, "invariants for cyclic derived subgroup (Garcia-Lucas+del Rio+Stanojkovski)");
    bins := RefineBins( p^n, bins, CyclicDerivedInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    #NEW NEW
    # refine by maximal abelian direct factor, GarciaLucas24  
    Info(InfoModIsom, 1, "refine by maximal abelian direct factor (Garcia-Lucas)"); 
    bins := RefineBins( p^n, bins, MaximalAbelianDirectFactor, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    #NEW NEW
    # refine by lattice of canonical subgroups, GarciaLucas24
    Info(InfoModIsom, 1, "refine by lattice of canonical subgroups (Garcia-Lucas)"); 
    bins := RefineBins( p^n, bins, NormalSubgroupsInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW and CHANGED
    # refine by known theoretical results
    Info(InfoModIsom, 1, "checking if covered by theory");
    for i in [1..Size(bins)] do
        if ForAll(bins[i], IsInt) then
            bins[i] := Filtered(bins[i], x -> not IsCoveredByTheory(SmallGroup(p^n, x)));
        else
            bins[i] := Filtered(bins[i], x -> not IsCoveredByTheory(x));
        fi; 
    od;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;


    #NEW and CHANGED
    # refine by dimension of cohomology  
    if coh then
        Info(InfoModIsom, 1, "refine by second cohomology group"); 
        bins := RefineBins( p^n, bins, DimensionTwoCohomology, false );
        bins := Filtered( bins, x -> Length(x)>1 );
        if bins <> [ ] then
            if InfoLevel(InfoModIsom) >= 1 then
                bins := RefineBins( p^n, bins, DimensionSecondHochschild, true );
            else
                bins := RefineBins( p^n, bins, DimensionSecondHochschild, false );
            fi;
        fi;
        bins := Filtered( bins, x -> Length(x)>1 );
        Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
        if Length(bins)=0 then return bins; fi;
    fi;    

    # CHANGED
    # refine by conjugacy classes
    Info(InfoModIsom, 1, "refine by conjugacy classes (Parmenter+Polcino-Milies,Kuelshammer,Roggenkamp/Wursthorn)");
    if InfoLevel(InfoModIsom) >= 1 then
        bins := RefineBins( p^n, bins, ConjugacyClassInfo, true );
    else
        bins := RefineBins( p^n, bins, ConjugacyClassInfo, false );
    fi;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # refine by subgroups
    Info(InfoModIsom, 1, "refine by elem-ab subgroups (Quillen)");
    if InfoLevel(InfoModIsom) >= 1 then
        bins := RefineBins( p^n, bins, SubgroupsInfo, true );
    else
        bins := RefineBins( p^n, bins, SubgroupsInfo, false );
    fi;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups \n");
    if Length(bins)=0 then return bins; fi;

    return bins;
end );

BindGlobal("BinsByGT", function(L...)
local p, n;

    p := L[1];
    n := L[2];
    if Length(L) = 2 then
        return MIPBinsByGTInternal(p, n, [1..NumberSmallGroups(p^n)], true);
    elif Length(L) = 3 then
        return MIPBinsByGTInternal(p, n, L[3], true);
    elif Length(L) = 4 then
        return MIPBinsByGTInternal(p, n, L[3], L[4]);
    fi;
end);

#########
### Variations of the main function for all fields, not just the prime field. Input as before
BindGlobal( "MIPBinsByGTAllFieldsInternal", function( p, n, L, coh )
local bins, cent, i;

    if not IsPrimeInt(p) then
        Error("<p> must be a prime integer");
    fi;
    
    if not IsList(L) or not (ForAll(Flat(L), IsInt) or ForAll(Flat(L), IsGroup)) then
        Error("<L> must be a list of integers or a list of groups");
    fi;
    bins := [ShallowCopy(L)];

    # refine by abelian invariants
    Info(InfoModIsom, 1, "refine by abelian invariants of group");
    bins := RefineBins( p^n, bins, AbelianInvariants, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # refine by center
    Info(InfoModIsom, 1, "refine by abelian invariants of center (Sehgal/Ward)");
    cent := function(G) return AbelianInvariants(Center(G)); end;
    bins := RefineBins( p^n, bins, cent, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    ### NEW. 
    # refine by center and derived subgroup interplay
    Info(InfoModIsom, 1, "refine by abelian invariants of center with derived subgroup (Sandling)");
    bins := RefineBins( p^n, bins, CenterDerivedInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # CHANGED
    # refine by Jennings series
    Info(InfoModIsom, 1, "refine by Jennings series");
    bins := RefineBins( p^n, bins, JenningsInfoAllFields, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Jennings series of derived group
    Info(InfoModIsom, 1, "refine by Jennings series of derived group (Sandling)");
    bins := RefineBins( p^n, bins, JenningsDerivedInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW
    # refine by Nilpotency class
    # From Baginski-Konovalov 05
    Info(InfoModIsom, 1, "refine by nilpotency class for some cases (Baginski+Konovalov)"); 
    bins := RefineBins( p^n, bins, NilpotencyClassInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW NEW
    # refine by Omega and Agemo from MargolisSakuraiStanojkovski23 (Theorem B) and Garcial-Lucas22 (Examples 3.7)
    Info(InfoModIsom, 1, "refine by Omega and Agemo (Margolis+Stanojkovski+Sakurai/Garcia-Lucas)");  
    bins := RefineBins( p^n, bins, AgemoInvariantAllM, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    if bins <> [ ] then  
        bins := RefineBins( p^n, bins, OmegaInvariantAllM, false );
        bins := Filtered( bins, x -> Length(x)>1 );
    fi;
    if bins <> [ ] then
        bins := RefineBins( p^n, bins, AgemoCenterInvariantAllM, false );
        bins := Filtered( bins, x -> Length(x)>1 );
    fi;
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW NEW
    # refine by results on p odd with cyclic derived from GarciaLucas-del Rio-Stanojkovski 22 (Theorem 4.2 and Corollaries D, E) and GarciaLucas-del Rio 24 (Theorem A)
    Info(InfoModIsom, 1, "lower invariants for cyclic derived (Garcia-Lucas+del Rio+Stanojkovski)");
    bins := RefineBins( p^n, bins, CyclicDerivedInfoAllFields, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    #NEW NEW
    # refine by maximal elementary abelian direct factor, MargolisSakuraiStanojkovski23  
    Info(InfoModIsom, 1, "refine by maximal elementary abelian direct factor (Margolis+Sakurai+Stanojkovski)"); 
    bins := RefineBins( p^n, bins, MaximalElementaryAbelianDirectFactor, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    #NEW NEW
    # refine by lattice of canonical subgroups, GarciaLucas24
    Info(InfoModIsom, 1, "refine by lattice of canonical subgroups (Garcia-Lucas)"); 
    bins := RefineBins( p^n, bins, NormalSubgroupsInfo, false );
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # NEW and CHANGED
    # refine by known theoretical results
    Info(InfoModIsom, 1, "checking if covered by theory");
    for i in [1..Size(bins)] do
        if ForAll(bins[i], IsInt) then
            bins[i] := Filtered(bins[i], x -> not IsCoveredByTheoryAllFields(SmallGroup(p^n, x)));
        else
            bins[i] := Filtered(bins[i], x -> not IsCoveredByTheoryAllFields(x));
        fi; 
    od;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;


    #NEW and CHANGED
    # refine by dimension of cohomology  
    if coh then
        Info(InfoModIsom, 1, "refine by second cohomology group"); 
        bins := RefineBins( p^n, bins, DimensionTwoCohomology, false );
        bins := Filtered( bins, x -> Length(x)>1 );
        if bins <> [ ] then
            if InfoLevel(InfoModIsom) >= 1 then
                bins := RefineBins( p^n, bins, DimensionSecondHochschild, true );
            else
                bins := RefineBins( p^n, bins, DimensionSecondHochschild, false );
            fi;
        fi;
        bins := Filtered( bins, x -> Length(x)>1 );
        Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
        if Length(bins)=0 then return bins; fi;
    fi;

    # CHANGED
    # refine by conjugacy classes
    Info(InfoModIsom, 1, "refine by conjugacy classes (Parmenter+Polcino-Milies,Kuelshammer,Roggenkamp/Wursthorn)");
    if InfoLevel(InfoModIsom) >= 1 then
        bins := RefineBins( p^n, bins, ConjugacyClassInfo, true );
    else
        bins := RefineBins( p^n, bins, ConjugacyClassInfo, false );
    fi;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups");
    if Length(bins)=0 then return bins; fi;

    # refine by subgroups
    Info(InfoModIsom, 1, "refine by elem-ab subgroups (Quillen)");
    if InfoLevel(InfoModIsom) >= 1 then
        bins := RefineBins( p^n, bins, SubgroupsInfo, true );
    else
        bins := RefineBins( p^n, bins, SubgroupsInfo, false );
    fi;
    bins := Filtered( bins, x -> Length(x)>1 );
    Info(InfoModIsom, 1, Length(bins)," bins with ",Length(Flat(bins))," groups \n");
    if Length(bins)=0 then return bins; fi;

    return bins;
end );

BindGlobal("BinsByGTAllFields", function(L...)
local p, n;

    p := L[1];
    n := L[2];
    if Length(L) = 2 then
        return MIPBinsByGTAllFieldsInternal(p, n, [1..NumberSmallGroups(p^n)], true);
    elif Length(L) = 3 then
        return MIPBinsByGTAllFieldsInternal(p, n, L[3], true);
    elif Length(L) = 4 then
        return MIPBinsByGTAllFieldsInternal(p, n, L[3], L[4]);
    fi;
end);


#####################################################################################################
##
##  MIPSplitGroupsByGroupTheoreticalInvariants ( L )
##
#### The functions, so that one does not need to input the orders of the groups. 
#### Input: list of groups
#### Output: the bins which one gets from L after applying all the implemented group-theoretical invariants
BindGlobal( "MIPSplitGroupsByGroupTheoreticalInvariants", function( L )
local G, s, p, n;
    for G in L do
        if not IsGroup(G) then
            Error("Input must be a list of groups.");
        fi;
    od;
    s := Size(L[1]);
    if not IsPrimePowerInt(s) then 
        Error("Order of groups in the input must be a prime power order.");
    fi;
    for G in L do
        if not Size(G) = s then
            Error("Input must consists of groups of the same order.");
        fi;
    od;
    p := PrimeDivisors(s)[1];
    n := Log(s, p);
    return MIPBinsByGTInternal(p, n, L, true);
end );

## similar as before without applying 2-cohomology
BindGlobal( "MIPSplitGroupsByGroupTheoreticalInvariantsNoCohomology", function( L )
local G, s, p, n;
    for G in L do
        if not IsGroup(G) then
            Error("Input must be a list of groups.");
        fi;
    od;
    s := Size(L[1]);
    if not IsPrimePowerInt(s) then 
        Error("Order of groups in the input must be a prime power order.");
    fi;
    for G in L do
        if not Size(G) = s then
            Error("Input must consists of groups of the same order.");
        fi;
    od;
    p := PrimeDivisors(s)[1];
    n := Log(s, p);
    return MIPBinsByGTInternal(p, n, L, false);
end );


## similar as before, only utilizing invariants which apply over all fields
BindGlobal( "MIPSplitGroupsByGroupTheoreticalInvariantsAllFields", function( L )
local G, s, p, n;
    for G in L do
        if not IsGroup(G) then
            Error("Input must be a list of groups.");
        fi;
    od;
    s := Size(L[1]);
    if not IsPrimePowerInt(s) then 
        Error("Order of groups in the input must be a prime power order.");
    fi;
    for G in L do
        if not Size(G) = s then
            Error("Input must consists of groups of the same order.");
        fi;
    od;
    p := PrimeDivisors(s)[1];
    n := Log(s, p);
    return MIPBinsByGTAllFieldsInternal(p, n, L, true);
end );


## similar as before, only invariants over all fields and no 2-cohomology
BindGlobal( "MIPSplitGroupsByGroupTheoreticalInvariantsAllFieldsNoCohomology", function( L )
local G, s, p, n;
    for G in L do
        if not IsGroup(G) then
            Error("Input must be a list of groups.");
        fi;
    od;
    s := Size(L[1]);
    if not IsPrimePowerInt(s) then 
        Error("Order of groups in the input must be a prime power order.");
    fi;
    for G in L do
        if not Size(G) = s then
            Error("Input must consists of groups of the same order.");
        fi;
    od;
    p := PrimeDivisors(s)[1];
    n := Log(s, p);
    return MIPBinsByGTAllFieldsInternal(p, n, L, false);
end );



