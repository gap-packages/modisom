#######################
## these are functions which in previous versions were inside detbins.gi
## they compute some ring-theoretical invariants, so have been separated from the group-theoretical invariants in detbins.gi
## the functions have not been changed. In particular, they can functions at the end of the file only apply to ids, not to lists of groups.

#############################################################################
##
#F Center and Commutator 
##
BindGlobal( "CoeffsCenterBasis", function(A)
    local G, F, e, c, b, i, j, k;

    # set up
    G := UnderlyingMagma(A);
    F := LeftActingDomain(A);
    e := Elements(G);

    # compute conjugacy classes (except 1^G)
    c := Orbits(G, e); c := Filtered(c, x -> x[1] <> e[1]);

    # a basis of Z(A) cap I
    b := [];
    for i in [1..Length(c)] do

        # set up
        b[i] := List( [1..Length(e)], x -> Zero(F) );

        # determine coeffs of class-sum
        for j in [1..Length(c[i])] do
            k := Position(e, c[i][j]);
            b[i][k] := One(F);
        od;

        # if class-length is > 1, then it is in I - otherwise enforce
        if Length(c[i]) = 1 then b[i][1] := -One(F); fi;
    od;
    return b;
end );

BindGlobal( "CoeffsCommutatorBasis", function(A)
    local G, F, e, c, b, i, j, k, h, v;

    # set up
    G := UnderlyingMagma(A);
    F := LeftActingDomain(A);
    e := Elements(G);

    # compute conjugacy classes of length > 1
    c := Orbits(G, e); 
    c := Filtered(c, x -> Length(x)>1);

    # get basis
    b := [];
    for i in [1..Length(c)] do
        k := Position(e, c[i][1]);
        for j in [2..Length(c[i])] do
            v := List([1..Length(e)], x -> Zero(F));
            h := Position(e,c[i][j]);
            v[k] := One(F);
            v[h] := -One(F);
            Add( b, v );
        od;
    od;
 
    # that's it
    return b;
end );

#############################################################################
##
#F Compute powers for ideals in series
##
BindGlobal( "LinearPowersBySeries", function( A, bases )
    local l, r, i, w, p;

    l := Length(bases);
    r := List([1..l], x -> []);
    p := PrimePGroup(UnderlyingMagma(A));

    for i in [1..l] do

        # compute powers iteratedly
        w := StructuralCopy(bases[i]);
        while Length(w) > 0 do
            w := List(w, x -> (x*Basis(A))^p);
            w := MyTriangulizedBaseMat(List(w, x -> Coefficients(Basis(A),x)));
            Add(r[i], w);
        od;
        Unbind(r[i][Length(r[i])]);

        # if there are no powers, then we are done
        if Length(r[i]) = 0 then return Filtered(r, x -> Length(x)>0); fi;
    od;

    return Filtered(r, x -> Length(x)>0); 
end );

#############################################################################
##
#F Compute maximal abelian factors
##
BindGlobal( "IsAbFac", function( I, J )
    local b, i, j;

    b := Basis(I);
    for i in [1..Length(b)] do
        for j in [i+1..Length(b)] do
            if not (b[i]*b[j] - b[j]*b[i]) in J then 
                return false;
            fi;
        od;
    od;
    return true;
end );

BindGlobal( "MaximalAbelianFactors", function( A, ids )
    local f, i, n;

    n := Length(ids);
    f := [];
    for i in [1..n] do

        # find largest index j with ids[i]/ids[j] abelian
        f[i] := First(Reversed([i+1..n]), x -> IsAbFac(ids[i],ids[x]));

        # if j = n, then we are done
        if f[i] = n then 
            return Concatenation( f, List([i+1..n], x -> n) );
        fi;
    od;
end );

#############################################################################
##
#F power map on center
##
BindGlobal( "PowerMapCenter", function(A)
    local s, c, d, l, t, n, k, r;

    # set up
    s := CoeffsWeightedBasisOfRad(A);
    c := CoeffsCenterBasis(A);
    d := Dimension(A)-1;
    l := s.weights[d];

    # compute series of Z(A) cap I^n 
    t := [];
    for n in [1..l] do
        k := Position(s.weights,n);
        t[n] := SumIntersectionMat(s.basis{[k..d]}, c)[2];
    od;

    # compute power maps
    r := LinearPowersBySeries( A, t );

    # return dimensions
    return List(r, x -> List(x, Length) );
end );

#############################################################################
##
#F power map on abelian factors
##
BindGlobal( "PowerMapAbelian", function(A)
    local s, d, l, ids, bas, res, bcm, bct, b, I, f, r, cm, ct, i, j, k, n;

    # set up
    s := CoeffsWeightedBasisOfRad(A);
    d := Dimension(A)-1;
    l := s.weights[d];

    # get ideals and bases
    ids := [];
    bas := [];
    for n in [1..l] do
        k := Position(s.weights, n);
        b := s.basis{[k..d]};
        I := SubalgebraNC( A, b*Basis(A), "basis"); 
        SetIsIdealInParent(I, true);
        Add( ids, I );
        Add( bas, b );
    od;
    Add( ids, SubalgebraNC( A, [] ) );
    Add( bas, [] );

    # determine maximal abelian factors
    #Print("get maximal abelian factors \n");
    f := MaximalAbelianFactors( A, ids );

    # get power maps
    #Print("get linear powers \n");
    r := LinearPowersBySeries( A, bas );

    # some more infos
    cm := CoeffsCommutatorBasis(A);
    ct := CoeffsCenterBasis(A);

    # evaluate result
    res := [];
    for i in [1..Length(r)] do
        res[i] := [];
        for j in [i+1..f[i]] do
            res[i][j] := [];
            for k in [1..Length(r[i])] do
                #Print("  consider ",i," mod ",j," with ",k,"th powers \n");

                # take image
                b := MyTriangulizedBaseMat( Concatenation( r[i][k], bas[j] ) );

                # relate to comms
                bcm := RankMat(Concatenation( cm, b ));

                # relate to cent
                bct := RankMat(Concatenation( ct, b ));
            
                # store info
                res[i][j][k] := [Length(b), bcm, bct];
            od;
        od;
    od;

    return res;
end );
    
#############################################################################
##
#F power map on small factors
##
BindGlobal( "MyIsMember", function( bas, u )
    if Length(bas) = 0 then return (u = 0*u); fi;
    return not IsBool(SolutionMat(bas, u));
end );

BindGlobal( "PowerMapKernels", function( A, bas, n, m )
    local p, dim, siz, ppp, fac, k, q, l, v, w, u, h;

    p := PrimePGroup(UnderlyingMagma(A));

    # set up
    dim := Length(bas[n]) - Length(bas[n+m]);
    ppp := List( [1..dim], x -> p );
    fac := bas[n]{[1..dim]};
    siz := [];

    # loop
    k := 0;
    repeat
        k := k+1;
        q := p^k;
        l := Minimum(Length(bas), n*q+m);
        siz[k] := 0;
        for h in [0..p^dim-1] do
            v := CoefficientsMultiadic( ppp, h ) * fac;
            w := v * Basis(A);
            u := Coefficients(Basis(A), w^q);
            if MyIsMember(bas[l], u) then 
                siz[k] := siz[k]+1; 
            fi;
        od;
    until siz[k] = p^dim;
 
    return siz;
end );
            
BindGlobal( "PowerMapSmall", function(A)
    local s, d, l, p, bas, n, k, b, res, i, j, dim;

    # set up
    s := CoeffsWeightedBasisOfRad(A);
    d := Dimension(A)-1;
    l := s.weights[d];
    p := PrimePGroup(UnderlyingMagma(A));
    #Print(s.weights,"\n");

    # get ideals and bases
    bas := [];
    for n in [1..l] do
        k := Position(s.weights, n);
        b := s.basis{[k..d]};
        Add( bas, b );
    od;
    Add( bas, [] );

    # evaluate result
    res := List([1..l], x -> []);
    for i in [1..l] do
        for j in [i+1..l] do
            dim := Length(bas[i]) - Length(bas[j]);
            #Print("   consider ",i," with ",j," of dim ",dim,"\n");
            if p^dim <= POWER_LIMIT then 
                res[i][j-i] := PowerMapKernels( A, bas, i, j-i );
            fi;
        od;
    od;

    return res;
end );

BindGlobal( "PowerMapFLSmall", function(A)
    local s, d, l, p, bas, n, k, b, res, j, dim;

    # set up
    s := CoeffsWeightedBasisOfRad(A);
    d := Dimension(A)-1;
    l := s.weights[d];
    p := PrimePGroup(UnderlyingMagma(A));

    # get ideals and bases
    bas := [];
    for n in [1..l] do
        k := Position(s.weights, n);
        b := s.basis{[k..d]};
        Add( bas, b );
    od;
    Add( bas, [] );

    # evaluate result
    res := [];
    for j in [2..l] do
        dim := Length(bas[1]) - Length(bas[j]);
        if p^dim <= POWER_LIMIT then 
            res[j-1] := PowerMapKernels( A, bas, 1, j-1 );
        fi;
    od;

    return res;
end );

#############################################################################
##
#F Refine bins with ring theory
##
BindGlobal( "RefineBinByFunc", function( bin, obj, func )
    local res, inv, new, j, k;

    res := List( [1..Length(bin)], x -> func(obj[x]) );
    inv := Set(res);

    # split bin
    new := List(inv, x -> []);
    for j in [1..Length(res)] do
        k := Position(inv,res[j]);
        Add( new[k], bin[j] );
    od;
    return new;
end );

BindGlobal( "RefineBinByRT", function( bin, alg )
    local new, i, w, l, j, d, func;

    # translate
    new := [[1..Length(bin)]];

    # get weights
    w := WeightedBasisOfRad(alg[1]).weights;
    l := w[Length(w)];

    # 1. step: refine by p-powers on center
    Print("  refine by p-power map on center \n");
    for i in [1..Length(new)] do
        new[i] := RefineBinByFunc( new[i], alg{new[i]}, PowerMapCenter );
    od;
    new := Filtered( Concatenation(new), x -> Length(x) > 1 );
    Print("  -- ",Length(new)," sub-bins \n");
    if Length(new)=0 then return []; fi;

    # next step: refine by p-powers on small factors
    Print("  refine by p-power map on first layer small factors \n");
    POWER_LIMIT := 500;
    for i in [1..Length(new)] do
	new[i] := RefineBinByFunc( new[i], alg{new[i]}, PowerMapFLSmall );
    od;
    new := Filtered( Concatenation(new), x -> Length(x) > 1 );
    Print("  -- ",Length(new)," sub-bins \n");
    if Length(new)=0 then return []; fi;

#    # refine by p-powers on center
#    POWER_LIMIT := 2000;
#    Print("  refine by p-power map on all layers medium factors \n");
#    for i in [1..Length(new)] do
#	new[i] := RefineBinByFunc( new[i], alg{new[i]}, PowerMapSmall );
#    od;
#    new := Filtered( Concatenation(new), x -> Length(x) > 1 );
#    Print("  -- ",Length(new)," sub-bins \n");

    # translate remaining bins
    return List(new, x -> bin{x});
end );

BindGlobal( "RefineBinsByRT", function( p, n, bins )
    local news, grps, algs, i, j;

    news := [];
    for i in [1..Length(bins)] do

	Print("start bin ",i," of ",Length(bins));
	Print(" with ",Length(bins[i])," cands \n");

	# set up
	grps := List( bins[i], x -> SmallGroup(p^n,x) );
	algs := List( grps, x -> GroupRing(GF(p),x) );

	# refine this bin
	news[i] := RefineBinByRT( bins[i], algs );
    od;

    return Concatenation(news);
end );


