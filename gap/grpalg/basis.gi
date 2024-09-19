
BindGlobal( "CoeffsNatBasisOfAug", function(A)
    local d, F, I, i;
    d := Dimension(A);
    F := LeftActingDomain(A);
    I := IdentityMat(d, F){[2..d]};
    for i in [1..d-1] do I[i][1] := -One(F); od;
    return I;
end );

####################################################################################################################
# The pcgs of the return can be used as a pcgs of G. Weights are in the sense of Jennings basis, 
# i.e. in which maximal power of the augmentation ideal they still appear.
# Input: finite p-group
# Output: The pcgs of subsequent quotients of a Jennings series and their weights
BindGlobal("PcgsJenningsSeries", function( G )
local s, b, w, i, g;

     s := JenningsSeries(G);
     b := [];
     w := [];
     for i in [ 1 .. Length(s)-1 ] do
         g := ModuloPcgs(s[i],s[i+1]);
         Append(b,g);
         Append(w, List(g, x -> i));
     od;

     return rec( pcgs := PcgsByPcSequence(FamilyObj(One(G)), b), weights := w );
end);


BindGlobal( "WeightedBasisOfRad", function(A)
    local G, p, n, js, jb, jw, ww, wb, wc, df, i, h;

    # set up
    G := UnderlyingMagma(A);
    p := PrimePGroup(G);
    n := Length(Pcgs(G));

    # get a special basis of G
    js := PcgsJenningsSeries(G);
    jb := js.pcgs;
    jw := js.weights;
    df := List( jb, x -> (x-One(A)) );

    # set up for basis with weights
    wc := Tuples( [ 0 .. p-1 ], n ); RemoveSet( wc, 0*[1 .. n] );
    wc := List(wc, Reversed);
    ww := [];
    wb := [];

    # determine a basis with weights for A
    for i in [1..Length(wc)] do
        ww[i] := Sum( [1..n], x -> wc[i][x]*jw[x] );
        wb[i] := Product( [1..n], x -> df[x]^wc[i][x] );
    od;

    # sort and return
    h := Sortex(ww);
    return rec( basis := Permuted(wb, h),
                weights := ww,
                exps := Permuted(wc, h) );
end );

BindGlobal( "CoeffsWeightedBasisOfRad", function(A)
    local B;
    B := WeightedBasisOfRad(A);
    B.basis := List(B.basis, x -> Coefficients(Basis(A),x));
    return B;
end );


