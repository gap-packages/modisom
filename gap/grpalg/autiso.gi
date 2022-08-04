
BindGlobal( "CanoFormWithAutGroupOfRad", function(A)
    local T;
    T := TableByWeightedBasisOfRad(A);
    return CanoFormWithAutGroupOfTable(T);
end );

BindGlobal( "CanonicalFormOfRad", function(A)
    return CanoFormWithAutGroupOfRad(A).cano;
end );

BindGlobal( "AutomorphismGroupOfRad", function(A)
    local G, d, F, a, B;
    G := CanoFormWithAutGroupOfRad(A).auto;
    d := Length(G.one);
    F := G.field;
    a := Concatenation(G.glAutos, G.agAutos);
    B := Subgroup(GL(d, F), a);
    SetSize(B, G.size);
    return B;
end );

BindGlobal( "MatPlus", function( mat, F, k )
    local M, i, j;
    M := NullMat( Length(mat)+1, Length(mat)+1, F );
    M[1][k] := One(F);
    for i in [1..Length(mat)] do
        for j in [1..Length(mat)] do
            M[i+1][j+1] := mat[i][j];
        od;
    od;
    return M;
end );

