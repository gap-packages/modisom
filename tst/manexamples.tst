gap> START_TEST("ModIsom package: manexamples.tst");
gap> LoadPackage("modisom",false);
true

#####################################################################
# Chapter 2
#####################################################################
gap> T := rec( dim := 3, 
>               fld := GF(2), 
>               rnk := 2, 
>               wgs := [ 1, 1, 2 ],
>               wds := [ ,, [ 2, 1 ] ],
>               tab := [] );;
gap> T.tab[1] := [[0,0,0],[0,0,1]] * One(T.fld);;
gap> T.tab[2] := [[0,0,1],[0,0,0]] * One(T.fld);;
gap> GetEntryTable( T, 3, 1 );
[ 0*Z(2), 0*Z(2), 0*Z(2) ]

#####################################################################
gap> A := GroupRing(GF(2), SmallGroup(8,3));
<algebra-with-one over GF(2), with 3 generators>
gap> NilpotentTableOfRad(A);
rec( dim := 7, fld := GF(2), rnk := 2, 
  tab := 
    [ 
      [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ], 
      [ [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0, Z(2)^0, Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ],, 
      [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ] ], 
  wds := [ ,, [ 1, 2 ],, [ 1, 4 ], [ 2, 4 ], [ 1, 6 ] ], 
  wgs := [ 1, 1, 2, 2, 3, 3, 4 ] )

#####################################################################
# Chapter 3
#####################################################################
gap> A := GroupRing(GF(2), SmallGroup(8,3));;
gap> T := TableByWeightedBasisOfRad(A);;
gap> C := CanoFormWithAutGroupOfTable(T);;

# check that the canonical form is not equal to T
gap> CompareTables(C.cano, T);
false

# the order of the automorphism group
gap> C.auto.size;
512

# the entries of the canonical table as far as they are bounded
gap> C.cano.tab;
[ [ <a GF2 vector of length 7>, <a GF2 vector of length 7>, 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ], 
  [ <a GF2 vector of length 7>, <a GF2 vector of length 7>, 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ], 
  [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ], 
  [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ], 
  [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ], 
  [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
      [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ] ]

#####################################################################
# Chapter 4
#####################################################################
gap> bins := BinsByGT(2,6);
[ [ 156, 158, 160 ], [ 155, 157 ], [ 173, 176 ], [ 179, 180 ] ]

#####################################################################
gap> MIPSplitGroupsByAlgebras(2,6,bins[1]);
rec( bins := [  ], splits := [ [ 7, [ 156, 158, 160 ] ] ], time := 2243 )

#####################################################################
# Chapter 5
#####################################################################
gap> F := FreeAssociativeAlgebra(GF(2), 2);;
gap> g := GeneratorsOfAlgebra(F);;
gap> r := [g[1]^2, g[2]^2];;
gap> A := F/r;;
gap> NilpotentQuotientOfFpAlgebra(A,3);
rec( def := [ 1, 2 ], dim := 8, fld := GF(2), 
  img := [ <a GF2 vector of length 8>, <a GF2 vector of length 8> ], 
  mat := [ [  ], [  ] ], rnk := 2, 
  tab := 
    [ 
      [ <a GF2 vector of length 8>, <a GF2 vector of length 8>, 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
        , 
      [ <a GF2 vector of length 8>, <a GF2 vector of length 8>, 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
        , 
      [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ] ]
        , 
      [ [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ], 
          [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] 
         ] ], 
  wds := [ ,, [ 2, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 4 ], [ 2, 5 ], [ 1, 6 ] ], 
  wgs := [ 1, 1, 2, 2, 3, 3, 4, 4 ] )

#####################################################################
# Chapter 6
#####################################################################
gap> KuroshAlgebra(2,2,Rationals);
next step with dim 2
  got cover of dim 6
  induced autos 
    subspace has dim 0
    subspace has dim 3
  found subspace of dim 3 in 4
next step with dim 3
  got cover of dim 6
rec( bas := [ [ 1, 0, 0, 0 ], [ 0, 1, 1, 0 ], [ 0, 0, 0, 1 ], [ 0, 1, 0, 0 ] ]
    , com := false, dim := 3, fld := Rationals, rnk := 2, 
  tab := [ [ [ 0, 0, 0 ], [ 0, 0, -1 ], [ 0, 0, 0 ] ], 
      [ [ 0, 0, 1 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ] ], wds := [ ,, [ 2, 1 ] ], 
  wgs := [ 1, 1, 2 ] )
gap> STOP_TEST("manexamples.tst",10000);

