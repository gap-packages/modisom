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
#I  refine by abelian invariants of group (Sehgal/Ward)
#I  13 bins with 256 groups
#I  refine by abelian invariants of center (Sehgal/Ward)
#I  30 bins with 237 groups
#I  refine by abelian invariants of center with derived subgroup (Sandling)
#I  39 bins with 234 groups
#I  refine by small group algebra (Sandling/Baginski)
#I  34 bins with 122 groups
#I  refine by Jennings series (Passi+Sehgal/Ritter+Sehgal/Hertweck)
#I  36 bins with 120 groups
#I  refine by Jennings series of derived group (Sandling)
#I  37 bins with 118 groups
#I  refine by Baginski criteria
#I  35 bins with 114 groups
#I  refine by Baginski-Caranti invariant
#I  35 bins with 114 groups
#I  refine by nilpotency class for some cases (Baginski+Konovalov)
#I  35 bins with 114 groups
#I  refine by Omega and Agemo (Margolis+Stanojkovski+Sakurai/Garcia-Lucas)
#I  33 bins with 110 groups
#I  refine by lower central series (Margolis+Stanojkovski)
#I  33 bins with 110 groups
#I  invariants for cyclic derived subgroup (Garcia-Lucas+del Rio+Stanojkovski)
#I  33 bins with 110 groups
#I  refine by maximal abelian direct factor (Garcia-Lucas)
#I  33 bins with 110 groups
#I  refine by lattice of canonical subgroups (Garcia-Lucas)
#I  32 bins with 108 groups
#I  checking if covered by theory
#I  30 bins with 103 groups
#I  refine by second cohomology group
  start bin 1 of 27
  start bin 2 of 27
  start bin 3 of 27
  start bin 4 of 27
  start bin 5 of 27
  start bin 6 of 27
  start bin 7 of 27
  start bin 8 of 27
  start bin 9 of 27
  start bin 10 of 27
  start bin 11 of 27
  start bin 12 of 27
  start bin 13 of 27
  start bin 14 of 27
  start bin 15 of 27
  start bin 16 of 27
  start bin 17 of 27
  start bin 18 of 27
  start bin 19 of 27
  start bin 20 of 27
  start bin 21 of 27
  start bin 22 of 27
  start bin 23 of 27
  start bin 24 of 27
  start bin 25 of 27
  start bin 26 of 27
  start bin 27 of 27
#I  13 bins with 27 groups
#I  refine by conjugacy classes (Parmenter+Polcino-Milies,Kuelshammer,Roggenkamp/Wursthorn)
  start bin 1 of 13
  start bin 2 of 13
  start bin 3 of 13
  start bin 4 of 13
  start bin 5 of 13
  start bin 6 of 13
  start bin 7 of 13
  start bin 8 of 13
  start bin 9 of 13
  start bin 10 of 13
  start bin 11 of 13
  start bin 12 of 13
  start bin 13 of 13
#I  6 bins with 13 groups
#I  refine by elem-ab subgroups (Quillen)
  start bin 1 of 6
  start bin 2 of 6
  start bin 3 of 6
  start bin 4 of 6
  start bin 5 of 6
  start bin 6 of 6
#I  4 bins with 9 groups 

[ [ 156, 158, 160 ], [ 155, 157 ], [ 173, 176 ], [ 179, 180 ] ]


#####################################################################
gap> CheckBin(2,6,100,bins[1]);
compute tables through power series 
  determined table for 1
  determined table for 2

refine bin 
  weights yields bins [ [ 1, 2 ] ]
  layer 1 yields bins [ [ 1, 2 ] ]
  layer 2 yields bins [ [ 1, 2 ] ]
  layer 3 yields bins [ [ 1, 2 ] ]
  layer 4 yields bins [  ]

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

