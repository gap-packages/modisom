gap> START_TEST("cfstab.tst");

# pcgs of UT(d,q) as pairs [ mat, mat ], by superdiagonal; both components
# of a result must agree
gap> UTPairs := function( d, F )
>      local pcgs, s, i, lam, m;
>      pcgs := [];
>      for s in [1..d-1] do for i in [1..d-s] do for lam in Basis( F ) do
>        m := IdentityMat( d, F );; m[i][i+s] := lam;;
>        m := ImmutableMatrix( F, m );;
>        Add( pcgs, DirectProductElement( [ m, m ] ) );
>      od; od; od;
>      return pcgs;
>    end;;

# Compare with the orbit: the canonical form lies in it and is constant on
# it, its length is p^(|pcgs|-|stab|), and stab fixes the canonical form.
# GF(4) needs several pivots per coordinate.
gap> CheckSubspace := function( pcgs, U, F )
>      local p, one, C, orb;
>      p := Characteristic( F );
>      one := pcgs[1]^0;
>      C := SubspaceCanonicalForm( pcgs, one, U, F );
>      orb := Orbit( Group( List( pcgs, x -> x[2] ) ), U,
>                    OnSubspacesByCanonicalBasis );
>      return MyTriangulizedBaseMat( U * C.tran[2] ) = C.cano
>        and C.tran[1] = C.tran[2]
>        and ForAll( C.stab, x -> x[1] = x[2]
>                and MyTriangulizedBaseMat( C.cano * x[2] ) = C.cano )
>        and Length( orb ) = p^( Length( pcgs ) - Length( C.stab ) )
>        and ForAll( orb, W ->
>              SubspaceCanonicalForm( pcgs, one, W, F ).cano = C.cano );
>    end;;
gap> F := GF(3);; pcgs := UTPairs( 4, F );;
gap> CheckSubspace( pcgs, MyTriangulizedBaseMat( Z(3)^0 * [[1,1,0,2],[0,0,1,1]] ), F );
true
gap> CheckSubspace( pcgs, MyTriangulizedBaseMat( Z(3)^0 * [[1,0,2,1]] ), F );
true
gap> CheckSubspace( pcgs{[3..6]}, MyTriangulizedBaseMat( Z(3)^0 * [[1,2,0,0],[0,0,1,2],[0,0,0,1]] ), F );
true
gap> F := GF(4);; pcgs := UTPairs( 3, F );;
gap> CheckSubspace( pcgs, MyTriangulizedBaseMat( [[Z(4),Z(4)^2,Z(4)^0]] ), F );
true
gap> CheckSubspace( pcgs, MyTriangulizedBaseMat( [[Z(4)^0,Z(4),0*Z(4)],[0*Z(4),0*Z(4),Z(4)^0]] ), F );
true
gap> CheckSubspace( pcgs{[2,3,5,6]}, MyTriangulizedBaseMat( [[Z(4)^0,Z(4),Z(4)^2]] ), F );
true

# elements acting trivially stay in the stabilizer: act on F^4 / <e_4>
gap> F := GF(3);; pcgs := UTPairs( 4, F );;
gap> qpcgs := List( pcgs, x -> DirectProductElement( [ x[1],
>                   ImmutableMatrix( F, x[2]{[1..3]}{[1..3]} ) ] ) );;
gap> V := MyTriangulizedBaseMat( Z(3)^0 * [[1,2,1]] );;
gap> D := SubspaceCanonicalForm( qpcgs, qpcgs[1]^0, V, F );;
gap> qorb := Orbit( Group( List( qpcgs, x -> x[2] ) ), V, OnSubspacesByCanonicalBasis );;
gap> Length( qorb ) * 3^Length( D.stab ) = 3^Length( qpcgs );
true
gap> ForAll( D.stab, x -> MyTriangulizedBaseMat( D.cano * x[2]{[1..3]}{[1..3]} ) = D.cano );
true

# BlockCanonicalForm under the Borel subgroup of GL(4,3): constant on the
# orbit, and G becomes the stabilizer of the canonical form
gap> torus := List( [1..4], function( i ) local m;
>      m := IdentityMat( 4, F );; m[i][i] := -m[i][i];;
>      m := ImmutableMatrix( F, m );;
>      return DirectProductElement( [ m, m ] ); end );;
gap> borel := rec( field := F, one := pcgs[1]^0, glAutos := torus,
>      glPerms := [ (1,2), (3,4), (5,6), (7,8) ], glOrder := 16,
>      agAutos := pcgs, size := 16 * 3^6 );;
gap> U := MyTriangulizedBaseMat( Z(3)^0 * [[1,1,0,2],[0,1,1,1]] );;
gap> orb := Orbit( Group( List( Concatenation( torus, pcgs ), x -> x[2] ) ), U,
>                  OnSubspacesByCanonicalBasis );;
gap> G := ShallowCopy( borel );; cf := BlockCanonicalForm( G, U );;
gap> MyTriangulizedBaseMat( U * cf.tran[2] ) = cf.cano and cf.cano in orb;
true
gap> G.size * Length( orb ) = borel.size;
true
gap> ForAll( Concatenation( G.glAutos, G.agAutos ), x -> x[1] = x[2]
>            and MyTriangulizedBaseMat( cf.cano * x[2] ) = cf.cano );
true
gap> Size( Group( List( Concatenation( G.glAutos, G.agAutos ), x -> x[2] ),
>                 One( borel.one[2] ) ) ) = G.size;
true
gap> ForAll( orb, function( W ) local H;
>      H := ShallowCopy( borel );
>      return BlockCanonicalForm( H, W ).cano = cf.cano and H.size = G.size;
>    end );
true

# without gl part
gap> G := ShallowCopy( borel );; G.glAutos := [];; G.glPerms := [];; G.glOrder := 1;;
gap> cf := BlockCanonicalForm( G, U );;
gap> cf.cano = SubspaceCanonicalForm( pcgs, borel.one, U, F ).cano;
true
gap> G.size * Length( Orbit( Group( List( pcgs, x -> x[2] ) ), U,
>                            OnSubspacesByCanonicalBasis ) ) = 3^6;
true

#
gap> STOP_TEST("cfstab.tst", 1);
