gap> START_TEST("ModIsom package: kuroshalgebras.tst");

# Version 2.3.3 of the package contains the Kurosh algebras for 
# - n=2 and arbitrary d and F
# - (d,n) = (2,3) and arbitrary F
# - (d,n) = (3,3) and F = Q or |F| ∈ {2,3,4}
# - (d,n) = (4,3) and F = Q or |F| ∈ {2,3,4}
# - (d,n) = (2,4) and F = Q or |F| ∈ {2,3,4,9}
# - (d,n) = (2,5) and F = Q or |F| ∈ {2,3,4,5,8,9}.
#
#
# Check that the library is accessible by KuroshAlgebraByLib. We use
# separate line for each case to simplify reading in case of test errors 
# 
# These two cases are to be addressed later after dealing with infolevels
# - n=2 and arbitrary d and F
# - (d,n) = (2,3) and arbitrary F
#
# (d,n) = (3,3) and F = Q or |F| ∈ {2,3,4}
gap> K:=KuroshAlgebraByLib(3,3,Rationals);;
gap> K:=KuroshAlgebraByLib(3,3,GF(2));;
gap> K:=KuroshAlgebraByLib(3,3,GF(3));;
gap> K:=KuroshAlgebraByLib(3,3,GF(4));;

# (d,n) = (4,3) and F = Q or |F| ∈ {2,3,4}
gap> K:=KuroshAlgebraByLib(4,3,Rationals);;
gap> K:=KuroshAlgebraByLib(4,3,GF(2));;
gap> K:=KuroshAlgebraByLib(4,3,GF(3));;
gap> K:=KuroshAlgebraByLib(4,3,GF(4));;

# (d,n) = (2,4) and F = Q or |F| ∈ {2,3,4,9}
gap> K:=KuroshAlgebraByLib(2,4,Rationals);;
gap> K:=KuroshAlgebraByLib(2,4,GF(2));;
gap> K:=KuroshAlgebraByLib(2,4,GF(3));;
gap> K:=KuroshAlgebraByLib(2,4,GF(4));;
gap> K:=KuroshAlgebraByLib(2,4,GF(9));;

# (d,n) = (2,5) and F = Q or |F| ∈ {2,3,4,5,8,9}.
gap> K:=KuroshAlgebraByLib(2,5,Rationals);;
gap> K:=KuroshAlgebraByLib(2,5,GF(2));;
gap> K:=KuroshAlgebraByLib(2,5,GF(3));;
gap> K:=KuroshAlgebraByLib(2,5,GF(4));;
gap> K:=KuroshAlgebraByLib(2,5,GF(5));;
gap> K:=KuroshAlgebraByLib(2,5,GF(8));;
gap> K:=KuroshAlgebraByLib(2,5,GF(9));;

#
gap> STOP_TEST("manexamples.tst",10000);

