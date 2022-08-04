BindGlobal( "Kur_3_3_Q", function() 
local A, i; 
A := rec(); 
A.rnk := 3; 
A.dim := 77; 
A.com := false; 
A.fld := Rationals; 
A.wgs := [ 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
  5, 5, 5, 5, 5 ]; 
A.wds := [ ,,, [ 1, 1 ], [ 2, 1 ], [ 3, 1 ], [ 1, 2 ], [ 2, 2 ], [ 3, 2 ], 
  [ 1, 3 ], [ 2, 3 ], [ 3, 3 ], [ 2, 4 ], [ 3, 4 ], [ 1, 5 ], [ 2, 5 ], 
  [ 3, 5 ], [ 1, 6 ], [ 2, 6 ], [ 3, 6 ], [ 2, 7 ], [ 3, 7 ], [ 3, 8 ], 
  [ 1, 9 ], [ 2, 9 ], [ 3, 9 ], [ 2, 10 ], [ 3, 10 ], [ 3, 11 ], [ 1, 13 ], 
  [ 2, 13 ], [ 3, 13 ], [ 1, 14 ], [ 2, 14 ], [ 3, 14 ], [ 2, 15 ], 
  [ 3, 15 ], [ 3, 16 ], [ 1, 17 ], [ 2, 17 ], [ 3, 17 ], [ 2, 18 ], 
  [ 3, 18 ], [ 3, 19 ], [ 2, 21 ], [ 3, 21 ], [ 3, 22 ], [ 1, 23 ], 
  [ 2, 23 ], [ 3, 23 ], [ 2, 24 ], [ 3, 24 ], [ 3, 25 ], [ 3, 27 ], 
  [ 3, 28 ], [ 3, 29 ], [ 2, 30 ], [ 3, 30 ], [ 3, 31 ], [ 2, 33 ], 
  [ 3, 33 ], [ 3, 34 ], [ 2, 36 ], [ 3, 36 ], [ 3, 37 ], [ 2, 38 ], 
  [ 3, 38 ], [ 2, 39 ], [ 3, 39 ], [ 3, 40 ], [ 3, 43 ], [ 3, 44 ], 
  [ 2, 48 ], [ 3, 48 ], [ 3, 49 ], [ 3, 52 ], [ 3, 53 ] ]; 
A.tab := [  ]; 
A.tab[1] := [ [ 1, 4, 1 ], [ 2, 7, 1 ], [ 3, 10, 1 ], [ 5, 15, 1 ], 
  [ 6, 18, 1 ], [ 7, 13, -1 ], [ 7, 15, -1 ], [ 8, 16, -1 ], [ 8, 21, -1 ], 
  [ 9, 24, 1 ], [ 10, 14, -1 ], [ 10, 18, -1 ], [ 11, 17, -1 ], 
  [ 11, 19, -1 ], [ 11, 22, -1 ], [ 11, 24, -1 ], [ 11, 27, -1 ], 
  [ 12, 20, -1 ], [ 12, 28, -1 ], [ 13, 30, 1 ], [ 14, 33, 1 ], 
  [ 15, 30, -1 ], [ 16, 31, -1 ], [ 16, 36, -1 ], [ 17, 39, 1 ], 
  [ 18, 33, -1 ], [ 19, 32, -1 ], [ 19, 34, -1 ], [ 19, 37, -1 ], 
  [ 19, 39, -1 ], [ 19, 42, -1 ], [ 20, 35, -1 ], [ 20, 43, -1 ], 
  [ 21, 31, 1 ], [ 22, 32, 2 ], [ 22, 37, 1 ], [ 22, 39, 1 ], [ 23, 48, 1 ], 
  [ 24, 32, -1 ], [ 24, 39, -1 ], [ 25, 38, -1 ], [ 25, 40, -1 ], 
  [ 25, 46, -1 ], [ 25, 48, -1 ], [ 25, 51, -1 ], [ 26, 47, -1 ], 
  [ 26, 52, -1 ], [ 27, 32, -1 ], [ 27, 34, 1 ], [ 27, 37, -1 ], 
  [ 27, 39, -1 ], [ 28, 35, 1 ], [ 29, 41, -1 ], [ 29, 44, -1 ], 
  [ 29, 47, 1 ], [ 29, 54, -1 ], [ 31, 57, -1 ], [ 32, 58, -1 ], 
  [ 34, 60, -1 ], [ 35, 61, -1 ], [ 38, 59, -1 ], [ 39, 58, 1 ], 
  [ 40, 59, 1 ], [ 40, 64, -1 ], [ 40, 68, -1 ], [ 41, 65, -1 ], 
  [ 41, 69, -1 ], [ 44, 62, -1 ], [ 44, 65, 1 ], [ 44, 69, 1 ], 
  [ 45, 63, 1 ], [ 46, 59, 1 ], [ 46, 64, 1 ], [ 47, 69, -1 ], 
  [ 49, 73, -1 ], [ 50, 74, -1 ], [ 51, 59, -1 ], [ 51, 64, -1 ], 
  [ 52, 65, -1 ], [ 53, 67, -1 ], [ 53, 70, -1 ], [ 53, 74, 1 ], 
  [ 54, 62, 1 ], [ 54, 65, 1 ], [ 55, 71, 1 ], [ 56, 76, 1 ] ]; 
A.tab[2] := [ [ 1, 5, 1 ], [ 2, 8, 1 ], [ 3, 11, 1 ], [ 4, 13, 1 ], 
  [ 5, 16, 1 ], [ 6, 19, 1 ], [ 7, 21, 1 ], [ 9, 25, 1 ], [ 10, 27, 1 ], 
  [ 11, 23, -1 ], [ 11, 25, -1 ], [ 12, 26, -1 ], [ 12, 29, -1 ], 
  [ 13, 31, 1 ], [ 14, 34, 1 ], [ 15, 36, 1 ], [ 17, 40, 1 ], [ 18, 42, 1 ], 
  [ 19, 38, -1 ], [ 19, 40, -1 ], [ 20, 41, -1 ], [ 20, 44, -1 ], 
  [ 21, 45, 1 ], [ 22, 38, 2 ], [ 22, 40, 1 ], [ 22, 46, 1 ], [ 23, 49, 1 ], 
  [ 24, 51, 1 ], [ 25, 49, -1 ], [ 26, 50, -1 ], [ 26, 53, -1 ], 
  [ 27, 48, -1 ], [ 27, 51, -1 ], [ 28, 41, 2 ], [ 28, 44, 1 ], 
  [ 28, 54, 1 ], [ 29, 50, 1 ], [ 30, 57, 1 ], [ 32, 59, -1 ], [ 33, 60, 1 ], 
  [ 35, 62, -1 ], [ 36, 63, 1 ], [ 37, 59, 1 ], [ 37, 64, 1 ], [ 38, 66, 1 ], 
  [ 39, 68, 1 ], [ 40, 66, -1 ], [ 41, 67, -1 ], [ 41, 70, -1 ], 
  [ 42, 59, 1 ], [ 42, 68, -1 ], [ 43, 62, 1 ], [ 43, 69, -1 ], 
  [ 44, 67, 1 ], [ 47, 70, -1 ], [ 48, 73, 1 ], [ 50, 75, -1 ], 
  [ 51, 73, -1 ], [ 52, 70, 1 ], [ 52, 74, -1 ], [ 54, 67, 1 ], 
  [ 54, 70, 1 ], [ 55, 72, 1 ], [ 56, 77, 1 ] ]; 
A.tab[3] := [ [ 1, 6, 1 ], [ 2, 9, 1 ], [ 3, 12, 1 ], [ 4, 14, 1 ], 
  [ 5, 17, 1 ], [ 6, 20, 1 ], [ 7, 22, 1 ], [ 8, 23, 1 ], [ 9, 26, 1 ], 
  [ 10, 28, 1 ], [ 11, 29, 1 ], [ 13, 32, 1 ], [ 14, 35, 1 ], [ 15, 37, 1 ], 
  [ 16, 38, 1 ], [ 17, 41, 1 ], [ 18, 43, 1 ], [ 19, 44, 1 ], [ 21, 46, 1 ], 
  [ 22, 47, 1 ], [ 23, 50, 1 ], [ 24, 52, 1 ], [ 25, 53, 1 ], [ 27, 54, 1 ], 
  [ 28, 55, 1 ], [ 29, 56, 1 ], [ 30, 58, 1 ], [ 31, 59, 1 ], [ 33, 61, 1 ], 
  [ 34, 62, 1 ], [ 36, 64, 1 ], [ 37, 65, 1 ], [ 38, 67, 1 ], [ 39, 69, 1 ], 
  [ 40, 70, 1 ], [ 42, 69, -1 ], [ 43, 71, 1 ], [ 44, 72, 1 ], 
  [ 45, 66, -1 ], [ 46, 67, -1 ], [ 48, 74, 1 ], [ 49, 75, 1 ], 
  [ 51, 67, 1 ], [ 51, 74, -1 ], [ 52, 76, 1 ], [ 53, 77, 1 ], [ 54, 72, -1 ] 
 ]; 
for i in [1..A.rnk] do 
    A.tab[i] := ListToMat(A.tab[i],A.dim,A.fld); 
od; 
return A; 
end ); 
