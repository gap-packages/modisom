#############################################################################
##Functions to check conjecture on Jennings bound from MM20
##
## input: prime and JenningsBound between p-groups
## output: conjectured Jennings bound from MM20
BindGlobal( "ConjecturedJenBound", function(p, jb)
if p = 2 then
  return 2*jb;
else
 return 2*jb - ((p-1)/2);
fi;
end );

## 
## 
BindGlobal( "JenningsBoundConjecture", function(p, n, tup)
local s, l, ids, jb, b, pair, jbp, bp, a, lp;
if IsString(tup.splits[1]) or IsString(tup.splits[1][1]) then
  Print("Data is given in format which does not allow the check. \n");
  return fail;
else
  s := tup.splits;
  l := Maximum(List(s, x -> x[1]));
  ids := Concatenation(List(s, x -> x[2]));
  jb := JenningsBound(p,n,ids);
  b := ConjecturedJenBound(p,jb);
  if l <= b then  
    return true;
  else
    for pair in Combinations(ids, 2) do
    jbp := JenningsBoundPairwise(p,n,pair)[1][3];
    bp := ConjecturedJenBound(p, jbp);
      if l > bp then
        a := MIPBinSplit(p,n,false,jbp,1,pair);
        lp := Maximum(List(a.splits, x -> x[1]));
        if lp > bp then
          return false;
        else 
          return true;
        fi;
      fi;
    od;
  fi;
fi;
end );

BindGlobal( "JenningsBoundConjectureIsStrict", function(p, n, tup)
local s, l, ids, jb, b;
if IsString(tup.splits[1]) or IsString(tup.splits[1][1]) then
  Print("Data is given in format which does not allow the check. \n");
  return fail;
else
  s := tup.splits;
  l := Maximum(List(s, x -> x[1]));
  ids := Concatenation(List(s, x -> x[2]));
  if Size(ids) > 3 then
    Print("This function is currently only implemented for pairs of groups \n");
    return fail;
  fi;
  jb := JenningsBound(p,n,ids);
  b := ConjecturedJenBound(p,jb);
  if l = b then  
    return true;
  else 
    return false;
  fi;
fi;
end );
