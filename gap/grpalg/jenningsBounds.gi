###########################################################################################################################################
##########################################################################################################################################3
#############################################################################
##
#O JenningsBound( p, n, bin )
##
## input prime p, number n and list bin of groups of order p^n or group ids of order p^n.
## output: smallest number r such that for some of the groups G and H in bin the quotients of G and H modulo the (r+1)-th dimension subgroup are not isomorphic
BindGlobal( "JenningsBound", function( p, n, bin )
local grps, js, i, ids, l, j;

    if not IsPrimeInt(p) then
        Error("<p> must be a prime integer");
    fi;

    # Check whether we have a list of integers or a list of groups
    if ForAll(bin, IsInt) then
        grps := List(bin, x -> SmallGroup(p^n, x));
    elif ForAll(bin, IsGroup) then
        grps := bin;
    else
        Error("<L> must be a list of integers or a list of groups");
    fi;

    js := List(grps, x -> JenningsSeries(x));
	
#    for i in Reversed([1..Length(js[1])]) do
#        ids := List(js, x -> IdSmallGroup(x[1]/x[i]));
#        if Length(Set(ids))=1 then
#            return i;
#        fi;
#    od;

    l := Minimum(List(js, Length));
	
    for i in [1..l] do
        if i = l then # this assumes that input groups are pairwise not isomorphic
          return l-1;
        fi;
        if ID_AVAILABLE(Size(js[1][1]/js[1][i])) <> fail then
            ids := List(js, x -> IdSmallGroup(x[1]/x[i]));
        else
            ids := [1];
            for j in [2..Size(grps)] do
                if IsomorphismGroups(js[1][1]/js[1][i], js[j][1]/js[j][i]) = fail then
                    ids := [1,2];
                    break;
                fi;
            od;
        fi;
        if Length(Set(ids))>1 then
            return i-1;
        fi;
    od;

end );


#############################################################################
##
#O JenningsBoundPairwise( p, n, bin )
##
## input prime p, number n and list bin of groups of order p^n or group ids of order p^n.
## output: smallest number r such that for all the groups G and H in bin the quotients of G and H modulo the (r+1)-th dimension subgroup are not isomorphic
BindGlobal( "JenningsBoundPairwise", function( p, n, bin )
local grps, js, i, ids, c, ret, idav;

    if not IsPrimeInt(p) then
        Error("<p> must be a prime integer");
    fi;

    if ID_AVAILABLE(p^n) <> fail then
        idav := true;
    else
        idav := false;	      
    fi;

    # Check whether we have a list of integers or a list of groups
    if ForAll(bin, IsInt) then
        grps := List(bin, x -> SmallGroup(p^n, x));
    elif ForAll(bin, IsGroup) then
        grps := bin;
    else
        Error("<L> must be a list of integers or a list of groups");
    fi;

    ret := [];

    for c in Combinations([1..Length(grps)], 2) do
        if idav then
            Add(ret, [IdSmallGroup(grps[c[1]]), IdSmallGroup(grps[c[2]]), JenningsBound(p, n, grps{c})]);
        elif "IdGroup" in KnownAttributesOfObject(grps[c[1]]) and "IdGroup" in KnownAttributesOfObject(grps[c[2]]) then
            Add(ret, [IdSmallGroup(grps[c[1]]), IdSmallGroup(grps[c[2]]), JenningsBound(p, n, grps{c})]);
        else
            Add(ret, [grps[c[1]], grps[c[2]], JenningsBound(p, n, grps{c})]);
        fi;
    od;

    return ret;
end );



