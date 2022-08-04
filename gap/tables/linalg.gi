BindGlobal( "SiftInto", function( B, c )
    local dep, d, i;

    # catch some simple cases
    if c = 0*c then return false; fi;
    if Length(B) = 0 then B[1] := NormedRowVector(c); return true; fi;
    if Length(B) = Length(B[1]) then return false; fi;

    # sift
    dep := List(B, PositionNonZero);
    while true do
        d := PositionNonZero(c);
        if d > Length(c) then return false; fi;
        i := Position(dep,d);
        if IsBool(i) then
            B[Length(B)+1] := NormedRowVector(c);
            return true;
        fi;
        AddRowVector( c, B[i], -c[d] );
    od;

end );

BindGlobal( "OrderByDepth", function( U )
    SortParallel(List(U, PositionNonZero), U);
    return U;
end );

