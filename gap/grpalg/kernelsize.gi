# Computes the size of the kernel of the map I^n/I^(n+m) -> I^n*p^l/I^n*p^l+m; x -> x^(p^l), cf. Hertweck-Soriano Section 2.7
BindGlobal("KernelSizePowerMap", function(T, n, m, l)
local Res, p, max, v0, pos, tup, s, v, w, i, posm, t, count, k;
    Res := 0;
    k := T.fld;
    p := Characteristic(k);
    max := Maximum(T.wgs);
    if max < n*(p^l)+m-1 then
        Print("Weight of table is too small, recompute for weight ", n*(p^l)+m, "\n");
        return fail;
    fi;
    pos := Filtered([1..Size(T.wgs)], x -> n <= T.wgs[x] and T.wgs[x] <= n+m-1);
    s := Size(pos);
    posm := Filtered([1..Size(T.wgs)], x -> T.wgs[x] < n*(p^l)+m);
    t := Size(posm);
    v0 := ListWithIdenticalEntries(t, 0*Z(p));
    count := 0;
    #Act on "projective" space to save time, hence the s-1
    for tup in IteratorOfTuples(k, s-1) do
        count := count + 1;
        if count mod 1000 = 0 and InfoLevel(InfoModIsom) >= 1 then
            Print("Step ", count, " of ", NrTuples(k,s-1), "\n");
        fi;
        #Case last entry is 0
        v := ListWithIdenticalEntries(T.dim, 0*Z(p));
        for i in [1..s-1] do
            v[pos[1]-1+i] := tup[i]*(Z(p)^0);
        od;
        w := v;
        for i in [1..p^l-1] do
            w := MultByTable(T, w, v);
        od;
        if w{[1..t]} = v0 then
            Res := Res + 1; 
        fi;
        #Case last entry is 1
        v := ListWithIdenticalEntries(T.dim, 0*Z(p));
        v[pos[1]-1+s] := Z(p)^0;
        for i in [1..s-1] do
            v[pos[1]-1+i] := tup[i]*Z(p);
        od;
        w := v;
        for i in [1..(p^l)-1] do
            w := MultByTable(T, w, v);
        od;
        if w{[1..t]} = v0 then
        #As any multipel of w is also in the kernel
            Res := Res + Size(k)-1; 
        fi;
    od; 
    return Res;
end);
