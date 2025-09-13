## convert the element el of the group algebra kG to an element of the table T. Of course T and kG should have the same underlying group
MIPElementAlgebraToTable := function(el, kG, T)
local k, maxweight, i, lim, sup, coeffs, v0, res, exps, nonzeroexps, sum, S, initprod, prod, s, es, innersum, fac, l;
    k := T.fld;
    if Augmentation(el) <> Zero(k) then
        Error("<el> must be an element in the augmentation ideal");
    fi;

    maxweight := Maximum(T.wgs);
    # we find the limit of the group elements we need for the case T does not see the whole group basis
    if maxweight < Maximum(T.pre.jen.weights) then
        for i in [1..Size(T.pre.jen.pcgs)] do
            if T.pre.jen.weights[i] > maxweight then
                lim := i - 1;
            fi;
        od;
    else
        lim := Size(T.pre.jen.pcgs);
    fi; 
    sup := Support(el);
    coeffs := CoefficientsBySupport(el);
    v0 := ListWithIdenticalEntries(T.dim, Zero(k));
    res := ShallowCopy(v0);
    ## for each element g1^e1...gm^em we convert g1^e1...gm^em to an expression in the Jennings basis using the formulas known for this
    for i in [1..Size(sup)] do
        exps := ExponentsOfPcElement(T.pre.jen.pcgs, sup[i]){[1..lim]};
        nonzeroexps := PosNonzero(exps);
        sum := ShallowCopy(v0);
        # run over all supsets of non-zero exponents. Each gives a summand
        for S in Combinations( nonzeroexps ) do
            initprod := false;
            prod := ShallowCopy(v0);
            # this is the factor (g_s - 1)^es which is converted using binomial coefficients
            for s in S do
                es := exps[s];
                innersum := ShallowCopy(v0);
                innersum[T.pre.poswone[s]] := One(k)*Binomial(es, 1);
                fac := ShallowCopy(v0); # corresponding to g_s - 1
                fac[T.pre.poswone[s]] := One(k);
                for l in [2..es] do
                    fac := MultByTable(T, fac, fac); # obtain (g_s - 1)^l
                    innersum := innersum + One(k)*Binomial(es, l)*fac;
                od;
                if initprod then 
                    prod := MultByTable(T, prod, innersum);
                else
                    prod := innersum;
                    initprod := true;
                fi;
            od;
            sum := sum + prod;
        od;
        res := res + coeffs[i]*sum;
    od;
    return res;
end;

### convert from element of table to a representative in the group algebra
MIPElementTableToAlgebra := function(v, T, kG)
local k, res, i, exps, prod, j, G, iota, basisels;
    k := T.fld; 
    res := Zero(kG);
    G := UnderlyingMagma(kG);
    iota := Embedding(G, kG);
    basisels := List(T.pre.jen.pcgs, x -> x^iota - One(kG));
    for i in [1..T.dim] do
        if v[i] <> Zero(k) then
            exps := T.pre.exps[i];
            prod := One(kG);
            for j in [1..Size(exps)] do
                prod := prod*(basisels[j]^exps[j]);
            od;
            prod := v[i]*prod;
            res := res + prod;
        fi;
    od;
    return res;
end;
