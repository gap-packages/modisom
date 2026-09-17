#############################################################################
##
#F VectorCanonicalForm( pcgs, v, F, l, base )
##
## Pcgs consists of 2-tuples, component 2 acts per right multiplication.
## Computes modulo base{[l+1..n]} or mod [] if base=fail
##
BindGlobal( "VectorCanonicalForm", function( pcgs, v, F, l, base )
    local p, f, d, o, stab, tran, cano, indu, tail, B, project, projMat, zero,
          echelon, echelonCoeffs, pivots, residue, residueCoeffs, pivot,
          coeff, lead, i, j, k, e, ec, w, wc, b, s, t; 

    # the trivial case is not supported
    if Length(pcgs) = 0 then return fail; fi;

    # set up
    p := Characteristic(F);
    f := DegreeOverPrimeField(F);
    d := Length(v);
    o := IdentityMat(d,F);

    # get a basis of F over its prime field
    B := Basis(F);
    zero := Zero( GF(p) );

    # the projection onto the first l coordinates modulo base, as a
    # matrix: it is applied to every tail in every round below
    projMat := IndVectorMatrix( l, base, F );
    if projMat = fail then
        project := g -> IndVector( g, l, base );
    else
        project := g -> g * projMat;
    fi;

    # init
    stab := ShallowCopy(pcgs);
    tran := pcgs[1]^0;
    cano := ShallowCopy(v);
    indu := project( cano );

    # get tails
    tail := List( stab, x -> project(cano*(x[2] - o)));

    # use induction on natural flag
    for i in [2..l] do

        # catch relevant entry
        e := List( tail, x -> x[i] );
        ec := List(e, x -> Coefficients(B,x));
        w := indu[i];
        wc := Coefficients(B,w);

        # compute stabilizer: sift the entries into an echelon basis,
        # keeping each basis vector as a combination of the entries
        # chosen so far.  SolutionMat would echelonise the chosen ones
        # again for every entry.
        #   b              the chosen elements, as indices into stab and ec
        #   echelon        the echelon basis, with leading entries 1
        #   pivots         the column of each leading entry
        #   echelonCoeffs  each basis vector in terms of the chosen entries
        b := [];
        echelon := [];
        pivots := [];
        echelonCoeffs := [];
        for j in Reversed([1..Length(e)]) do

            # residue = ec[j] - residueCoeffs * (the chosen entries)
            residue := ec[j];
            residueCoeffs := ListWithIdenticalEntries( Length(b), zero );
            for k in [1..Length(echelon)] do
                coeff := residue[pivots[k]];
                if coeff <> zero then
                    residue := residue - coeff * echelon[k];
                    residueCoeffs := residueCoeffs
                                     + coeff * echelonCoeffs[k];
                fi;
            od;
            pivot := PositionNonZero( residue );

            if pivot > Length( residue ) then

                # the entry lies in the span of the chosen ones
                s := List( residueCoeffs, IntFFE );
                for k in Reversed([1..Length(s)]) do
                    if s[k]<>0 then 
                        stab[j] := stab[j]*stab[b[k]]^(-s[k] mod p);
                    fi;
                od;
            else

                # the entry enlarges the span
                lead := residue[pivot];
                Add( b, j );
                Add( echelon, residue / lead );
                Add( pivots, pivot );
                echelonCoeffs := List( echelonCoeffs,
                                       x -> Concatenation( x, [zero] ) );
                Add( echelonCoeffs, Concatenation(
                         List( residueCoeffs, x -> -x/lead ), [lead^-1] ) );
            fi;
        od;

        # compute minimal element
        t := CoeffsMinimalElement(wc, ec{b});

        # get transversal element
        for k in Reversed([1..Length(t)]) do
            if t[k]<>0 then 
                tran := tran * stab[b[k]]^t[k];
            fi;
        od;

        # set up for next round
        if t <> 0*t then
            cano := v * tran[2];
            indu := project( cano );
        fi;
        if Length(b)>0 then 
            stab := stab{Difference([1..Length(e)], b)};
            tail := List( stab, x -> project(cano*(x[2] - o)));
        fi;
    od;

    if MIP_CHECK_CNF then 
        if ForAny( stab, x -> 
            IndVector(cano*x[2],l,base) <> IndVector(cano,l,base) ) then 
            Error("stabilizer does not stabilize in vector cano form");
        fi;
    fi;

    return rec( cano := cano, stab := stab, tran := tran );
end );


#############################################################################
##
#F SubspaceCanonicalForm( pcgs, id, base, F )
##
## Assumes that base is echelonised.
## Pcgs consists of 2-tuples, component 2 acts per right multiplication.
##
BindGlobal( "SubspaceCanonicalForm", function( pcgs, id, base, F )
    local d, l, I, stab, cano, tran, n, c, b, f;

    # a preliminary check
    if Length(pcgs) = 0 or Length( base ) = 0 then 
        return rec( cano := base, stab := pcgs, tran := id );
    fi;

    # set up
    d := Length(base[1]);
    l := Length(base);
    I := IdentityMat(d,F);

    # use induction on length of base
    for n in Reversed([1..l]) do

        # the first basis vector 
        if n = l then 

            # get cano form
            c := VectorCanonicalForm( pcgs, base[n], F, d, fail );
            #Print("    ag-stab length: ",Length(pcgs)-Length(c.stab),"\n");

            stab := c.stab;
            tran := c.tran;
            cano := MyTriangulizedBaseMat( base * tran[2] );

        # the others 
        elif Length(stab) > 0 then 

            # get basis through cano{[n+1..l]}
            f := BaseSteinitzVectors( I, cano{[n+1..l]} ).factorspace;
            b := Concatenation( f, cano{[n+1..l]});

            # get cano form
            c := VectorCanonicalForm( stab, cano[n], F, d-l+n, b );
            #Print("    ag-stab length: ",Length(stab)-Length(c.stab),"\n");

            # translate result
            stab := c.stab;
            tran := tran * c.tran;
            cano := MyTriangulizedBaseMat( base * tran[2] );
        fi;

    od;

    if MIP_CHECK_CNF then 
        if not ForAll( stab, x -> cano = MyTriangulizedBaseMat(cano*x[2]) ) then 
            Error("stabilizer does not stabilize in subspace cano form");
        fi;
    fi;

    return rec( cano := cano, stab := stab, tran := tran );
end );
        
