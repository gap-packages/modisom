#############################################################################
##
##  Compares two outputs of bench.g, see there.
##
SetPrintFormattingStatus( "*stdout*", false );

BenchRead := function( file )
    local res, line, w;
    res := [];
    for line in SplitString( StringFile( file ), "\n" ) do
        w := Filtered( SplitString( line, " " ), x -> x <> "" );
        if Length( w ) = 4 then
            Add( res, rec( name := w[1], time := Int( w[2] ), size := w[3],
                           cano := w[4] ) );
        fi;
    od;
    return res;
end;

BenchRatio := function( x, y )
    local q;
    if x = 0 then return ""; fi;
    q := String( QuoInt( 100 * y, x ) + 1000 );
    return Concatenation( String( Int( q{[1..Length( q )-2]} ) - 10 ), ".",
                          q{[Length( q )-1, Length( q )]} );
end;

a := BenchRead( before );;
b := BenchRead( after );;
Print( String( "case", -16 ), String( "before", 10 ), String( "after", 10 ),
       String( "ratio", 8 ), "\n" );
for x in a do
    y := First( b, z -> z.name = x.name );
    if y = fail then continue; fi;
    Print( String( x.name, -16 ), String( x.time, 10 ), String( y.time, 10 ),
           String( BenchRatio( x.time, y.time ), 8 ) );
    if x.size <> y.size then Print( "  ORDER DIFFERS" ); fi;
    if x.cano <> y.cano then Print( "  canonical form differs" ); fi;
    Print( "\n" );
od;
QUIT;
