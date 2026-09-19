#############################################################################
##
##  Benchmarks for the canonical form of a modular group algebra.
##
##  From the package directory:
##
##      gap -q --packagedirs . dev/benchmarks/bench.g > before.txt
##      ... change the code ...
##      gap -q --packagedirs . dev/benchmarks/bench.g > after.txt
##      gap -q -c 'before := "before.txt";; after := "after.txt";;' \
##          dev/benchmarks/compare.g
##
##  A line of output is a name, the CPU time in ms, the order of the
##  automorphism group found and a checksum of the canonical table.  The
##  order must not change.  The checksum may, if the canonical form is
##  chosen in another way, but then it changes for every case.
##
##  Add  -c 'BENCH_SLOW := true;;'  for the cases that take minutes.
##
LoadPackage( "modisom" );
SetInfoLevel( InfoModIsom, 0 );
SetPrintFormattingStatus( "*stdout*", false );   # no line breaks
if not IsBound( BENCH_SLOW ) then BENCH_SLOW := false; fi;

# name, prime, order and id of the group, and the power of the radical
# that is factored out (0 for none)
BENCH_CASES := [
  [ "32-49",        2,  32,     49, 0 ],
  [ "64-200",       2,  64,    200, 0 ],
  [ "64-220",       2,  64,    220, 0 ],
  [ "64-240",       2,  64,    240, 0 ],
  [ "64-264",       2,  64,    264, 0 ],
  [ "64-266",       2,  64,    266, 0 ],
  [ "27-3",         3,  27,      3, 0 ],
  [ "81-7",         3,  81,      7, 0 ],
  [ "256-26973/5",  2, 256,  26973, 5 ],
  [ "512-406945/4", 2, 512, 406945, 4 ],
];

BENCH_SLOW_CASES := [
  [ "256-26973/6",  2, 256,  26973, 6 ],
  [ "128-2326",     2, 128,   2326, 0 ],
  [ "256-542",      2, 256,    542, 0 ],
  [ "256-6732",     2, 256,   6732, 0 ],
];

if BENCH_SLOW then Append( BENCH_CASES, BENCH_SLOW_CASES ); fi;

BenchRun := function( case )
    local FG, T, t, R, tab;
    FG := GroupRing( GF( case[2] ), SmallGroup( case[3], case[4] ) );
    if case[5] = 0 then
        T := NilpotentTableOfRad( FG );
    else
        T := TableOfRadQuotient( FG, case[5] );
    fi;
    Reset( GlobalMersenneTwister, 1 );
    Reset( GlobalRandomSource, 1 );
    t := Runtime();
    R := CanoFormWithAutGroupOfTable( T );
    t := Runtime() - t;

    # the entries as integers: their printed form depends on how the
    # vectors are stored
    tab := List( R.cano.tab, x -> List( x, y -> List( y, IntFFE ) ) );
    Print( case[1], " ", t, " ", R.auto.size, " ",
           HexSHA256( String( tab ) ){[1..16]}, "\n" );
end;

for case in BENCH_CASES do BenchRun( case ); od;
QUIT;
