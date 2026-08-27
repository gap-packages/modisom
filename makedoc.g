#############################################################################
##
##  makedoc.g
##
##  Builds the package documentation with AutoDoc/GAPDoc.
##
#############################################################################

LoadPackage("AutoDoc");

# Run this from the package's root directory: gap makedoc.g
AutoDoc(rec(
    autodoc := rec(scan_dirs := []),
    gapdoc := rec(main := "main", files := []),
    extract_examples := true,
    scaffold := rec(
        includes := [
            "intro.xml",
            "tables.xml",
            "autiso.xml",
            "modisom.xml",
            "nilquo.xml",
            "relfree.xml"
        ],
        entities := rec(
            ModIsom := "<Package>ModIsom</Package>",
        ),
        bib := "modisom.bib",
    ),
));

QuitGap();
