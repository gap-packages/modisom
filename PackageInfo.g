#############################################################################
##  
##  PackageInfo.g for the package `modisom'                      Bettina Eick
##  
SetPackageInfo( rec(
PackageName := "ModIsom",
Subtitle := "Computing automorphisms and checking isomorphisms for modular group algebras of finite p-groups",
Version := "2.0",
Date := "01/09/2011",

Persons := [
  rec( 
    LastName      := "Eick",
    FirstNames    := "Bettina",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "beick@tu-bs.de",
    WWWHome       := "http://www.icm.tu-bs.de/~beick",
    PostalAddress := Concatenation( [
                       "Institut Computational Mathematics\n",
                       "Pockelsstrasse 14, 38106 Braunschweig\n",
                       "Germany" ] ),
    Place         := "Braunschweig",
    Institution   := "TU Braunschweig"
  ) ],

Status := "dev",

PackageWWWHome := "http://www-public.tu-bs.de:8080/~beick/soft/modisom/",

ArchiveFormats := ".tar.gz",
ArchiveURL     := Concatenation( ~.PackageWWWHome, "modisom-", ~.Version ),
README_URL     := Concatenation( ~.PackageWWWHome, "README" ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),

AbstractHTML := 
  "The <span class=\"pkgname\">ModIsom</span> package contains various methods for computing with nilpotent associative algebras. In particular, it contains a method to determine the automorphism group and to test isomorphis of such algebras over finite fields and of modular group algebras of finite p-groups, it contains a nilpotent quotient algorithm for finitely presented associative algebras and a method to determine Kurosh algebras.",

PackageDoc := rec(
  BookName  := "ModIsom",
  ArchiveURLSubset := ["doc", "htm"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing with nilpotent associative algebras",
  Autoload  := true
),


Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [["Polycyclic", ">=1.0"]], 
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

BannerString := "Loading ModIsom 2.0... \n",
AvailabilityTest := ReturnTrue,
Autoload := false,
Keywords := ["modular isomorphism problem",
             "automorphism group", 
             "isomorphism testing",
             "nilpotent algebras",
             "nilpotent quotient",
             "Kurosh algebras"]

));
