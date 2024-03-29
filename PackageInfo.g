#############################################################################
##  
##  PackageInfo.g for the package `modisom'                      Bettina Eick
##  
SetPackageInfo( rec(
PackageName := "ModIsom",
Subtitle := "Computing automorphisms and checking isomorphisms for modular group algebras of finite p-groups",
Version := "2.5.4",
Date := "27/02/2023", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec( 
    LastName      := "Eick",
    FirstNames    := "Bettina",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "beick@tu-bs.de",
    WWWHome       := "http://www.iaa.tu-bs.de/beick",
    PostalAddress := Concatenation(
               "Institut Analysis und Algebra\n",
               "TU Braunschweig\n",
               "Universitätsplatz 2\n",
               "D-38106 Braunschweig\n",
               "Germany" ),
    Place         := "Braunschweig",
    Institution   := "TU Braunschweig"
  ),
  rec(
    LastName      := "Konovalov",
    FirstNames    := "Olexandr",
    IsAuthor      := false,
    IsMaintainer  := true,
    Email         := "obk1@st-andrews.ac.uk",
    WWWHome       := "https://olexandr-konovalov.github.io/",
    PostalAddress := Concatenation( [
                     "School of Computer Science\n",
                     "University of St Andrews\n",
                     "Jack Cole Building, North Haugh,\n",
                     "St Andrews, Fife, KY16 9SX, Scotland" ] ),
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
  ) ],

Status := "accepted",
CommunicatedBy := "Olexandr Konovalov (St Andrews)",
AcceptDate := "11/2013",

PackageWWWHome  := "https://gap-packages.github.io/modisom/",
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/modisom",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/modisom-", ~.Version ),
ArchiveFormats := ".tar.gz",

AbstractHTML := 
  "The <span class=\"pkgname\">ModIsom</span> package contains various methods for computing with nilpotent associative algebras. In particular, it contains a method to determine the automorphism group and to test isomorphis of such algebras over finite fields and of modular group algebras of finite p-groups, and it contains a nilpotent quotient algorithm for finitely presented associative algebras and a method to determine Kurosh algebras.",

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
  GAP := ">=4.7",
  NeededOtherPackages := [["Polycyclic", ">=1.0"]], 
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

BannerString := Concatenation( "Loading ModIsom ", ~.Version, "... \n"),
AvailabilityTest := ReturnTrue,
Autoload := false,
TestFile := "tst/testall.g",
Keywords := ["modular isomorphism problem",
             "automorphism group", 
             "isomorphism testing",
             "nilpotent algebras",
             "nilpotent quotient",
             "Kurosh algebras"]

));
