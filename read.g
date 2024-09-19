#############################################################################
##
#W    read.g                 The ModIsom package                 Bettina Eick
##

#############################################################################
##
#R  Global vars
##
SMTX.RAND_ELM_LIMIT := 5000;
MIP_GLLIMIT := infinity;
COVER_LIMIT := 100;
POWER_LIMIT := 1000;

# some flags for the algorithm
if not IsBound(USE_PARTI) then USE_PARTI := false; fi;
if not IsBound(USE_CHARS) then USE_CHARS := false; fi;
if not IsBound(USE_MSERS) then USE_MSERS := false; fi;

# checking modes for the package
if not IsBound(CHECK_AUT) then CHECK_AUT := false; fi;
if not IsBound(CHECK_STB) then CHECK_STB := false; fi;
if not IsBound(CHECK_CNF) then CHECK_CNF := false; fi;
if not IsBound(CHECK_NQA) then CHECK_NQA := false; fi;

# store info
if not IsBound(STORE) then STORE := true; fi;
if not IsBound(COVER) then COVER := true; fi;
if not IsBound(ALLOW) then ALLOW := true; fi;

#############################################################################
##
#R  Read the install files.
##
ReadPackage("modisom", "gap/cfstab/general.gi");
ReadPackage("modisom", "gap/cfstab/pgroup.gi");
ReadPackage("modisom", "gap/cfstab/orbstab.gi");
ReadPackage("modisom", "gap/cfstab/check.gi");

ReadPackage("modisom", "gap/tables/sparse.gi");
ReadPackage("modisom", "gap/tables/linalg.gi");
ReadPackage("modisom", "gap/tables/weight.gi");
ReadPackage("modisom", "gap/tables/tables.gi");
ReadPackage("modisom", "gap/tables/genalg.gi");
ReadPackage("modisom", "gap/tables/basic.gi");
ReadPackage("modisom", "gap/tables/quots.gi");
ReadPackage("modisom", "gap/tables/cover.gi");
ReadPackage("modisom", "gap/tables/isom.gi");

ReadPackage("modisom", "gap/nilalg/basic.gi");
ReadPackage("modisom", "gap/nilalg/liecl.gi");

ReadPackage("modisom", "gap/autiso/chains.gi");
ReadPackage("modisom", "gap/autiso/fprint.gi");
ReadPackage("modisom", "gap/autiso/tstepc.gi");
ReadPackage("modisom", "gap/autiso/initial.gi");

ReadPackage("modisom", "gap/autiso/check.gi");
ReadPackage("modisom", "gap/autiso/induc.gi");
ReadPackage("modisom", "gap/autiso/autiso.gi");

ReadPackage("modisom", "gap/grpalg/basis.gi");
ReadPackage("modisom", "gap/grpalg/tables.gi");
ReadPackage("modisom", "gap/grpalg/autiso.gi");
ReadPackage("modisom", "gap/grpalg/check.gi");
ReadPackage("modisom", "gap/grpalg/head.gi");

ReadPackage("modisom", "gap/grpalg/jenningsBounds.gi");
ReadPackage("modisom", "gap/grpalg/detbins.gi");
ReadPackage("modisom", "gap/grpalg/detbinsRT.gi");
ReadPackage("modisom", "gap/grpalg/collect.gi");
ReadPackage("modisom", "gap/grpalg/chkbins.gi");
ReadPackage("modisom", "gap/grpalg/kernelsize.gi");
ReadPackage("modisom", "gap/grpalg/tabletoalgebraandback.gi");
ReadPackage("modisom", "gap/grpalg/jenningsConjecture.gi");

ReadPackage("modisom", "gap/algnq/exam.gi");
ReadPackage("modisom", "gap/algnq/algnq.gi");

ReadPackage("modisom", "gap/rfree/kurosh.gi");
ReadPackage("modisom", "gap/rfree/polyid.gi");
ReadPackage("modisom", "gap/rfree/engel.gi");

#ReadPackage("modisom", "gap/group/group.gi");

ReadPackage("modisom", "gap/pilib/kur_3_3_Q.gi");
ReadPackage("modisom", "gap/pilib/kur_3_3_3.gi");
ReadPackage("modisom", "gap/pilib/kur_3_3_2.gi");
ReadPackage("modisom", "gap/pilib/kur_3_3_4.gi");
ReadPackage("modisom", "gap/pilib/kur_4_3_Q.gi");
ReadPackage("modisom", "gap/pilib/kur_4_3_3.gi");
ReadPackage("modisom", "gap/pilib/kur_4_3_2.gi");
ReadPackage("modisom", "gap/pilib/kur_4_3_4.gi");
ReadPackage("modisom", "gap/pilib/kur_2_4_Q.gi");
ReadPackage("modisom", "gap/pilib/kur_2_4_3.gi");
ReadPackage("modisom", "gap/pilib/kur_2_4_9.gi");
ReadPackage("modisom", "gap/pilib/kur_2_4_2.gi");
ReadPackage("modisom", "gap/pilib/kur_2_4_4.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_Q.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_5.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_3.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_9.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_2.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_4.gi");
ReadPackage("modisom", "gap/pilib/kur_2_5_8.gi");
ReadPackage("modisom", "gap/pilib/readlib.gi");

