# CHANGES to the 'ModIsom' GAP package

## 3.2.0 (unreleased)

 - prefixed the global variables defined in `read.g` with `MIP_`, e.g.
   `COVER_LIMIT` is now `MIP_COVER_LIMIT`
 - `IsCoveredByTheory` and `IsCoveredByTheoryAllFields` now do what their
   names suggest

## 3.1.0 (2026-01-09)

 - fixed `JenningsInfoAllFields`, which used an invariant valid only over
   `GF(p)`, namely the quotient by the second dimension subgroup
   (Passi-Sehgal 72); this made the results of versions 3.0.0 and 3.0.1
   mathematically wrong
 - `JenningsBoundPairwise` now does what it always claimed to do; the old
   function is now called `JenningsBoundPairwiseWithIds`
 - fixed a bug in `MIPBinsByGTInternal` and `MIPBinsByGTAllFieldsInternal`
   which could cause break loops, as invariants were applied to empty bins
 - added the invariant `G/Agemo(G')` for `G'` cyclic and `G` 2-generated
   (Margolis-Sakurai 25) to `CyclicDerivedInfo` and
   `CyclicDerivedInfoAllFields`
 - `IsCoveredByTheoryAllFields` now covers metacyclic groups
   (Margolis-Sakurai 25, Theorem 3.6)
 - `IsCoveredByTheory` now covers Roehl 90, Theorem 3.1.2, using
   Minac-Rogelstad-Tan 2016, Proposition 3.4 for the dimensions of
   dimension subgroup quotients of free pro-p groups
 - `IsCoveredByTheory` and `SandlingInfo` now cover Brenner-Garcia-Lucas 24,
   Theorems A and B
 - `MIPSplitGroupsByAlgebras` now returns 1 instead of 0 when the first
   layers of the group algebras are not isomorphic, for consistency with
   its other return values
 - changed the output of `NormalSubgroupsInfo` so that one can tell which
   entry belongs to which subgroup; the invariant itself is unchanged
 - fixed several invariant functions which raised errors for groups they
   were not written for, such as abelian groups
 - extended the test files
 - minor manual changes

## 3.0.1 (2024-12-28)

 - made LAGUNA a required package, which avoids warnings when ModIsom is
   loaded

## 3.0.0 (2024-09-23)

This release merges the functionality of ModIsomExt, a package by Leo
Margolis and Tobias Moede, into ModIsom; see
<https://doi.org/10.1016/j.jaca.2022.100001>. ModIsomExt installed its
functions as methods; here they are all `BindGlobal`.

New group-theoretic invariants and splitting functions (`detbins.gi`):

 - `RefineBins`, which refines lists of groups and not only of ids,
   `ConjugacyClassInfo` (with the parameter of Parmenter-Polcino
   Milies), `JenningsInfo` (with Hertweck's `G/D4(G)` and better id
   recognition), `SandlingInfo` (with Baginski and Margolis-Moede on
   small group rings), `BaginskiInfo` (Baginski 99), `CenterDerivedInfo`,
   `FrattiniInfo`, `JenningsDerivedInfo`, `BaginskiCarantiInfo`,
   `NilpotencyClassInfo` and `DimensionTwoCohomology`, all from ModIsomExt
 - `JenningsInfoAllFields`, using only subsequent Jennings quotients
 - `DimensionSecondHochschild`, the dimension of the second Hochschild
   cohomology group
 - `IsCoveredByTheory` extended by various new results, and the new
   `IsCoveredByTheoryAllFields`
 - `Theorem41MS22` for the result of Margolis-Stanojkovski
 - `MaximalAbelianDirectFactor` and `MaximalElementaryAbelianDirectFactor`
   for cancelling factors
 - `AgemoInvariantAllM`, `OmegaInvariantAllM`, `AgemoCenterInvariantAllM`
   and `NormalSubgroupsInfo`, based on Garcia-Lucas and
   Margolis-Sakurai-Stanojkovski
 - `CyclicDerivedInfo` and `CyclicDerivedInfoAllFields`, based on
   Garcia-Lucas-Del Rio-Stanojkovski and Garcia-Lucas-Del Rio
 - `MIPBinsByGTInternal` extended by the new invariants and an option to
   turn off cohomology computations; `BinsByGT` now takes 2-4 arguments
   and uses the new invariants
 - `MIPBinsByGTAllFieldsInternal` and `BinsByGTAllFields` for all fields
 - `MIPSplitGroupsByGroupTheoreticalInvariants` and its
   `NoCohomology`, `AllFields` and `AllFieldsNoCohomology` variants, which
   take just a list of groups

Algebra side (`chkbins.gi`, formerly `checkbin.gi`):

 - `MIPBinSplit` from ModIsomExt, with a fix of the main function and a
   change allowing the field to be increased
 - `MIPSplitGroupsByAlgebras`, which takes a list of groups and optionally
   a number by which to increase the field
 - reports on time and similar data, as in ModIsomExt

New files:

 - `collect.gi` from ModIsomExt: computing with Jennings bases in Loewy
   series quotients of the augmentation ideal, plus the new function
   `ModIsomTable`
 - `jenningsBounds.gi` and `jenningsConjecture.gi`: the Jennings bound
   functions, which in ModIsomExt lived in `detbins.gi`
 - `kernelsize.gi`: kernel sizes of power maps
 - `tabletoalgebraandback.gi`: `MIPElementTableToAlgebra` and
   `MIPElementAlgebraToTable`, converting elements between the group
   algebra and the table
 - `detbinsRT.gi`: the ring-theoretic functions, from the second half of
   the old `checkbin.gi`, now `chkbins.gi`

Other changes:

 - fixed `PowerBasisWeights` in `fprint.gi`, as already done in ModIsomExt
 - documented the new functions, the changes to existing ones, and a number
   of functions that were undocumented before

## 2.5.4 (2023-02-27)

## 2.5.3 (2022-08-09)

## 2.5.2 (2022-04-29)

## 2.5.1 (2020-01-28)

## 2.5.0 (2019-02-19)

## 2.4.0 (2018-09-16)
