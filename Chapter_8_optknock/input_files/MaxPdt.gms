*** this program finds out the thoretical maximum of pdt that can be achieved

$INLINECOM /*  */
$set myroot iAF1260/iAF1260


options
        limrow = 10000
        limcol = 10000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        sysout = off
        solprint = on;


Sets

***     Metabolites
i
$include "%myroot%_cmp.txt"

***     Reactions
j
$include "%myroot%_rxnnames.txt"

*** Reactions shut off during aerobic glucose condition
offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

*** Exchange Reactions used for metabolite uptake
exchange(j)
$include "%myroot%_source_M9.txt"

;
***************************************************
Parameters

*** Reaction Type
rxntype(j)
$include "%myroot%_rxntype.txt"

*** Lower/Upper flux bounds
LB(j) , UB(j)

*** Stoichiometry Matrix
S(i,j)
$include "%myroot%_sij.txt"

;


****************************************************
Variables

*** Maximizing function
z

*** Reaction Fluxes
v(j)

;
*************************************************** 
Equations

Obj
Stoic   stoichiometric limitation 

;

Obj..			z =e= v('EX_ttdca(e)') ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;


**************************************************
*** Get bounds for Fluxes

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = 1000;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -1000;
UB(j)$(rxntype(j) = 1) = 1000;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = 1000;
LB(j)$(exchange(j)) = -1000;

***shutting off biomass core rxn flux

LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

***shutting off reactions in aerobic glucose conditions

LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** putting ATP Maintenance value

LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** Glucose and oxygen uptake

LB('EX_glc(e)') = -100;
LB('EX_o2(e)') = -200;

v.lo(j) = LB(j);
v.up(j) = UB(j);

Model maxpdt /all/ ;

maxpdt.optfile = 1;

Solve maxpdt using lp maximizing z;

file ff /MaxPdt.txt/;
put ff;
put "Maximum product yield is : " z.l:0:8//;
loop(j,
	put j.tl:0:30,"    ",v.l(j):0:8/;
);
putclose ff;
