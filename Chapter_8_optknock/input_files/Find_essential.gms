*** this program finds out the rxns which are essential to biomass production in the network

$INLINECOM /*  */
$set myroot iAF1260/iAF1260
$set condition aero
$set biom 0.093052
$set reg without


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

*** Reactions shut off during aerobic glucose condition
offanaero(j)
$include "%myroot%_offanaero.txt"

*** Exchange Reactions used for metabolite uptake
exchange(j)
$include "%myroot%_source_M9.txt"

*** a dummy set
dummy(j)
;
***************************************************
Parameters

*** Reaction Type
rxntype(j)
$include "%myroot%_rxntype.txt"

*** Lower/Upper flux bounds
LB(j) , UB(j)

*** min/max flux bounds
high(j)

*** Stoichiometry Matrix
S(i,j)
$include "%myroot%_sij.txt"

c(j)

x

;

c(j) = 1;

***************************************************
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
con1a,con1b

;

Obj..			z =e= v('Ec_biomass_iAF1260_WT_59p81M') ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;
con1a(j)..		v(j)  =g=  LB(j)*c(j) ;
con1b(j)..		v(j)  =l=  UB(j)*c(j) ;

**************************************************
*** Get bounds for Fluxes

*scalar m /1000/;
scalar M /1000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = M;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -M;
UB(j)$(rxntype(j) = 1) = M;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = M;
LB(j)$(exchange(j)) = -1000;

*** setting core biomass flux to 0
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** setting LB wild biomass flux to 10% of theoretical maximum
*LB('Ec_biomass_iAF1260_WT_59p81M') = %biom%;

***shutting off reactions in anaerobic glucose conditions

*LB(j)$(offanaero(j)) = 0;
*UB(j)$(offanaero(j)) = 0;

*LB(j)$(offaeroglucose(j)) = 0;
*UB(j)$(offaeroglucose(j)) = 0;

*** putting ATP Maintenance value

LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

***glucose and oxygen uptake rates

LB('EX_glc(e)') = -10;
LB('EX_o2(e)') = -20;

v.lo(j) = LB(j);
v.up(j) = UB(j);

Model cycles /all/ ;

cycles.optfile = 1;

for(x = 1 to card(j),

	c(j)$(ord(j) ne x) = 1;
	c(j)$(ord(j) eq x) = 0;

	Solve cycles using lp maximizing z;
	If(cycles.modelstat eq 1,
		high(j)$(c(j) =0) = z.l;
	);
);

file ff /Find_essential_%condition%_%reg%.txt/;
put ff;
put "/"//;
put "*** Rxns which are essential"//;
loop(j,
	if(high(j) lt %biom%,
		put "'"j.tl:30:0"'"/;
	);
);
put "/";
putclose ff;
