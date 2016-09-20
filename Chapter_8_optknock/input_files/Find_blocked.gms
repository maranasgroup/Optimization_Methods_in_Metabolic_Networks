*** this program finds out the rxns which are blocked in the network

$INLINECOM /*  */
$set myroot iAF1260/iAF1260
$set condition aero
$set biom 0.3037
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
low(j) , high(j)

*** Stoichiometry Matrix
S(i,j)
$include "%myroot%_sij.txt"

x

;

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
;

Obj..			z =e= sum(j$dummy(j),v(j)) ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;

**************************************************
*** Get bounds for Fluxes

*scalar m /100/;
scalar M /100000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = M;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -M;
UB(j)$(rxntype(j) = 1) = M;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = M;
LB(j)$(exchange(j)) = -100;

***shutting off biomass core rxn flux

LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

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

	dummy(j)$(ord(j) eq x) = yes;
	dummy(j)$(ord(j) ne x) = no;

	Solve cycles using lp maximizing z;
	If(cycles.modelstat eq 1,
		high(dummy) = z.l;
	);

	Solve cycles using lp minimizing z;
        If(cycles.modelstat eq 1,
		low(dummy) = z.l;
	);
);

scalar epsilon /0.000001/;

file ff /Find_blocked_%condition%_%reg%.txt/;
put ff;
put "/"//;
put "*** Rxns which are blocked in the model"//;
loop(j,
	if(high(j) lt epsilon and low(j) gt -epsilon,
		put "'"j.tl:30:0"'"/;
	);
);
put "/";
putclose ff;
