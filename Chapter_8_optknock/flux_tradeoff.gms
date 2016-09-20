*** This code finds the tradeoff between product and biomass

****************************************************
$INLINECOM /*  */
$set myroot input_files/iAF1260/iAF1260
$set myroot2 input_files/
$set pdt ac
$set value 154.285
$set max_val 171.428
$set biom 0.3037459

Sets

*** compounds
i
$include "%myroot%_cmp.txt"

*** rxns
j
$include "%myroot%_rxnnames.txt"

*** Reactions shut off during anerobic glucose condition
anaerobic(j)
$include "%myroot%_offanaero.txt"

*** Reactions shut off during aerobic glucose condition
offaeroglu(j)
$include "%myroot%_offaeroglu.txt"

*** rxns providing nutrients
exchange(j)
$include "%myroot%_source_M9.txt"

*** constraints
cons(j)
/
'EX_glc(e)'
'ATPM'
/

*** objective
obj(j) / 'Ec_biomass_iAF1260_WT_59p81M' /

*** biomass
bio(j) / 'Ec_biomass_iAF1260_WT_59p81M' /

;

****************************************************

Parameters

*** stoichiometric matrix
s(i,j)
$include "%myroot%_sij.txt"

*** rxn type
rxntype(j)
$include "%myroot%_rxntype.txt"

*** flux bounds
UB(j)
LB(j)
;

****************************************************

Variables

*** objective
z , zd

*** rxn fluxes
v(j)

*** dual for s(i,j)
lamda(i)

*** dual for constraint
mu(j)

negative variable

*** dual for biomass
phi(j)

*** dual for LB(j)
deltal(j)

positive variable

*** dual for UB(j)
deltau(j)
;

Parameter c /0/;
****************************************************

Equations

zprimal
primal1
primal2
;

*** primal

zprimal..               z =e= v('EX_%pdt%(e)');
*zprimal..               z =e= v('Ec_biomass_iAF1260_WT_59p81M');
primal1(i)..            sum(j, s(i,j) * v(j)) =e= 0;
*primal2..               v('Ec_biomass_iAF1260_WT_59p81M') =e= .884419*c*0.1;
primal2..               v('Ec_biomass_iAF1260_WT_59p81M') =e= .22007*c*0.1;
*primal2..               v('Ec_biomass_iAF1260_WT_59p81M') =e= .317949*c*0.1;

****************************************************

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
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** setting LB wild biomass flux to 10% of theoretical maximum
*LB('Ec_biomass_iAF1260_WT_59p81M') = 0.09305;

*** shutting off reactions in aerobic glucose conditions

*LB(j)$anaerobic(j) = 0;
*UB(j)$anaerobic(j) = 0;

LB(j)$offaeroglu(j) = 0;
UB(j)$offaeroglu(j) = 0;

LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** uptake

LB('Ex_glc(e)') = -10;
LB('EX_o2(e)') = -20;

*** bounds

*LB('ALCD2x_f') = 0;
*UB('ALCD2x_f') = 0;

LB('TPI_f') = 0;
UB('TPI_f') = 0;

*LB('PSP_L') = 0;
*UB('PSP_L') = 0;

*LB('GHMT2r_f') = 0;
*UB('GHMT2r_f') = 0;

*LB('GHMT2r_f') = 0;
*UB('GHMT2r_f') = 0;

LB('ATPS4rpp_f') = 0;
UB('ATPS4rpp_f') = 0;

*LB('ENO_f') = 0;
*UB('ENO_f') = 0;

*LB('SUCDi') = 0;
*UB('SUCDi') = 0;

*LB('PGI_f') = 0;
*UB('PGI_f') = 0;

*LB('MDH_f') = 0;
*UB('MDH_f') = 0;

*LB('GLDBRAN2') = 0;
*UB('GLDBRAN2') = 0;

*LB('LPLIPAL2E160') = 0;
*UB('LPLIPAL2E160') = 0;


v.lo(j) = LB(j);
v.up(j) = UB(j);
****************************************************

model primal
/
zprimal
primal1
primal2
/;

primal.optfile = 1;
primal.holdfixed = 1;
scalar yield;


*solve primal using lp maximizing z;
*solve primal using lp minimizing z;


file ff /flux_tradeoff.csv/;
ff.pc=5;
put ff;

for(c = 0 to 10 by 1,
    solve primal using lp maximizing z;
    put z.l,;
    solve primal using lp minimizing z;
    put z.l:0:8,v.l('Ec_biomass_iAF1260_WT_59p81M'):0:8/;
);


/*
c = 10;
solve primal using lp minimizing z;
If(primal.modelstat eq 1,
	put "/"//;
	loop(j$(rxntype(j) ne 2),
		put "'"j.tl:0:30"'  " v.l(j):0:8/;
	);
	put /"/"
);

putclose ff;
*/
