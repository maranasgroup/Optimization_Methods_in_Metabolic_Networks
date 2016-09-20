*************** OptForce for FA synthesis in E. coli ****************
*** OptForce set for all chain lengths
****************************************************************************

$INLINECOM /*  */
$set myroot1 input_files/iAF1260/iAF1260
$set myroot2 input_files/
$set pdt ac
$set value 154.285
$set max_val 171.428
$set biom 0.09305
$set condition aero
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

*******************************************
Sets

*** Metabolites
i
$include "%myroot1%_cmp.txt"

*** Reactions
j
$include "%myroot1%_rxnnames.txt"

*** fluxes providing nutrients
exchange(j)
$include "%myroot1%_source_M9.txt"

*** Reactions shut off during anaerobic glucose condition
anaerobic(j)
$include "%myroot1%_offanaero.txt"

*** Reactions shut off during anaerobic glucose condition
offaeroglucose(j)
$include "%myroot1%_offaeroglu.txt"

*** Rxns encoded by invivo essential genes
essential(j)
$include "%myroot2%Find_essential_%condition%_%reg%.txt"

*** Product
biom(j) / 'Ec_biomass_iAF1260_WT_59p81M' /


*** rxns with no gene associated
nogene(j)
$include "%myroot1%_nogene_rxns.txt"

*** Exchange rxns
excg(j)
$include "%myroot1%_ex_pp_rxns.txt"

$ONEMPTY
*** Blocked rxns
blocked(j) 
$include "%myroot2%Find_blocked_%condition%_%reg%.txt"

intervention(j)

index /1*1000/
;

*****************************************************************************************

Parameters

*** stoichiometry matrix
s(i,j)
$include "%myroot1%_sij.txt"

*** rxn type
rxntype(j)
$include "%myroot1%_rxntype.txt"

*** flux bounds
LB(j) , UB(j)


*** matrix to store Interventions
matrix(index,j)
;

matrix(index,j) = 0;


********************************
Variables

*** objectives
z
zprimal
zdual

*** flux values
v(j)

*** duals for Sij
lambda(i)

*** for linearizing y(j)*delta(j)
w1l(j), w1u(j)

;

positive variables

*** duals for bounds
delta1l(j), delta1u(j)
;

binary variables
y0(j)
;

**************** INITIALIZING PARAMETRIC VARIABLES AND SETS ****************************

*** Scalars

scalar counter /1/;
scalar M /2000/;
scalar c /0.0001/;
scalar k ;
scalar n /0/;

****************************************************************************************

Equations

primal1
primal7
primal8
primal

dual
dual1
dual2

outer
outer1

outer3
outer4
outer5
outer6
outer7
outer8
outer9
outer10
outer11

outer24
;

*** outer

outer..			z =e= v('EX_%pdt%(e)');
outer1..		sum(j, y0(j)) =e= k;
outer3..		v('Ec_biomass_iAF1260_WT_59p81M') =e= sum(j, delta1u(j)*UB(j) - w1u(j)*UB(j)) - sum(j, delta1l(j)*LB(j) - w1l(j)*LB(j));

outer24(index)..	sum(j, matrix(index,j) * y0(j)) =l= k-1;

outer4(j)..		w1l(j) =l= M*y0(j);
outer5(j)..		w1l(j) =g= -M*y0(j);
outer6(j)..        	w1l(j) =l= delta1l(j) + M*(1-y0(j) );
outer7(j)..		w1l(j) =g= delta1l(j) - M*(1-y0(j) );

outer8(j)..             w1u(j) =l= M*y0(j);
outer9(j)..             w1u(j) =g= -M*y0(j);
outer10(j)..            w1u(j) =l= delta1u(j) + M*(1-y0(j) );
outer11(j)..            w1u(j) =g= delta1u(j) - M*(1-y0(j) );


*** primal

primal..		zprimal =e= v('Ec_biomass_iAF1260_WT_59p81M');
primal1(i)..		sum(j, (S(i,j)*v(j))) =e= 0;
primal7(j)..		v(j) =g= LB(j)*(1-y0(j) );
primal8(j)..		v(j) =l= UB(j)*(1-y0(j) );

*** dual

dual..			zdual =e= sum(j, delta1u(j)*UB(j) - delta1u(j)*y0(j)*UB(j)) - sum(j, delta1l(j)*LB(j) - delta1l(j)*y0(j)*LB(j));

dual1(j)$biom(j)..		sum(i,lambda(i)*S(i,j)) + delta1u(j) - delta1l(j)  =e= 1;
dual2(j)$(not biom(j))..	sum(i,lambda(i)*S(i,j)) + delta1u(j) - delta1l(j)  =e= 0;

***************************************************************
*** Get bounds for Fluxes

scalar max /1000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = max;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -max;
UB(j)$(rxntype(j) = 1) = max;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = max;
LB(j)$(exchange(j)) = -max;

*** setting core biomass flux to 0
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** setting LB wild biomass flux to 10% of theoretical maximum
LB('Ec_biomass_iAF1260_WT_59p81M') = %biom%;

***shutting off reactions in aerobic glucose conditions

*LB(j)$(offanaero(j)) = 0;
*UB(j)$(offanaero(j)) = 0;

*LB(j)$(offaeroglucose(j)) = 0;
*UB(j)$(offaeroglucose(j)) = 0;

*** putting ATP Maintenance value

LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** Glucose and oxygen uptake

LB('Ex_glc(e)') = -10;
LB('EX_o2(e)') = -20;

*** fix binaries
y0.fx(j)$essential(j) = 0;
y0.fx(j)$blocked(j) = 0;
y0.fx(j)$exchange(j) = 0;
y0.fx(j)$nogene(j) = 0;

y0.fx('ATPM') = 0;
y0.fx('Ec_biomass_iAF1260_WT_59p81M') = 0;

*** Fixing certain interventions

y0.fx('PGCD') = 0;
y0.fx('PSERT') = 0;
y0.fx('ETOHtex_f') = 0;

y0.fx('H2Otpp_f') = 0;
y0.fx('H2Otex_f') = 0;
y0.fx('O2tex_f') = 0;
y0.fx('O2tpp_f') = 0;
y0.fx('CO2tex_f') = 0;

********************************

model bilevel
/
primal1
primal7
primal8

dual1
dual2

outer
outer1

outer3
outer4
outer5
outer6
outer7
outer8
outer9
outer10

outer24
/
;

options iterlim = 1000000;
bilevel.optfile = 1;

scalar cutoff /0/;
scalar cutoff1 /0/;
scalar yield;

file forced /Optknock_%reg%.txt/;
put forced;
put '****************************************'//;

for (k = 1 to 1 by 1,
	counter = 1;
	n = 0;
	matrix(index,j) = 0;
	put 'Number of Interventions (k) =' k//;
	put '****************************************'//;

	while( counter = 1,
        	solve bilevel using mip maximizing z;
        	n = n + 1;
		if(n eq 1, cutoff = z.l);
        	if(bilevel.modelstat eq 1,
			if( z.l <= cutoff1 + .001,
				counter = 0;
				put "Objective value/better solution not found"/;
				cutoff1 = cutoff;
			);
                	put 'Iteration No. =' n//;
                	put 'v(%pdt%) = 'z.l:0:8/;

                	put /'knockouts :'//;
                	loop (j$(y0.l(j) = 1),
                        	put "'"j.tl:0"'   "v.l(j):0:8/;
                	);
                	put /'********************'//;

			matrix(index, j)$(ord(index) = n and (y0.l(j) = 1)) = 1;
        	);

        	if(bilevel.modelstat ne 1 or n >= 6,
                	counter = 0;
                	put "No feasible solution achieved/max n reached, modelstat : "bilevel.modelstat//;
			cutoff1 = cutoff;
        	);
	);

	put /'****************************************'/;
);

putclose forced;

