***Finds upper and lower bounds to fluxes of Wild type based on experimental data

$INLINECOM /*  */
$set myroot iAF1260/iAF1260
$set data 2

options limrow = 1000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000;

Sets

***     Metabolites
i
$include "%myroot%_cmp.txt"

***     Reactions
j
$include "%myroot%_rxnnames.txt"

*** Exchange Reactions used for metabolite uptake
exchange(j)
$include "%myroot%_source_M9.txt"

*** fluxes shut off during Aerobic Glucose condition
offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

*** biomass precursors
biom(i)
$include "%myroot%_biom_precursors.txt"

*** Dummy set for metabolites
dummy(i)

*** Biomass precursor set
bio(j)
/
'EX_adocbl(e)'
'EX_ala-L(e)'
'EX_arg-L(e)'
'EX_asn-L(e)'
'EX_asp-L(e)'
'EX_ca2(e)'
'EX_cl(e)'
'EX_cobalt2(e)'
'EX_colipa(e)'
'EX_cu2(e)'
'EX_cys-L(e)'
'EX_enter(e)'
'EX_fe2(e)'
'EX_fe3(e)'
'EX_gln-L(e)'
'EX_glu-L(e)'
'EX_gly(e)'
'EX_gthrd(e)'
'EX_h2o(e)'
'EX_his-L(e)'
'EX_ile-L(e)'
'EX_k(e)'
'EX_leu-L(e)'
'EX_lys-L(e)'
'EX_met-D(e)'
'EX_mg2(e)'
'EX_mn2(e)'
'EX_mobd(e)'
'EX_nh4(e)'
'EX_phe-L(e)'
'EX_pheme(e)'
'EX_pro-L(e)'
'EX_ptrc(e)'
'EX_ser-L(e)'
'EX_so4(e)'
'EX_spmd(e)'
'EX_thr-L(e)'
'EX_trp-L(e)'
'EX_tyr-L(e)'
'EX_val-L(e)'
'EX_zn2(e)'
/
;
***************************************************
Parameters

*** Stoichiometry Matrix
s(i,j)
$include "%myroot%_sij.txt"

*** Dummy Stoichiometry Matrix
s_dummy(i,j)

*** Rection Type
rxntype(j)
$include "%myroot%_rxntype.txt"

*** low bound of flux
LB(j)

*** upper bound of flux
UB(j)

***max biomass flux for a precursor
high(i)

x

;

*** rxns modified to break cycles
*$include "%myroot%_rxn_modification.txt"

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

*** experimental flux and biomass constraints
$include "exp_%data%_rxnflux1.txt"

;

Obj..                           	z =e= v('Ec_biomass_iAF1260_WT_59p81M');
Stoic(i)..                 		sum(j, s_dummy(i,j) * v(j) )  =e=  0 ;

$include "exp_%data%_rxnflux2.txt"

**************************************************
*** Get bounds for Fluxes

scalar vmax /1000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = vmax;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -vmax;
UB(j)$(rxntype(j) = 1) = vmax;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = vmax;
LB(j)$(exchange(j)) = -vmax;

*** setting core biomass flux to 0
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** setting LB for biomass flux equal to max achievable s.t. exp flux constraints(max is 9.623195)
*LB('Ec_biomass_iAF1260_WT_59p81M') = 5.004;

***shutting off reactions in aerobic glucose conditions

LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** putting ATP Maintenance value

LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** Glucose and oxygen uptake

LB('Ex_glc(e)') = -100;
LB('EX_o2(e)') = -200;

v.lo(j) = LB(j);
v.up(j) = UB(j);

Model biomass
/
all
/ ;
biomass.optfile = 1;

s_dummy(i,j) = s(i,j);
alias(i,i1);
alias(j,j1);

Solve biomass using lp maximizing z ;

/*
loop(biom(i1),
	s_dummy(i,'Ec_biomass_iAF1260_WT_59p81M') = 0;
	s_dummy(i1,'Ec_biomass_iAF1260_WT_59p81M') = s(i1,'Ec_biomass_iAF1260_WT_59p81M');
	
***putting back ATP GAM reaction
*	s_dummy('atp[c]','Ec_biomass_iAF1260_WT_59p81M') = s('atp[c]','Ec_biomass_iAF1260_WT_59p81M');
*	s_dummy('adp[c]','Ec_biomass_iAF1260_WT_59p81M') = s('adp[c]','Ec_biomass_iAF1260_WT_59p81M');
*	s_dummy('h[c]','Ec_biomass_iAF1260_WT_59p81M') = s('h[c]','Ec_biomass_iAF1260_WT_59p81M');
*	s_dummy('pi[c]','Ec_biomass_iAF1260_WT_59p81M') = s('pi[c]','Ec_biomass_iAF1260_WT_59p81M');
*       s_dummy('ppi[c]','Ec_biomass_iAF1260_WT_59p81M') = s('ppi[c]','Ec_biomass_iAF1260_WT_59p81M');

     	Solve biomass using lp maximizing z ;
     	If(biomass.modelstat eq 1,
        	high(i1) = z.l;
     	);
);

file file1 /Biomass_precursor.txt/;
put file1;
put "***max biomass fluxes for each precursor"//;
loop(biom(i),
    put i.tl,"    ",high(i):0:8/;
);
putclose file1;
*/
