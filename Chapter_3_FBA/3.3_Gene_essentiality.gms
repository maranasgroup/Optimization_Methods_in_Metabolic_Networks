******************** Example 3.3: Code for Gene Essentiality Analysis in E. coli model ********************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (Gene_essentilaity.gms) is to idenitify the set of reactions in E. coli
*       that must be active to ensure biomass production
*       In this sample code, we set the minimum threshold of biomass production at 10% of the theoretical
*       maximum biomass flux (i.e., 0.0930523 hr-1 (under aerobic) and 0.0231304 hr-1 (under anaerobic)
*       for 10 mmol/gDWhr glucose uptake)
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	Gene_essentiality.txt - GAMS output text file that   	*
*			 displays the set of in silico          *
*			 essential reactions 		        *
***************************************************************** 

***** Specifying the directories where the root files are present
$INLINECOM /*  */
$set myroot iAF1260/iAF1260

***** Specifying the minimum essentilaity threshold of biomass
$set biom 0.093052
*** NOTE: Set the theshold value at 0.0231304 under anarobic conditions

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
        

******************************** SET DEFINITIONS ***************************************
*
*       i                       = set of all metabolites
*       j                       = set of all reactions
*       offaeroglucose(j)       = set of reactions that should be turned off under aerobic glucose conditions
*       offanaero(j)            = set of reactions that should be turned off under aerobic glucose conditions
*	exchange(j)		= set of reactions defined in the minimal M9 media for uptake
*
*****************************************************************************************
Sets

i
$include "%myroot%_cmp.txt"

j
$include "%myroot%_rxnnames.txt"

offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

offanaero(j)
$include "%myroot%_offanaero.txt"

exchange(j)
$include "%myroot%_source_M9.txt"
;

****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       rxntype(j)      = Specified whether a reaction is irreversible (0), reversible
*                         (1 - forward, 2 - backward reaction), or exchange reactions (4)                  
*       high(j)         = Stores the value for the maximum biomass flux when reaction j
*                         is removed
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       c(j)            = Indicates whether a reaction is on (1) or off (0)
*       x               = counter
*
*****************************************************************************************
Parameters

S(i,j)
$include "%myroot%_sij.txt"

rxntype(j)
$include "%myroot%_rxntype.txt"

LB(j) , UB(j)

high(j)

c(j) , x
;

***************************** VARIABLE DEFINITIONS **************************************
*
*       v(j)            = Flux of a reaction j (+/-)
*       z               = Objective value (+/-)
*
*****************************************************************************************
Variables

z
v(j)
;
***************************** EQUATION DEFINITIONS **************************************
*
*       Stoic(i)        = Stoichiometric Constraint
*       Obj             = Maximize the biomass flux
*	con1a/con1b     = Allow each reaction v(j) to vary between LB(j)and UB(j) when
*			  the corressponding c(j) = 1 and fix it to zero when c(j) = 0
*
*****************************************************************************************
Equations

Obj
Stoic
con1a,con1b   
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			z =e= v('Ec_biomass_iAF1260_WT_59p81M') ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;
con1a(j)..		v(j)  =g=  LB(j)*c(j) ;
con1b(j)..		v(j)  =l=  UB(j)*c(j) ;

*****************************************************************************************

***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS *************************
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

*** Turn off biomass core reaction (i.e., setting lower and upper bounds of the reactions to zero)
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** Turn off reactions that are not active in aerobic glucose conditions
LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;
*** NOTE: Turn off offanaero(j) set of reactions instead of offaeroglucose(j)
*** for simulations under anaerobic conditions 
*LB(j)$(offanaero(j)) = 0;
*UB(j)$(offanaero(j)) = 0;

*** setting NGAM value for ATP Maintenance reaction
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** UPTAKE CONDITIONS for glucose and oxygen
LB('EX_glc(e)') = -10;
LB('EX_o2(e)') = -20;
*** NOTE: Set lower bound of oxygen uptake to zero under anaerobic condition
*LB('EX_o2(e)') = 0;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j) = LB(j);
v.up(j) = UB(j);
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
Model gene_ess
/
Obj
Stoic
con1a, con1b
/;

gene_ess.optfile = 1;
*****************************************************************************************

******SOLVING THE Linear Programming (LP) problem turning of one reaction at a time******

*** loop through reach reaction and turn it off by setting c(j) = 0 for that reaction
*** Subsequently, maximixe the biomass flux and store its value in high(j)
*** card(j) = cardinality of set j (i.e., total of reeactions in set j)
*** ord(j) = ordinality of set j (i.e., the position set j the iteration is currently in)

for(x = 1 to card(j),

	c(j)$(ord(j) ne x) = 1;
	c(j)$(ord(j) eq x) = 0;

	Solve gene_ess using lp maximizing z;
	If(gene_ess.modelstat eq 1,
		high(j)$(c(j) =0) = z.l;
	);
);
*****************************************************************************************

*******************EXPORTING THE RESULTS IN A TXT FILE***********************************

*** loop through each reaction j and check of maximum biomass achievable when the said
*** reaction is tuened off id greater that the minimum threshold

file ff /Gene_essentiality.txt/;
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
*****************************************************************************************