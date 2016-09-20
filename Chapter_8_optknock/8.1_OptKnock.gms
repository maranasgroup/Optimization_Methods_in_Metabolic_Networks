********************************* Example 8.1: Code for OptKnock algorithm ********************************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (OptKnock.gms) is to identify reaction removals that redirects flux
*       towards the desired chemical with simultaneous maximization of biomass flux. This code identifies
*       upto a pair of reaction removals that maximizes acetate flux in E.coli under aerobic conditions
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	OptKnock.txt -  GAMS output text file that displays  	*
*			the maximum achievable acetate flux     *
*                       with two reaction removals along with   *
*                       reaction removal strategies             *
*                       biomass                                 *
***************************************************************** 

***** Specifying the directories where the root files are present
$INLINECOM /*  */
$set myroot1 input_files/iAF1260/iAF1260
$set myroot2 input_files/
$set condition aero
$set reg with

***** Specifying the minimum threshold of biomass
$set biom 0.09305

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
*	exchange(j)		= set of reactions defined in the minimal M9 media for uptake
*       essential(j)            = set of essential in silico essential reactions under aerobic condition
*       biom(j)                 = set containing reaction for biomass flux
*       nogene(j)               = set of reactions that do not have any gpr associations
*       excg(j)                 = set of transport reactions
*       blocked(j)              = set of blocked reactions
*       index                   = set to store alternate interventions
*
*****************************************************************************************
Sets

i
$include "%myroot1%_cmp.txt"

j
$include "%myroot1%_rxnnames.txt"

exchange(j)
$include "%myroot1%_source_M9.txt"

offaeroglucose(j)
$include "%myroot1%_offaeroglu.txt"

essential(j)
$include "%myroot2%Find_essential_%condition%_%reg%.txt"

biom(j) / 'Ec_biomass_iAF1260_WT_59p81M' /

nogene(j)
$include "%myroot1%_nogene_rxns.txt"

*** Exchange rxns
excg(j)
$include "%myroot1%_ex_pp_rxns.txt"

$ONEMPTY
blocked(j) 
$include "%myroot2%Find_blocked_%condition%_%reg%.txt"

index /1*100/
;

****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       rxntype(j)      = Specified whether a reaction is irreversible (0), reversible
*                         (1 - forward, 2 - backward reaction), or exchange reactions (4)                  
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       low(j)/high(j)  = Stores the minimum/maximum flux od each reaction 
*       epsilon         = A small value
*       x               = counter
*       alt(index,j)    = Stores alternate solutions of reaction removals. For each
*                         index, a value of 1 indicates reaction is removed, and a
*                         value of 0 indicates reaction is not removed
*
*****************************************************************************************

Parameters

s(i,j)
$include "%myroot1%_sij.txt"

rxntype(j)
$include "%myroot1%_rxntype.txt"

*** flux bounds
LB(j) , UB(j)


*** matrix to store Interventions
alt(index,j)
;

***************************** VARIABLE DEFINITIONS **************************************
*       REAL VARIABLES:
*       v(j)            = Flux of a reaction j (+/-)
*       z               = Objective value (+/-)
*       lambda(i)       = Dual variable for Primal constraints on network stoichiometry
*       wL(j)/wU(j)     = Variable for lineraing product of a dual and binary variales
*       POSITIVE VARIABLES:
*       muL(j)/muU(j)   = Dual variables for Primal constraints on limits of reaction fluxes
*       BIANRY VARIABLES
*       y(j)            = Indicates 1 for reaction removal, 0 otherwise
*
******************************************************************************************
Variables
z
v(j)
lambda(i)
wL(j), wU(j)
;

Positive variables
muL(j), muU(j)
;

binary variables
y(j)
;

**************** INITIALIZING PARAMETRIC VARIABLES AND SETS ****************************

*** Scalars

scalar counter /1/;
scalar M /2000/;
scalar k ;
scalar n /0/;
scalar cutoff /0/;
scalar cutoff1 /0/;
alt(index,j) = 0;

****************************************************************************************

***************************** EQUATION DEFINITIONS **************************************
*
*       OUTER LEVEL CONSTRAINTS:
*       Obj                             = Maximize flux of target chemical
*       Primal_dual                     = Strong duality cosntraint to equate primal objective equal to dual objective
*       Int_cut                         = Integer cut constraint to find alternate solutions
*       Tot                             = Total number of interventions
*       Lin_muLa/b/c/d(j)               = Constraints to linearize muL(j)*y(j) with wL(j)
*       Lin_muUa/b/c/d(j)               = Constraints to linearize muU(j)*y(j) with wL(j)
*       PRIMAL CONSTRAINTS:
*       Stoic(i)                        = Stoichiometric Constraint
*       Flux_limL(j)/Flux_limU(j)       = Lower/Upper limits on flux v(j)
*       DUAL CONSTRAINTS:
*       dual1a                          = Dual constraint for primal biomass reaction variable
*       dual1b                          = Dual constraint for primal variables other than bioamss
*
*****************************************************************************************
Equations

Obj
Primal_dual
Int_cut
Tot
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd

Stoic
Flux_limL,Flux_limU

dual1a,dual1b
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			z =e= v('EX_ac(e)');
Tot..		        sum(j, y(j)) =e= k;
Primal_dual..		v('Ec_biomass_iAF1260_WT_59p81M') =e= sum(j, muU(j)*UB(j) - wU(j)*UB(j)) - sum(j, muL(j)*LB(j) - wL(j)*LB(j));
Int_cut(index)..	sum(j$(alt(index,j) = 0),y(j)) + sum(j$(alt(index,j) = 1),1-y(j)) =g= 1;

Lin_muLa(j)..		wL(j) =l= M*y(j);
Lin_muLb(j)..		wL(j) =g= -M*y(j);
Lin_muLc(j)..        	wL(j) =l= muL(j) + M*(1-y(j) );
Lin_muLd(j)..		wL(j) =g= muL(j) - M*(1-y(j) );

Lin_muUa(j)..           wU(j) =l= M*y(j);
Lin_muUb(j)..           wU(j) =g= -M*y(j);
Lin_muUc(j)..           wU(j) =l= muU(j) + M*(1-y(j) );
Lin_muUd(j)..           wU(j) =g= muU(j) - M*(1-y(j) );

Stoic(i)..		sum(j, (S(i,j)*v(j))) =e= 0;
Flux_limL(j)..		v(j) =g= LB(j)*(1-y(j) );
Flux_limU(j)..		v(j) =l= UB(j)*(1-y(j) );

dual1a(j)$biom(j)..		sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j)  =e= 1;
dual1b(j)$(not biom(j))..	sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j)  =e= 0;

*****************************************************************************************

***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS *************************

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

*** Turn off biomass core reaction (i.e., setting lower and upper bounds of the reactions to zero)
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*** Turn off reactions that are not active in aerobic glucose conditions
LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** setting lower bound for wild biomass flux to 10% of theoretical maximum
LB('Ec_biomass_iAF1260_WT_59p81M') = %biom%;

*** setting NGAM value for ATP Maintenance reaction
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** UPTAKE CONDITIONS for glucose and oxygen
LB('Ex_glc(e)') = -10;
LB('EX_o2(e)') = -20;

*** Fix binary variables of reactions whose removal is not allowed/have no effect on solution
y.fx(j)$essential(j) = 0;
y.fx(j)$blocked(j) = 0;
y.fx(j)$exchange(j) = 0;
y.fx(j)$nogene(j) = 0;
y.fx('ATPM') = 0;
y.fx('Ec_biomass_iAF1260_WT_59p81M') = 0;
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
model bilevel
/
Obj
Primal_dual
Int_cut
Tot
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd
Stoic
Flux_limL,Flux_limU
dual1a,dual1b
/
;

bilevel.optfile = 1;
bilevel.holdfixed = 1;
*****************************************************************************************

******SOLVING THE Linear Programming (LP) problem and EXPORTING RESULTS IN A TXT FILE****

*** 1. Iterate through total number of interventions starting from a single reaction removal
*** 2. Maximize the product flux and print its value and the interventions
*** 3. Store the intervention combination in alt(index,j)
*** 4. If an alternate solution with non-zero flux of product exists, find it;
***    else, move to higher number of intervention

file forced /Optknock.txt/;
put forced;
put '****************************************'//;

for (k = 1 to 2 by 1,
	counter = 1;
	n = 0;
	alt(index,j) = 0;
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
                	put 'v(ac) = 'z.l:0:8/;

                	put /'knockouts :'//;
                	loop (j$(y.l(j) = 1),
                        	put "'"j.tl:0"'   "v.l(j):0:8/;
                	);
                	put /'********************'//;

			alt(index, j)$(ord(index) = n and (y.l(j) = 1)) = 1;
        	);

        	if(bilevel.modelstat ne 1,
                	counter = 0;
                	put "No feasible solution achieved, modelstat : "bilevel.modelstat//;
			cutoff1 = cutoff;
        	);
	);

	put /'****************************************'/;
);

putclose forced;
*****************************************************************************************
