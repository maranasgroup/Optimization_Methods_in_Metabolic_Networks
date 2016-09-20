********************************* Example 11.1: Code for k_OptKnock algorithm ********************************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (k_OptKnock.gms) is to identify enzymatic changes (for reactions with
*       kinetic information) and reaction removals (for reactions with sole stoichiomatric description)
*       that achieves a biomass-coupled redirection flux towards a desired chemical consistent with core
*       kinetic description and stoichiometric constraints of the entire network. This code identifies
*       interventions that overproduces acetate in E. coli under aerobic glucose-fed conditions
***********************************************************************************************************
*       NOTE:
*
*       Using k-OptKnock, the genome-scale stoichiometry matrix is divided into two parts: reactions with stoichiometric
*       information only (Jstoic), and those having additional kinetic information (Jkin). The kinetic information
*       was extracted from the kinetic model of E. coli central metabolism developed in Khodayari et al (2014, PMID:).
*       The number of reactions for kinetic representation is a compromise between reduction of solution space using
*       kinetic data and run time for solving the nonlinear expressions of mass conservations. Upon exclusion of the
*       exchange/transport reactions and elimination of reactions not involved in succinate synthesis (such as glycogen
*       pathway), a subset of the kinetic model was selected containing 36 reactions and 31 metabolites. The resulting
*       model includes reactions from glycolysis/gluconeogenesis, PP pathway, TCA cycle, anaplerotic reactions, glyoxylate
*       shunt and ED pathway with available experimental data during model parameterization. This model was finally
*       supplemented with the stoichiometric iAF1260 model of E. coli (Feist et al., 2007, PMID:).
*
*       We solve a bilevel optimization formulation where we maximize the target flux by gradually increasing the total
*       number (k) of enzymatic interventions (for reactions in Jkin) and/or reaction deletions (for reactions in Jstoic)
*       Starting from a single intervention, we stop this procedure when the target flux does not improve appreciably with
*       additional interventions. The optimization formulations for identification of the interventions were altered from
*       the procedure described in the formulation in the book [k-OptKnock] to incorporate the kinetic information of each
*       reaction in Jkin as a function of the decomposed expressions of its elementary steps (see Khodayari et al for details)
*       instead of directly manipulating the reaction enzyme activities (vmax). Additional constraints were imposed to express
*       the flux of each reaction in Jkin as the difference of the forward and reverse reactions for each elementary step.
*       The sum of individual enzyme fractions e is represented by m(j) (i.e., normalized total enzyme concentration) that is
*       equal to one in the reference (wild-type) strain, but varies when up/down-regulated in mutant strains. Here, we
*       allowed the m of each reaction to vary between zero (i.e., removal of its activity) and a ten-fold up-regulation in
*       its expression to account for individual enzymatic perturbations in mutant strains. Likewise, the same limits of
*       variation were set for the individual enzyme fractions e for each reaction.The metabolite concentrations were allowed
*       to vary within five-fold from their steady-state values in the reference phenotype.
*
*       The outer problem worked on reactions in Jkin (using binary variables y_kin) whose enzymatic activity (i.e., m) must be
*       altered from their reference level to achieve the required flux redistribution towards succinate. The lower and upper
*       bounds on m(i.e., m_lb and m_ub) are functions of y_kin and the maximum ten-fold change. If selected for intervention,
*       m is allowed to vary from zero (m_lb=0 for y_kin=1) to a tenfold upregulation (m_ub=10 for y_kin=1) from its reference
*       expression (m=1). Otherwise, m is set to one. This information is relayed to the inner problem which maximizes biomass
*       while identifying reaction removals in Jstoic. It is noted that the flux of reaction sin Jkin are unaffected by the
*       inner problem as they are fixed by the kineticcexpressions in the outer poblem (constraint Con_kin in the code).
*
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	k_OptKnock.txt -GAMS output text file that displays  	*
*			the maximum achievable acetate flux     *
*                       with suggested interventions along with *
*                       alternate strategies of coupling acetate*
*                       production with biomass                 *
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
        optCR = 1E-4
        optCA = 1E-4
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
*       kin_rxns(j)             = set of reactions with kinetic information
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

kin_rxns(j)
/
$include "%myroot2%kin_rxn_list.txt"
/

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

*** kinetic parameters
$include "Reaction_parameter_list.txt"

*** matrix to store Interventions
alt(index,j)
;

***************************** VARIABLE DEFINITIONS **************************************
*
*       REAL VARIABLES:
*       v(j)            = Flux of a reaction j (+/-)
*       z               = Objective value (+/-)
*       lambda(i)       = Dual variable for Primal constraints on network stoichiometry
*       delta(j)        = Dual variable for Primal contraints that fix the flux for kinetic
                          part of the model
*       wL(j)/wU(j)     = Variable for lineraing product of a dual and binary variales
*       m(j)            = Total Enzyme activity for reactions with kinetic expressions.
*                         Alternate representation for v_max in formulation in the book
*       u(j)            = Flux of reaction in Jkin
*       In addition, there are lists of additional variables:
*       conc_list       : list of metabolite species concentration and enzyme fractions
*       elem_step_name_list     : list of elementary step fluxes for each reaction
*
*       POSITIVE VARIABLES:
*       muL(j)/muU(j)   = Dual variables for Primal constraints on limits of reaction fluxes
*
*       BIANRY VARIABLES
*       y_kin(j)        = Identifies interventions in the kinetic part of the model.
*                         Indicates 1 for enzyme that is perturbed from its reference activity,
*                         0 otherwise
*       y_stoic(j)      = Identifies interventions in the stoichiometric part of the model.
*                         Indicates 1 for reaction removal, 0 otherwise
*
******************************************************************************************
Variables
z
v(j)
lambda(i)
wL(j), wU(j)
m(j)

*** metabolite/complex concentrations
$include "%myroot2%conc_list.txt"

***elementary_step_names
$include "%myroot2%elem_step_name_list.txt"
;

Positive variables
muL(j), muU(j)
;

binary variables
y_kin(j), y_stoic(j)
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
*       Primal_dual                     = Strong duality constraint to equate primal objective
*                                         equal to dual objective
*       Int_cut                         = Integer cut constraint to find alternate solutions
*       Tot                             = Total number of interventions
*       VmaxLimL(j)/VmaxLimL(j)         = Lower/Upper limits on enzyme activity m(j)
*       Lin_muLa/b/c/d(j)               = Constraints to linearize muL(j)*y_stoic(j) with wL(j)
*       Lin_muUa/b/c/d(j)               = Constraints to linearize muU(j)*y_stoic(j) with wL(j)
*       In addition, there are lists of additional constraints
*       Conservation_of_mass_list       = Equates the sum of individual enzyme fractions e to
*                                         m
*       elem_kin_expressions            = Equates each elemetary step flux to reaction flux u(j)
*       Elementary_rxn_list             = kinetic expression of each elementary step
*
*       PRIMAL CONSTRAINTS:
*       Stoic(i)                        = Stoichiometric Constraint
*       Flux_limL(j)/Flux_limU(j)       = Lower/Upper limits on flux v(j)
*       Con_kin(j)                      = Fixes the flux of reaction in Jkin in the inner
*                                         problem from u(j)value in the outer problem
*       DUAL CONSTRAINTS:
*       dual1a                          = Dual constraint for primal biomass reaction variable
*       dual1b                          = Dual constraint for primal variables other than bioamss
*       dual1c                          = Dual constraint for primal variables in Jkin
*
*****************************************************************************************
Equations

Obj
Primal_dual
Int_cut
Tot
Vmax_LimL
Vmax_LimU
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd

$include "%myroot2%Conservation_of_mass_list_1.txt"
$include "%myroot2%elem_kin_expressions_1.txt"
$include "%myroot2%Elementary_rxn_list_1.txt"

Stoic
Flux_limL,Flux_limU

dual1a,dual1b,dual1c
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			        z =e= v('EX_ac(e)');
Tot..		                sum(j, y_kin(j) + y_stoic(j)) =e= k;
Primal_dual..		        v('Ec_biomass_iAF1260_WT_59p81M') =e= sum(j, muU(j)*UB(j) - wU(j)*UB(j)) - sum(j, muL(j)*LB(j) - wL(j)*LB(j)) +sum(j$kin_rxns(j), delta(j)*u(j)) ;
Int_cut(index)..	        sum(j$(alt(index,j) = 0),y_kin(j) + y_stoic(j)) + sum(j$(alt(index,j) = 1),1-y_kin(j)-y_stoic(j)) =g= 1;

Vmax_LimL(j)$kin_rxns(j)..	m(j) =g= 1-y_kin(j);
Vmax_LimU(j)$kin_rxns(j)..	m(j) =l= 1+9*y_kin(j) ;

$include "%myroot2%Conservation_of_mass_list_2.txt"
$include "%myroot2%elem_kin_expressions_2.txt"
$include "%myroot2%Elementary_rxn_list_2.txt"

Lin_muLa(j)..		        wL(j) =l= M*y_stoic(j);
Lin_muLb(j)..		        wL(j) =g= -M*y_stoic(j);
Lin_muLc(j)..        	        wL(j) =l= muL(j) + M*(1-y_stoic(j) );
Lin_muLd(j)..		        wL(j) =g= muL(j) - M*(1-y_stoic(j) );

Lin_muUa(j)..                   wU(j) =l= M*y_stoic(j);
Lin_muUb(j)..                   wU(j) =g= -M*y_stoic(j);
Lin_muUc(j)..                   wU(j) =l= muU(j) + M*(1-y_stoic(j) );
Lin_muUd(j)..                   wU(j) =g= muU(j) - M*(1-y_stoic(j) );

Stoic(i)..		        sum(j, (S(i,j)*v(j))) =e= 0;
Flux_limL(j)..		        v(j) =g= LB(j)*(1-y_stoic(j) );
Flux_limU(j)..		        v(j) =l= UB(j)*(1-y_stoic(j) );
Con_kin(j)$(kin_rxns(j))..      v(j) =e= u(j);

dual1a(j)$biom(j)..		sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j)  =e= 1;
dual1b(j)$(not biom(j) and not kin_rxns(j))..	sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j)  =e= 0;
dual1c(j)$kin_rxns(j)..	        sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j) + delta(j)  =e= 0;

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
y_kin.fx(j)$(not kin_rxns(j)) = 0;
y_stoic.fx(j)$kin_rxns(j) = 0;
y_stoic.fx(j)$essential(j) = 0;
y_stoic.fx(j)$blocked(j) = 0;
y_stoic.fx(j)$exchange(j) = 0;
y_stoic.fx(j)$nogene(j) = 0;
y_stoic.fx('ATPM') = 0;
y_stoic.fx('Ec_biomass_iAF1260_WT_59p81M') = 0;

********************************

*** setting limits for fold change of total enzyme activity
m.lo(j)$kin_rxns(j) = 0;
m.up(j)$kin_rxns(j) = 10;

*** fixing enzyme conc. for rxns without kin expressions
m.fx(j)$(not kin_rxns(j)) = 1;
********************************

*** bounds for metabolite/enzyme conc.
$include "%myroot2%conc_limits.txt"

*** bounds for elementary steps
$include "%myroot2%elem_step_limits.txt"
********************************

*** fixing a random value to the flux and enzyme activity

v.l(j) = uniform(v.lo(j), v.up(j));
m.l(j) = uniform(m.lo(j), m.up(j));
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
model bilevel
/
Obj
Primal_dual
Int_cut
Tot
Vmax_LimL
Vmax_LimU
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd
Stoic
Flux_limL,Flux_limU
Con_kin
dual1a,dual1b,dual1c

$include "%myroot2%Conservation_of_mass_list_1.txt"
$include "%myroot2%elem_kin_expressions_1.txt"
$include "%myroot2%Elementary_rxn_list_1.txt"
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

for (k = 1 to 8 by 1,
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
                	loop (j$(y_kin.l(j) = 1 or y_stoic.l(j) = 1),
                        	put "'"j.tl:0"'   "v.l(j):0:8/;
                	);
                	put /'********************'//;

			alt(index, j)$(ord(index) = n and (y_kin.l(j) = 1 or y_stoic.l(j) = 1)) = 1;
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
