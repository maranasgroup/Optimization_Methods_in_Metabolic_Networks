********************************* Example 10.1: Code for k_FBA algorithm ********************************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (k_FBA.gms) is to identify the contraction in feasible flux space of a
*       genome-scale model due to introduction of additional kinetic information for a subset of reactions. 
*       This code contrasts the trade-off between maximum acetate production and biomass production in E. coli 
*       under aerobic glucose-fed conditions with a purely stoichiometry-dependent model and when supplemented
*		with kinetic information for E. coli central metabolism.
***********************************************************************************************************
*       NOTE:
*
*       Using k-FBA, the genome-scale stoichiometry matrix is divided into two parts: reactions with stoichiometric
*       information only (Jstoic), and those having additional kinetic information (Jkin). The kinetic information
*       was extracted from the kinetic model of E. coli central metabolism developed in Khodayari et al (2014, PMID:).
*       The number of reactions for kinetic representation is a compromise between reduction of solution space using
*       kinetic data and run time for solving the nonlinear expressions of mass conservations. Upon exclusion of the
*       exchange/transport reactions and elimination of reactions not involved in acetate synthesis (such as glycogen
*       pathway), a subset of the kinetic model was selected containing 36 reactions and 31 metabolites. The resulting
*       model includes reactions from glycolysis/gluconeogenesis, PP pathway, TCA cycle, anaplerotic reactions, glyoxylate
*       shunt and ED pathway with available experimental data during model parameterization. This model was finally
*       supplemented with the stoichiometric iAF1260 model of E. coli (Feist et al., 2007, PMID:).
*
*       We solve a nonlinear prgramming (NLP) problem where we maximize the target flux at different fixed biomass levels,
*       while allowing the enzymatic activity of all reactions with kinetic information to vary. Note that the 
*		procedure to incorporate the kinetic information of each reaction in Jkin described in the formulation in the book
*       [k-MPY] was altered here. Instead of directly manipulating the reaction enzyme activities (vmax), they are expressed 
*       as a function of the decomposed expressions of their elementary steps (see Khodayari et al for details). Additional constraints were imposed to express
*       the flux of each reaction in Jkin as the difference of the forward and reverse reactions for each elementary step.
*       The sum of individual enzyme fractions e is represented by m(j) (i.e., normalized total enzyme concentration) that is
*       equal to one in the reference (wild-type) strain, but varies when up/down-regulated in mutant strains. Here, we
*       allowed the m of each reaction to vary between zero (i.e., removal of its activity) and a two-fold up-regulation in
*       its expression to account for individual enzymatic perturbations in mutant strains. Likewise, the same limits of
*       variation were set for the individual enzyme fractions e for each reaction.The metabolite concentrations were allowed
*       to vary within five-fold from their steady-state values in the reference phenotype.
*		The fluxes of all other reactions limited solely by stoichiometry are allowed to vary between their designated lower
*		and upper bounds. 
*
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	k_MPY_kin.csv -GAMS output text file that displays  			*
*			the maximum achievable acetate flux at incremental 
*			biomass levels with additional kinetic constraints
*	k_MPY_stoic.csv -GAMS output text file that displays  			*
*			the maximum achievable acetate flux at incremental 
*			biomass levels with only stoichiometry-limits
***************************************************************** 
$INLINECOM /*  */
$set myroot iAF1260/iAF1260
*$set val 135.83
$set val 220
$set max_val 281
$set pdt ac
$set pdt2 ac
$set biom_stoic 9.62319
$set biom_kin 8.21621


options limrow = 1000
        optCR = 1E-9
        optCA = 1E-9
        iterlim = 100000
        decimals = 8
        reslim = 300
        work = 5000000
	nlp = conopt
	lp = cplex; 

******************************** SET DEFINITIONS ***************************************
*
*       i                       = set of all metabolites
*       j                       = set of all reactions
*       offaeroglucose(j)       = set of reactions that should be turned off under aerobic glucose conditions
*		exchange(j)				= set of reactions defined in the minimal M9 media for uptake
*       kin(j)             		= set of reactions with kinetic information
*       v_blocked(j)              = set of blocked reactions
*
*****************************************************************************************
Sets

***     Metabolites
i
$include "%myroot%_cmp.txt"

***     Reactions
j
$include "%myroot%_rxnnames.txt"

exchange(j)
$include "%myroot%_source_M9.txt"

offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

kin(j)
/
$include "kin_rxn_list.txt"
/

v_blocked(j)
$include "findblocked.txt"
;

****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       rxntype(j)      = Specified whether a reaction is irreversible (0), reversible
*                         (1 - forward, 2 - backward reaction), or exchange reactions (4)                  
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction flux
*       In addition, there are lists of additional variables:
*       Reaction_parameter_list.txt: list of kinetic parameters and their values
*****************************************************************************************
Parameters

S(i,j)
$include "%myroot%_sij.txt"

rxntype(j)
$include "%myroot%_rxntype.txt"

*** kinetic parameters
$include "Reaction_parameter_list.txt"

LB(j),UB(j)
;
Parameter c /0/;

***************************** VARIABLE DEFINITIONS **************************************
*
*       REAL VARIABLES:
*       v(j)            = Flux of a reaction j (+/-)
*       z               = Objective value (+/-)
*       m(j)            = Total Enzyme activity for reactions with kinetic expressions.
*                         Alternate representation for v_max in formulation in the book
*       In addition, there are lists of additional variables:
*       conc_list       : list of metabolite species concentration and enzyme fractions
*       elem_step_name_list     : list of elementary step fluxes for each reaction
*
******************************************************************************************
Variables
z
v(j)
m(j)
$include "conc_list.txt"
$include "elem_step_name_list.txt"
;

***************************** EQUATION DEFINITIONS **************************************
*
*       OUTER LEVEL CONSTRAINTS:
*       Obj                             = Maximize flux of target chemical
*       Stoic(i)                        = Stoichiometric Constrain
*       Con_bio_kin						= Fixes the flux biomass reaction at a specified fraction
*										  of its maximum theoretical on inclusion of kinetic information
*       Con_bio_stoic					= Fixes the flux biomass reaction at a specified fraction
*										  of its maximum theoretical in stoichiometry-only simulation
*                                         problem from u(j) value in the outer problem
*       In addition, there are lists of additional constraints
*       Conservation_of_mass_list       = Equates the sum of individual enzyme fractions e to
*                                         m
*       elem_kin_expressions            = Equates each elemetary step flux to reaction flux u(j)
*       Elementary_rxn_list             = kinetic expression of each elementary step
*
*****************************************************************************************
Equations

Obj
Stoic   stoichiometric limitation
Con_bio_kin, Con_bio_stoic

$include "Conservation_of_mass_list_1.txt"
$include "elem_kin_expressions_1.txt"
$include "Elementary_rxn_list_1.txt"
;

Obj..                          z =e= v('EX_%pdt2%(e)');
*Obj..                          z =e= v('Ec_biomass_iAF1260_WT_59p81M');
Stoic(i)..                 	sum(j, S(i,j) * v(j) )  =e=  0 ;
Con_bio_kin..               v('Ec_biomass_iAF1260_WT_59p81M') =e= %biom_kin%*0.1*c;
Con_bio_stoic..             v('Ec_biomass_iAF1260_WT_59p81M') =e= %biom_stoic%*0.1*c;

$include "Conservation_of_mass_list_2.txt"
$include "elem_kin_expressions_2.txt"
$include "Elementary_rxn_list_2.txt"

***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS, ENZYMES AND METABOLITES*****
*** Get bounds for Fluxes

scalar bound /1000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = bound;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -bound;
UB(j)$(rxntype(j) = 1) = bound;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = bound;
LB(j)$(exchange(j)) = -bound;

*** setting core biomass flux to 0
LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

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

*setting bounds for total enzyme conc.
m.lo(j)$kin(j) = 0;
m.up(j)$kin(j) = 2;

m.fx(j)$(not kin(j)) = 1;

*** bounds for metabolite/enzyme conc.
$include "conc_limits.txt"

*** bounds for elementary steps
$include "elem_step_limits.txt"
********************************

*** fixing a random value to the flux variables

v.l(j) = uniform(v.lo(j), v.up(j));
m.l(j) = uniform(m.lo(j), m.up(j));
*c.l(i) = uniform(c.lo(i), c.up(i));

*removing other glucose transporters
v.fx('GLCt2pp') = 0;
v.fx('GLCDpp') = 0;
v.fx('GLCabcpp') = 0;

v.fx(j)$v_blocked(j) = 0;
v.fx(j)$offaeroglucose(j) = 0;

************************ DECLARING THE MODELS with constraints****************************
Model flux_to_kin
/
Obj
Stoic
Con_bio_kin
$include "Conservation_of_mass_list_1.txt"
$include "elem_kin_expressions_1.txt"
$include "Elementary_rxn_list_1.txt"
/ ;

flux_to_kin.optfile = 1;
flux_to_kin.holdfixed = 1;

Model flux_to_stoic
/
Obj
Stoic
Con_bio_stoic
/ ;

flux_to_stoic.optfile = 1;
flux_to_stoic.holdfixed = 1;


******SOLVING THE NonLinear Programming (NLP) problem for trade-off plot with kinclusion of constraints and EXPORTING RESULTS IN A CSV FILE****

*** loop through incremental pertage of biomass flux and store the minimum and maximum
*** flux of acetate

file f1 /k_MPY_kin.csv/;
f1.pc=5;
put f1;

for(c = 0 to 10 by 1,
    solve flux_to_kin using nlp maximizing z;
    put z.l,;
    solve flux_to_kin using nlp minimizing z;
    put z.l:0:8,v.l('Ec_biomass_iAF1260_WT_59p81M'):0:8/;
);
*****************************************************************************************

******SOLVING THE Linear Programming (LP) problem for trade-off plot under stoichiometry-only constraints and EXPORTING RESULTS IN A CSV FILE****

*** loop through incremental pertage of biomass flux and store the minimum and maximum
*** flux of acetate

file f2 /k_MPY_stoic.csv/;
f2.pc=5;
put f2;

for(c = 0 to 10 by 1,
    solve flux_to_stoic using lp maximizing z;
    put z.l,;
    solve flux_to_stoic using lp minimizing z;
    put z.l:0:8,v.l('Ec_biomass_iAF1260_WT_59p81M'):0:8/;
);
*****************************************************************************************
