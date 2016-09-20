******************** Example 3.5: Code for Flux Tradeoff Analysis in E. coli model ********************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (Flux_tradeoff.gms) is to generate the tradeoff plot between ethanol
*       and biomass (i.e., maximum and minimum ethanol production at incremental levels of biomass 
*       production) under anaerobic glucose conditions
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	Flux_tradeoff.txt - GAMS output text file that   	*
*			 displays the range of succiante flux   *
*			 at incremental levels of biomass       *
***************************************************************** 

***** Specifying the directories where the root files are present
$INLINECOM /*  */
$set myroot iAF1260/iAF1260

***** Specifying the theoretical maximum biomass flux
$set biom 0.16204927

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
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       c               = Indicates the percentage of biomass flux with respect to its
*                         theoretical maximum (0% to 100%)
*
*****************************************************************************************
Parameters

S(i,j)
$include "%myroot%_sij.txt"

rxntype(j)
$include "%myroot%_rxntype.txt"

LB(j) , UB(j)

c
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
*	Con_bio         = Setting the biomass flux to a value ranging uniformly from 0%
*                         to 100% of its theoretical maximum
*
*****************************************************************************************
Equations

Obj
Stoic
Con_bio   
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			z =e= v('EX_succ(e)') ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;
Con_bio..               v('Ec_biomass_iAF1260_WT_59p81M') =l= %biom%*0.1*c;

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
*LB(j)$(offaeroglucose(j)) = 0;
*UB(j)$(offaeroglucose(j)) = 0;
*** NOTE: Turn off offanaero(j) set of reactions instead of offaeroglucose(j)
*** for simulations under anaerobic conditions 
LB(j)$(offanaero(j)) = 0;
UB(j)$(offanaero(j)) = 0;

*** setting NGAM value for ATP Maintenance reaction
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** UPTAKE CONDITIONS for glucose and oxygen
LB('EX_glc(e)') = -10;
*LB('EX_o2(e)') = -20;
*** NOTE: Set lower bound of oxygen uptake to zero under anaerobic condition
LB('EX_o2(e)') = 0;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j) = LB(j);
v.up(j) = UB(j);

v.fx('LDH_D_f') = 0;
v.fx('LDH_D2') = 0;
v.fx('L-LACD2') = 0;
v.fx('L-LACD3') = 0;
v.fx('PFL') = 0;
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
Model flux_to
/
Obj
Stoic
Con_bio
/;

flux_to.optfile = 1;

******SOLVING THE Linear Programming (LP) problem and EXPORTING RESULTS IN A TXT FILE****

*** loop through incremental pertage of biomass flux and store the minimum and maximum
*** flux of ethanol

file ff /Flux_tradeoff_anaero.txt/;
ff.pc=5;
put ff;

put "max_%pdt%","min_%pdt%","biomass"/;
for(c = 0 to 10 by 1,
    solve flux_to using lp maximizing z;
    put z.l,;
    solve flux_to using lp minimizing z;
    put z.l:0:8,v.l('Ec_biomass_iAF1260_WT_59p81M'):0:8/;
);
*****************************************************************************************