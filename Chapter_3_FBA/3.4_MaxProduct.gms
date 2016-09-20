************************* Example 3.4: Code for maximizing Succinate flux in E. coli ********************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (MaxProduct.gms) is to maximize the succinate flux in E. coli under
*       (i) no constraint on biomass and (ii) minimum 10% of theoretical maximum biomass and glucose anaerobic conditions 
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	MaxProduct.txt - GAMS output text file that displays  	*
*			 the flux of succinate under two	*
*                        constraints on minimum biomass flux
***************************************************************** 

***** Specifying the directories where the root files are present
$INLINECOM /*  */
$set myroot iAF1260/iAF1260

***** Specifying the minimum threshold of biomass
$set biom 0.0231304

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
*       c               = Indicates if constraint on minimum biomass is present (1) or not (0) 
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
*       Obj             = Maximize the exchange flux of succinate
*       Con_bio         = Constraint on minimum production of biomass
*
*****************************************************************************************
Equations

Obj
Stoic
Con_bio 
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			z =e= v('EX_succ(e)');
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;
Con_bio..               v('Ec_biomass_iAF1260_WT_59p81M') =g= %biom%*c;

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

*** Turn off reactions that are not active under anaerobic conditions
LB(j)$(offanaero(j)) = 0;
UB(j)$(offanaero(j)) = 0;

*** setting NGAM value for ATP Maintenance reaction
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** UPTAKE CONDITIONS for glucose and oxygen
LB('EX_glc(e)') = -10;
LB('EX_o2(e)') = 0;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j) = LB(j);
v.up(j) = UB(j);
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
Model maxpdt
/
Obj
Stoic
Con_bio
/;

maxpdt.optfile = 1;

******SOLVING THE Linear Programming (LP) problem and EXPORTING RESULTS IN A TXT FILE****
file ff /MaxProduct.txt/;
put ff;

c = 0;
Solve maxpdt using lp maximizing z ;
put "Maximum Succinate Flux : " z.l:0:8,"Minimum Biomass Flux : 0"//;
c = 1;
Solve maxpdt using lp maximizing z ;
put "Maximum Succinate Flux : " z.l:0:8,"Minimum Biomass Flux : %biom%"//;

putclose ff;

*****************************************************************************************