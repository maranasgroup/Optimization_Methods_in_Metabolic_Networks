******************** Example 3.6: Code for Flux Variability Analysis in E. coli model ********************
*       Author: Anupam Chowdhury
***********************************************************************************************************
*       The objective of this code (Flux_variability_analysis.gms) is to generate the maximum and minimum 
*       flux of each reach at maximum level of biomass prodcution under aerobic glucose conditions.
*       Subsequently, we identify the reactions that are blocked, invariant and with variable range
***********************************************************************************************************

**** OUTPUT FILES ***********************************************
*	FVA.txt -        GAMS output text file that        	*
*			 displays blocked, invariant, and       *
*                        variable-range reactions with their    *
*                        flux ranges                            *
***************************************************************** 

***** Specifying the directories where the root files are present
$INLINECOM /*  */
$set myroot iAF1260/iAF1260

***** Specifying the theoretical maximum biomass flux
$set biom 0.93052299

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
*       dummy(j)                = a dummy set
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

dummy(j)
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
*
*****************************************************************************************
Parameters

S(i,j)
$include "%myroot%_sij.txt"

rxntype(j)
$include "%myroot%_rxntype.txt"

LB(j) , UB(j)

high(j) , low(j)

x

epsilon /1e-6/
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
*       Obj             = Maximize/Minimize flux of each reaction one at a time
*
*****************************************************************************************
Equations

Obj
Stoic
Con_bio   
;

* DEFNINING THE CONSTRAINTS *************************************************************
Obj..			z =e= sum(j$dummy(j),v(j));
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;

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

*** Set lower bound of biomass flux
LB('Ec_biomass_iAF1260_WT_59p81M') = %biom%;

*** Turn off reactions that are not active under aerobic conditions
LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** setting NGAM value for ATP Maintenance reaction
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

*** UPTAKE CONDITIONS for glucose and oxygen
LB('EX_glc(e)') = -10;
LB('EX_o2(e)') = -20;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j) = LB(j);
v.up(j) = UB(j);
*****************************************************************************************

************************ DECLARING THE MODEL with constraints****************************
Model fva
/
Obj
Stoic
/;

fva.optfile = 1;

******SOLVING THE Linear Programming (LP) problem turning of one reaction at a time*******

*** loop through each reaction; dummy(j) has one reaction in its set at each iteration
*** Next, maximize/minimize their flux and store them in high(j)/low(j)
for(x = 1 to card(j),

	dummy(j)$(ord(j) ne x) = no;
	dummy(j)$(ord(j) eq x) = yes;

	Solve fva using lp maximizing z;
	If(fva.modelstat eq 1,
		high(dummy) = z.l;
	);

	Solve fva using lp minimizing z;
        If(fva.modelstat eq 1,
		low(dummy) = z.l;
	);
);
*****************************************************************************************

*******************EXPORTING THE RESULTS IN A TXT FILE***********************************

*** loop through incremental pertage of biomass flux and store the mimum and maximum
*** flux of succinate

file ff /Flux_variability_analysis_aero.txt/;
put ff;

put "Blocked reactions"/;
loop(j$(high(j) lt epsilon and low(j) gt -epsilon),
        put j.tl:30:0, high(j):0:8,"  ", low(j):0:8/;
);

put /"Invariant reactions"/;
loop(j$((high(j) gt epsilon or low(j) lt -epsilon) and ((high(j) - low(j)) lt epsilon)),
        put j.tl:30:0, high(j):0:8,"  ", low(j):0:8/;
);

put /"Variable-range reactions"/;
loop(j$((high(j) - low(j)) gt epsilon),
        put j.tl:30:0, high(j):0:8,"  ", low(j):0:8/;
);

putclose ff
*****************************************************************************************