*** this program finds out the rxns which potentially form cycles in the network

$INLINECOM /*  */
$set myroot iAZ900/iAZ900


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


Sets

***     Metabolites
i
$include "%myroot%_cmp.txt"

***     Reactions
j
$include "%myroot%_rxnnames.txt"

*** Reactions shut off during aerobic glucose condition
offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

*** Exchange Reactions used for metabolite uptake
exchange(j)
$include "%myroot%_sourceflux_minimal.txt"

*** a dummy set
dummy(j)
;
***************************************************
Parameters

*** Reaction Type
rxntype(j)
$include "%myroot%_rxntype.txt"

*** Lower/Upper flux bounds
LB(j) , UB(j)

*** min/max flux bounds
low(j) , high(j)

*** Stoichiometry Matrix
S(i,j)

x

;

$include "%myroot%_sij.txt"

*** including rxn modifications to break cycles
*$include "%myroot%_rxn_modification.txt"

***************************************************
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
;

Obj..			z =e= sum(j$dummy(j),v(j)) ;
Stoic(i)..		sum(j, S(i,j) * v(j) )  =e=  0 ;

**************************************************
*** Get bounds for Fluxes

*scalar m /100/;
scalar M /100000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = M;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -M;
UB(j)$(rxntype(j) = 1) = M;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = M;
LB(j)$(exchange(j)) = -100;

***shutting off biomass core rxn flux

LB('biomass_core') = 0;
UB('biomass_core') = 0;

***shutting off reactions in aerobic glucose conditions

LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** putting ATP Maintenance value

LB('ATPM') = 1;
UB('ATPM') = 1;

*** Trace amount of essential nutrients that are present in experimental minimal medium (see the paper)

* 4-aminobenzoate
LB('EX_4abz(e)') = -0.5;

* biotin
LB('EX_btn(e)') = -0.5;

* inositol
LB('EX_inost(e)') = -0.5;

* nicotinate
LB('EX_nac(e)') = -0.5;

* pantothenate
LB('EX_pnto-R(e)') = -0.5;

* thiamin
LB('EX_thm(e)') = -0.5;

*** Glucose and oxygen uptake

LB('EX_glc(e)') = -100;
LB('EX_o2(e)') = -20;

v.lo(j) = LB(j);
v.up(j) = UB(j);

Model cycles /all/ ;

cycles.optfile = 1;

for(x = 1 to card(j),

	dummy(j)$(ord(j) eq x) = yes;
	dummy(j)$(ord(j) ne x) = no;

	Solve cycles using lp maximizing z;
	If(cycles.modelstat eq 1,
		high(dummy) = z.l;
	);

	Solve cycles using lp minimizing z;
        If(cycles.modelstat eq 1,
		low(dummy) = z.l;
	);
);

file ff /FindCycle.txt/;
put ff;
put "/"//;
put "*** Irreversible rxns unbounded"//;
loop(j,
	if(high(j) = M and rxntype(j) = 0 ,
		put "'"j.tl:30:0"'", low(j), "    " , high(j)/;
	);
);
put //"*** Reversible rxns unbounded"//;
loop(j,
        if((high(j) = M or low(j) = -M) and (rxntype(j) = 1 or rxntype(j) = 4),
                put "'"j.tl:30:0"'", low(j), "    " , high(j)/;
        );
);
put "/";
putclose ff;

/*
file file1 /Unbound_irr.txt/;
put file1;
put "/"//;
put "*** Irreversible Rxns unbounded"//;
loop(j,
        if(high(j) = M and rxntype(j) = 0,
                put "'"j.tl:30:0"'"/;
        );
);
put "/";
putclose file1;

file file2 /Unbound_rev_forward.txt/;
put file2;
put "/"//;
put "*** Reversible Rxns forward unbounded"//;
loop(j,
        if(high(j) = M and rxntype(j) = 1,
                put "'"j.tl:30:0"'"/;
        );
);
put "/";
putclose file2;

file file3 /Unbound_rev_backward.txt/;
put file3;
put "/"//;
put "*** Reversible Rxns backward unbounded"//;
loop(j,
        if(low(j) = -M and rxntype(j) = 1,
                put "'"j.tl:30:0"'"/;
        );
);
put "/";
putclose file3;

file file4 /Unbound_exchange.txt/;
put file4;
put "/"//;
put "*** Exchange Rxns unbounded"//;
loop(j,
        if(high(j) = M and rxntype(j) = 4,
                put "'"j.tl:30:0"'"/;
        );
);
put "/";
putclose file4;
*/
