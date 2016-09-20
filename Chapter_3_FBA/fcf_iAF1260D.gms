
$INLINECOM /*  */
$onlisting

OPTIONS 
       decimals = 8
       solprint = on
       reslim = 1000000
       iterlim = 10000000
       domlim = 10
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 1E-9
       work = 10000000
       mip = cplex
;

SETS
        i Set of metabolites 
$include metabolites.txt 

        j Set of reactions 
$include reactions.txt 

        reversible_f(j) Forward part of reversible rxns 
$include reversible_no_exchange_forward.txt 

        reversible_b(j) Backward part of reversible rxns 
$include reversible_no_exchange_backward.txt 

        irreversible(j)  Irreversible rxns
$include irreversible_reactions.txt 

        exchange_f(j) Exchange rxns
$include exchange_reactions_forward.txt 

        exchange_b(j) Exchange rxns
$include exchange_reactions_backward.txt 

	medium(j) Exchange rxns corresponding to compounds present in the growth medium 
$include minimal_medium.txt 

	regulation(j) Reactions that should off due to regulatory constraints 
$include regulated_reactions_glucose_aerobic_minimal.txt 

        blocked(j) The set of blocked rxns
$include blocked_reactions_glucose_aerobic_minimal.txt

        no_gpr(j)    Rxns with no GPR assocaition
$include noGPR_reactions.txt

        transport(j)  Transport reactions

        transport_like(j) Reactions whose subsystem is not Transport but transfer compounds between the cytosol and preplasm
$include transport_like_reactions.txt

        essential_inSilico(j)  In silico essential reactions
$include essential_reactions_in_silico_aerobic_minimal_glucose.txt

        essential_inVivo(j)  In vivo essential reactions
$include essential_reactions_in_vivo_aerobic_minimal_glucose.txt

        subs  Subsystem for reactions
$include subsystem_names.txt

        current_num(j) The current flux in numerator
 
        current_den(j) The current flux in denominator 

        iter /1*1000/

        fullyCoupled_counter(iter) 
;

* At first nothing is in eqncounter set
fullyCoupled_counter(iter)=no;

;

PARAMETERS
  UB(j) Lowerbound on reaction fluxes 

  LB(j) Upperbound on reaction fluxes 

  S(i,j) contains the Stoichiometric matrix of the metabolic model
$include S_matrix.txt 

  rxn_subs(j,subs)  Map of reactions and subsystems
$include reactions_subsystem.txt

  alreadyCoupled(j)

  fullyCoupled(iter,j)

  fullyCoupled_rep(iter,j) Representative of each fully coupled set

  Rmin

  Rmax

  counter1
  
  counter2

  done
;

***** Set the bounds *****
LB(j)$irreversible(j) = 0;
UB(j)$irreversible(j) = 1000;

* Reversible reactions 
LB(j)$reversible_f(j) = 0;
UB(j)$reversible_f(j) = 1000;
LB(j)$reversible_b(j) = 0;
UB(j)$reversible_b(j) = 1000;

* Exchange reactions
LB(j)$exchange_f(j) = 0;
UB(j)$exchange_f(j) = 1000;
LB(j)$exchange_b(j) = 0;
UB(j)$exchange_b(j) = 0;

***** Set the conditions for the growth medium *****
UB(j)$medium(j) = 1000;

UB('EX_glc(e)_b') = 1000;
UB('EX_o2(e)_b') = 1000;

* Turn off all reactions in the set regulation
UB(j)$(regulation(j)) = 0;

* Turn off the wild-type biomass equation
UB('Ec_biomass_iAF1260_WT_59p81M') = 0;


POSITIVE VARIABLES
        v(j)      Flux 
        r         Normalizaiton factor
;

VARIABLES
        z         Objective function 
;

v.lo(j)=LB(j);
v.up(j)=UB(j);

EQUATIONS
        obj              Objective function 
        massbalance(i)   Mass balance equations for each metabolite i
        den_const        The flux in denominator of the FCF
        glc_const        Constraints on glucose uptake
        o2_const         Constraint on oxygen uptake
        ATPM_const       Constraint on ATPM 
;

obj..             z=e= sum(j$current_num(j),v(j));
massbalance(i)..  sum(j,S(i,j)*v(j)) =e= 0;
den_const..       sum(j$current_den(j),v(j)) =e= 1;
glc_const..       v('EX_glc(e)_b') =l= 10*r;
o2_const..        v('EX_o2(e)_b') =l= 20*r;
ATPM_const..      v('ATPM') =e= 8.39*r;

************** Model definitions ********************
Model FCF 
/
  obj
  massbalance
  den_const
  glc_const
  o2_const
  ATPM_const
/;

FCF.optfile = 1;


FILE res /coupling_results.txt/;

transport(j)=no;
transport(j)$(rxn_subs(j,'Inorganic_Ion_Transport_and_Metabolism')=1)=yes;
transport(j)$(rxn_subs(j,'Transport_Inner_Membrane')=1)=yes;
transport(j)$(rxn_subs(j,'Transport_Outer_Membrane')=1)=yes;
transport(j)$(rxn_subs(j,'Transport_Outer_Membrane_Porin')=1)=yes;
transport(j)$(transport_like(j))=yes;


fullyCoupled_counter(iter)=no;
fullyCoupled(iter,j) = no;
fullyCoupled_rep(iter,j)=no;
alreadyCoupled(j)=0;

counter2=0;

alias(j,j1,j2);

LOOP(j1$([not blocked(j1)] and [alreadyCoupled(j1) = 0]),

  counter1 = 0;

  LOOP(j2$(ord(j2) > ord(j1)),
      current_num(j) = no;
      current_num(j1) = yes;

      current_den(j) = no;
      current_den(j2) = yes;

      SOLVE FCF USING LP MAXIMIZING z;
      if(FCF.modelstat = 1,
        Rmax = z.l;
      elseif (FCF.modelstat = 3),
        Rmax = -1;    
      );

      if([FCF.modelstat = 1] or [FCF.modelstat = 3],
         SOLVE FCF USING LP MINIMIZING z;
         if(FCF.modelstat = 1,
            Rmin = z.l;

            if([Rmin = 0] and [Rmax > 0],
               PUT res;
               PUT j1.tl:0," --> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
               alreadyCoupled(j2)=1;
            elseif ([Rmin > 0] and [Rmax > 0] and [Rmax > Rmin]), 
               PUT res;
               PUT j1.tl:0," <--> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
               alreadyCoupled(j2)=1;
            elseif ([Rmin > 0] and [Rmax > 0] and [Rmin = Rmax]), 
               PUT res;
               PUT j1.tl:0," <==> ",j2.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
               alreadyCoupled(j2)=1;
 
               if(counter1 = 0,
                  counter1 = counter1 + 1;
                  counter2 = counter2 + 1;
                  fullyCoupled_counter(iter)$(ord(iter)=counter2)=yes;
               );   

               fullyCoupled(iter,j1)$(ord(iter)=counter2) = 1;
               fullyCoupled(iter,j2)$(ord(iter)=counter2) = 1;

            elseif ([Rmin > 0] and [Rmax < 0]), 
               PUT res;
               PUT j2.tl:0," --> ",j1.tl:0,"  (Rmin , Rmax) = (",Rmin:0:8,",",Rmax:0:8,")"/;
            );

         else
           PUT res;
           PUT "Error! ",j1.tl:0,"  ",j2.tl:0/;
         );
     else
         PUT res;
         PUT "Error! ",j1.tl:0/;
        
     );
  );
);


************* Details of fully coupled sets *********
file res_fd /fullyCoupled_sets_details.txt/;


PUT "The total number of coupled sets is : ",card(fullyCoupled_counter):0:0//;

LOOP(fullyCoupled_counter,
    PUT res_fd;
    PUT /fullyCoupled_counter.tl:0,")"/;
    LOOP(j$(fullyCoupled(fullyCoupled_counter,j)=1),
        PUT res_fd;
        PUT j.tl:0,", ";
    );    
    PUT res_fd;
    PUT //;
);


************* Representatives of fully coupled sets *********
* This file contains the map between representatives of each fully coupled set and other 
* reactions in that set
file res_map /fullyCoupled_representatives_represented_map.txt/;

* This files contains the reactions in the fully coupled sets that are represented by another 
* reaction in that set
file res_nrp /fullyCoupled_representated.txt/;

LOOP(fullyCoupled_counter,
    done = 0;
    LOOP(j$([fullyCoupled(fullyCoupled_counter,j)=1] and [done = 0] and [not essential_inSilico(j)] and [not essential_inVivo(j)] and [not no_gpr(j)] and [not exchange_f(j)] and [not exchange_b(j)] and [not transport(j)]),
         fullyCoupled_rep(fullyCoupled_counter,j) = 1;
         done=1;
    );
);


alias(j,j1);

PUT res_map;
PUT "/"/;

PUT res_nrp;
PUT "/"/;


LOOP(fullyCoupled_counter,
  
    LOOP(j$(fullyCoupled_rep(fullyCoupled_counter,j) = 1), 
       LOOP(j1$([fullyCoupled(fullyCoupled_counter,j)=1] and [fullyCoupled_rep(fullyCoupled_counter,j) = 0]), 
           PUT res_map;
           PUT "'",j.tl:0,"'.'",j1.tl,"' 1"/;
     
           PUT res_nrp;
           PUT "'",j1.tl,"'"/;
       );
    );
    PUT /;
);

PUT res_map;
PUT "/"/;

PUT res_nrp;
PUT "/"/;
