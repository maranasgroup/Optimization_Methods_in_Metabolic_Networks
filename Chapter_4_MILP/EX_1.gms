$INLINECOM /*  */

options limrow = 1000
        optCR = 1E-9
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000
	mip = cplex;
        
***************************
Variables
z,x1

*;
*Binary Variables
y1,y2,y3;
***************************

Equations

Obj
con1
con2
;

Obj..       z =e= x1 - 4*y1 - 3*y2 - 4*y3;
con1..      x1 + y1 + y2 + y3 =g= 2;
con2..      x1 + 6*y1 + 4*y2 + 5*y3 =l= 13;
***************************
x1.lo = 0;
y1.lo = 0;
y2.lo = 0;
y3.lo = 0;
y1.up = 1;
y2.up = 1;
y3.up = 1;

y1.fx =1;
y2.fx =1;
y3.fx =0;
***************************

model example
/
Obj
con1
con2
/;

example.optfile = 1;
example.holdfixed = 1;

Solve example using mip minimizing z;

file f1 /EX_1.txt/;
put f1;
put "z =",z.l/;
put "x1 =",x1.l/;
put "y1 =",y1.l,"y2 =",y2.l,"y3 =",y3.l/;
putclose f1;