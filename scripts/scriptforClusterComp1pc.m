\[Phi] = 1; \[Omega] = 1;
fitUniversal[T_, Px_, Py_, \[Epsilon]_, Ex_, Ey_] := (-Ex)*Px -
   Ey*Py + 3.7309631540402794*^9*Px^2*Py^2 +
   1.7595573932168193*^9*(Px^4 + Py^4) -
   6.50780144218582*^8*(Px^4*Py^2 + Px^2*Py^4) -
      5.901102883438983*^8*(Px^6 + Py^6) +
   1.2712934148987874*^6*(Px^2 + Py^2)*(-120 +
      T) + (-2.510675608323508*^11*Px^2*Py^2 -
      4.17646670065401*^9*(Px^2 + Py^2) -

      3.685708175223359*^10*(Px^4 +
         Py^4))*\[Epsilon] + (6.486661547893217*^12*Px^2*Py^2 +
      1.2129917521706586*^11*(Px^2 + Py^2) +
      2.4051782603611133*^12*(Px^4 + Py^4))*\[Epsilon]^2;

solComp1pcAll[j_, k_] := ParallelTable[Table[Table[NSolve[{
       D[fitUniversal[T, Px, Py, -0.01., Ex, Ey], Px] == 0,
       D[fitUniversal[T, Px, Py, -0.01, Ex, Ey], Py],
       Px < 1,
       Py < 1,
       Px >= 0,
       Py >= 0},
             {Px, Py}], {T, k, j, \[Phi]}], {Ex, 0,
     150*10^5, \[Omega]*10^5}], {Ey, 0, 150*10^5, \[Omega]*10^5}];
Do[
 a = solComp1pcAll[n + 1, n];
 Export["/scratch/solComp1pc-" <> ToString[n] <> ".csv", a,
   "CSV"], {n, 1, 400, 1}
 ]
