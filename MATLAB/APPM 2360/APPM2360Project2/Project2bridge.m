%By:        Shane Billingsley
%Class:     APPM 2360 Intro to Diff Eqns with Linear Algebra
%Date:      Summer 2022


%% initialize
a=sin(pi/4);                  %matrix values
b=sin(pi/3);
c=cos(pi/3);

format compact;              %formatting
format short g;
%% solve matrix equation using mldivide
A = [-a 0 c 0 0 0 0;        %basic matrix for solving
      a 0 0 1 0 0 0;
      0 0 0 -1 1 0 0;
      0 0 -c 0 -1 0 0;
     -a -1 -b 0 0 0 0;
      a 0 0 0 0 -1 0;
      0 1 0 0 0 0 0;
      0 0 b 0 0 0 -1];
disp(rref(A))
disp(size(A))                 %check size
S =[0;0;0;0;0;0;-5000;0];    %solution matrix for solving
C = mldivide(A,S);           %solve matrix equation
disp(C);                     %display solution
%% put matrix in rref
Aa = [-a 0 c 0 0 0 0 0;        %augmented matrix for rref
      a 0 0 1 0 0 0 0;
      0 0 0 -1 1 0 0 0;
      0 0 -c 0 -1 0 0 0;
     -a -1 -b 0 0 0 0 0;
      a 0 0 0 0 -1 0 0;
      0 1 0 0 0 0 0 5000;
      0 0 b 0 0 0 -1 0];
B = rref(Aa);                  %put augmented matrix in rref
disp(Aa);                     %output matrices
disp(B);
%% solve modified matrix using mldivide and get determinant
D = [-a 0 c 0 0 0 0 0;      %new basic matrix for undetermined F3
      a 0 0 1 0 0 0 0;
      0 0 0 -1 1 0 0 0;
      0 0 -c 0 -1 0 0 0;
     -a -1 -b 0 0 0 0 0;
      a 0 0 0 0 -1 0 0;
      0 1 0 0 0 0 0 -1;
      0 0 b 0 0 0 -1 0];
G = [0;0;0;0;0;0;0;0];      %solve for undetermined F3
H = mldivide(D,G);
disp(H);
d = det(D);               %calculate determinant
disp(d);
disp(inv(D));
%% put modified matrix in rref
D = [-a 0 c 0 0 0 0 0 0;      %new augmented matrix for undetermined F3
      a 0 0 1 0 0 0 0 0;
      0 0 0 -1 1 0 0 0 0;
      0 0 -c 0 -1 0 0 0 0;
     -a -1 -b 0 0 0 0 0 0;
      a 0 0 0 0 -1 0 0 0;
      0 1 0 0 0 0 0 -1 0;
      0 0 b 0 0 0 -1 0 0];
F = rref(D);                 %put augmented matrix in rref
disp(F);