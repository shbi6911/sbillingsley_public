%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%Homework 1 Due 1/25/24

%% Question 4.8

%standard basis vectors
i_hat = [1;0;0]; j_hat = [0;1;0]; k_hat = [0;0;1];
%given points
A = [1;2;3];    B = [4;6;5];    C=[3;9;-2];
%unit vector from A to B
i_prime = (B-A)./norm(B-A);
%z prime vector from cross product
k_prime = (cross(B-A,C-A))./norm(cross(B-A,C-A));
%get y prime by crossing z and x
j_prime = cross(k_prime,i_prime);

%construct DCM matrix from standard to prime
Q_21 = [dot(i_prime,i_hat) dot(i_prime,j_hat) dot(i_prime,k_hat);
        dot(j_prime,i_hat) dot(j_prime,j_hat) dot(j_prime,k_hat);
        dot(k_prime,i_hat) dot(k_prime,j_hat) dot(k_prime,k_hat)];
%transpose to get prime to standard
Q_12 = Q_21';
%define given vector in prime
v_prime = [2;-1;3];
%convert v from prime to standard using transpose
v = Q_12*v_prime;
%disp(v);

%% Question 4.10

%define sines and cosines
c_a = cosd(40); s_a = sind(40);
c_b = cosd(25); s_b = sind(25);

%define rotation matrices
R1 = [1 0 0; 0 c_a s_a; 0 -s_a c_a];
R2 = [c_b 0 -s_b; 0 1 0; s_b 0 c_b];

%find rotation matrix
Q = R2*R1;
%disp(Q);

%Question 3
%define sines and cosines
c = cosd(45);   s=sind(45);
%define rotation matrices
R1 = [1 0 0; 0 c s; 0 -s c];
R2 = [c 0 -s; 0 1 0; s 0 c];
R3 = [c s 0; -s c 0; 0 0 1];

%first sequence
Q1 = R3*R2*R1;  %compute DCM
disp(Q1);
[V1,D1] = eig(Q1);  %find eigenvalues and eigenvectors
u1 = V1(:,1);       %extract axis of rotation
%disp(u1);
theta1 = acosd(real(D1(2,2)));
%disp(theta1);

%second sequence
Q2 = R1*R2*R3;  %compute DCM
%disp(Q2);
[V2,D2] = eig(Q2);  %find eigenvalues and eigenvectors
%disp(D2);
u2 = V2(:,3);       %extract axis of rotation
%disp(u2);
theta2 = acosd(real(D2(1,1)));
%disp(theta3);

%third sequence
Q3 = R3*R1*R3;  %compute DCM
%disp(Q3);
[V3,D3] = eig(Q3);  %find eigenvalues and eigenvectors
%disp(D3);
u3 = V3(:,3);       %extract axis of rotation
%disp(u3);
theta3 = acosd(real(D3(1,1)));
%disp(theta3);
%% Question 11.1
