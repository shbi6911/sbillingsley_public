%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

%a script to perform linear interpolation of an unknown value between two
%known values
%for a linear interpolation, you know one thing (the independent variable)
%and you want to know the corresponding dependent function value from
%tabular data.  This script asks for the known independent variable and the
%range of that variable and the associated dependent variable, and
%interpolates using the assumption that the function is linear between the
%given values

x=input("The thing you KNOW: ");
x_0=input("Range of thing you KNOW: first value: ");
x_1=input("Range of thing you KNOW: second value: ");
y_0=input("Range of thing you WANT to know: first value: ");
y_1=input("Range of thing you WANT to know: second value: ");

y=y_0+(x-x_0)*((y_1-y_0)/(x_1-x_0));

disp("The thing you want to know is:");
disp(y);