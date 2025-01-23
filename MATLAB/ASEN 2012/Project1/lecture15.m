%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

velocity = 0;
for i = 1:5
    velocity = velocity + (3-(0.01*velocity^2));
end
disp(velocity);