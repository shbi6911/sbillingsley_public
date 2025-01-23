%% Programming Homework 1 - (Shane Billingsley)

%% Task 5

% (modify this code for Task 5)
u_5 = [2.3];
[tout_5, xout_5] = ode45(@(t,x) monospinnerDynamics(t, x, u_5), [0 2], zeros(12, 1));
figure();
plotStateHistory(tout_5, xout_5);

%%
% Angular rate q, that is change in theta, initially grows the fastest.
% This is because the rotor thrust is initially oriented such that it
% causes a rotation about the j_hat axis in the body frame.  It rotates the
% monospinner about its own y-axis.  This initial rotation is quickly
% complicated by cross-coupling between angular acceleration terms.

%% Task 6

% (code for task 6 here)
u_6 = [0.0];
x_ic_6 = [0;0;0;0;0;0;1;1;-20;0;-5;0];
[tout_6, xout_6] = ode45(@(t,x) monospinnerDynamics(t, x, u_6),[0 5], x_ic_6);
figure();
plotStateHistory(tout_6, xout_6);
%%
%Our equations output inertial quantities defined in body-frame
%coordinates.  We see from the position that the monospinner follows a
%ballistic trajectory.  The oscillation of the velocity terms (particularly
%w, the z-velocity) indicate that it is rotating, and thus the body axes
%are rotating relative to the inertial frame.  This also causes the angular
%rates to change.  The angular momentum is constant in the inertial frame,
%but because it is defined in body coordinates, these change as the body
%axes rotate relative to inertial.


%% Task 7
%
%No aerodynamics forces or moments were modeled here.  The monospinner
%needs to exert enough thrust to counter gravity.  However, due to its
%construction, this creates moments about the body y and z axes.  In order
%to be stable in flight, these moments must be counteracted or the
%monospinner will rotate out of control.  Aerodynamic drag on the segments
%of the monospinner away from the CG should slow down the rotation by
%exerting a moment countering the rotation.  This is why it is necessary
%for the vehicle to be constructed in the way that it is, with large moment
%arms away from the CG.  In the video of its flight, we can see also that
%it flies with a constant rotation about its own z-axis.  This is necessary
%because of the angular momentum of the rotors.  This rotation will be kept
%stable by both drag and the moments of inertia.
