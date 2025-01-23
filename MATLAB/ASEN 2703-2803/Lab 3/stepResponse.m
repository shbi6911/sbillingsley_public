%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

num = 1;
dem = 1;
stepResp = tf(num,dem);
opts = stepDataOptions('StepAmplitude',0.5);
[x,t] = step(stepResp,opts);

plot(t,x);
xlabel("Time (s)");
ylabel("Angle (rad)");
xlim([-0.1 1]);
ylim([-0.1 1]);
hold on;