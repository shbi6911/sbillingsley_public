%provided function, not written by Shane Billingsley

% Automatic interpolation for Table A-15 (for liquid matter)
% ------------------------------------------------------------------------%
% Inputs:                                                                 %
%   variable    meaning                         unit                      %
%   ---------   ---------                       ------                    %
%   T           film temperature                C                         %
%                                                                         %
% Outputs:                                                                %
%   variable    meaning                         unit                      %
%   ---------   ---------                       ------                    %
%   k           thermal conductivity            W/m*K                     %
%   nu          kinematic viscosity             m^2/s                     %
%   Pr          Prandtl number                  none                      %
%   beta        volume expansion coefficient    1/K                       %
% ------------------------------------------------------------------------%

function [k, nu, Pr, beta] = get_A15_properties(T)
% Data from Table A-15 for liquid water [T rho k mu Pr beta]
A15_data = [0.01	999.8	0.561	0.001792	13.5   -0.000068
            5	    999.9	0.571	0.001519	11.2	0.000015
            10	    999.7	0.58	0.001307	9.45	0.0000733
            15	    999.1	0.589	0.001138	8.09	0.000138
            20	    998	    0.598	0.001002	7.01	0.000195
            25	    997	    0.607	0.000891	6.14	0.000247
            30	    996	    0.615	0.000798	5.42	0.000294
            35	    994	    0.623	0.00072	    4.83	0.000337
            40	    992.1	0.631	0.000653	4.32	0.000377
            45	    990.1	0.637	0.000596	3.91	0.000415
            50	    988.1	0.644	0.000547	3.55	0.000451
            55	    985.2	0.649	0.000504	3.25	0.000484
            60	    983.3	0.654	0.000467	2.99	0.000517
            65	    980.4	0.659	0.000433	2.75	0.000548
            70	    977.5	0.663	0.000404	2.55	0.000578
            75	    974.7	0.667	0.000378	2.38	0.000607
            80	    971.8	0.67	0.000355	2.22	0.000653
            85	    968.1	0.673	0.000333	2.08	0.00067
            90	    965.3	0.675	0.000315	1.96	0.000702
            95	    961.5	0.677	0.000297	1.85	0.000716
            100	    957.9	0.679	0.000282	1.75	0.00075
            110	    950.6	0.682	0.000255	1.58	0.000798
            120	    943.4	0.683	0.000232	1.44	0.000858
            130	    934.6	0.684	0.000213	1.33	0.000913
            140	    921.7	0.683	0.000197	1.24	0.00097
            150	    916.6	0.682	0.000183	1.16	0.001025
            160	    907.4	0.68	0.00017	    1.09	0.001145
            170	    897.7	0.677	0.00016	    1.03	0.001178
            180	    887.3	0.673	0.00015	    0.983	0.00121
            190	    876.4	0.669	0.000142	0.947	0.00128
            200	    864.3	0.663	0.000134	0.91	0.00135
            220	    840.3	0.65	0.000122	0.865	0.00152
            240	    813.7	0.632	0.000111	0.836	0.00172
            260	    783.7	0.609	0.000102	0.832	0.002
            280	    750.8	0.581	0.000094	0.854	0.00238
            300	    713.8	0.548	0.000086	0.902	0.00295];

% Separate table data into vectors for readability
T_data = A15_data(:, 1);
rho_data = A15_data(:, 2);
k_data = A15_data(:, 3);
mu_data = A15_data(:, 4);
Pr_data = A15_data(:, 5);
beta_data = A15_data(:, 6);
nu_data = mu_data./rho_data;

% Interpolation
k = interp1(T_data, k_data, T);
nu = interp1(T_data, nu_data, T);
Pr = interp1(T_data, Pr_data, T);
beta = interp1(T_data, beta_data, T);
end