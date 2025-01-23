%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

%preallocate results array
results = zeros(8,4);
counter = 1;    %lazy way to index, can't be bothered
for ii = 1:2        %I can use nested for loops and you can't stop me
    A = ii-1;
    for jj = 1:2
        O = jj-1;
        for kk = 1:2
            F = kk-1;
            results(counter,1:3) = [A O F];
            results(counter,4) = logicEval(A,O,F);
            counter = counter +1;
        end
    end
end
disp(["A" "O" "F" "output"]);
disp(results);


function s = logicEval(A,O,F)
%INPUTS     A,O,F       three Boolean input values
%OUTPUTS    output      output matching Boolean expression F+AO

D = A && O;
s = F | D;
end