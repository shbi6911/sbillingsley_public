%preallocate results array
results = zeros(16,5);
counter = 1;    %lazy way to index, can't be bothered
for ii = 1:2        %I can use nested for loops and you can't stop me
    A = ii-1;
    for jj = 1:2
        B = jj-1;
        for kk = 1:2
            C = kk-1;
            for hh = 1:2
                D = hh-1;
                results(counter,1:4) = [A B C D];
                results(counter,5) = logicEval(A,B,C,D);
                counter = counter +1;
            end
        end
    end
end
disp(["A" "B" "C" "D" "output"]);
disp(results);

function output = logicEval(A,B,C,D)
%INPUTS     A,B,C,D       three Boolean input values
%OUTPUTS    output      output matching Boolean expression
%           !A!B!C + !BC!D + A!B!C + !AB!CD

P = ~A & ~B & ~C;
Q = ~B & C & ~D;
R = A & ~B & ~C;
S = ~A & B & ~C & D;
output = P | Q | R | S;

end