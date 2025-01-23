%%printing plots
fig_start = figs.T2_Q1a(1);
fig_end = figs.T2_Q5f(end);
counter = 0;
for ii = fig_start:fig_end
    names = string(fieldnames(figs));
    probnum = floor(counter/6)+1;
    fignum = mod(counter,6)+1;
    filename = names(probnum) + "_" + string(fignum);
    disp(filename);
    saveas(figure(i),filename);
    counter = counter+1;
end