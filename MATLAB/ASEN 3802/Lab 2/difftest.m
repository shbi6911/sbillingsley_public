%This script REQUIRES that Lab2main be run first to generate temperature
%models
%Difftest generates 11 plots of different thermal diffusivity values
%against experimental data, to find which one has the best agreement with
%experiment

%This variable designates the material under investigation for this run.
mat="Brass26"; %Choose from {"steel","Alum25","Alum28","Brass26","Brass29"}

%loop through alpha values
for kk=1:11
    figure(kk);
    hold on;

    %loop to plot all thermocouples
    for ll=1:8
            plot(const.data.(mat).time,const.data.(mat).("th"+string(ll)),'r');
            plot(const.data.(mat).time,U2.(mat).exp.model2.("alpha"+string(kk))(ll,:),'k');
    end
    title(mat+" Alpha"+string(kk));
    xlabel ("Time (s)");
    ylabel ("Temperature (degrees C)");
    hold off;
end