
function [U] = temperature(t,N,mat)
%temperature This function calculates temperature in a bar across position
%and time
%INPUTS     t   :a 1D vector of time values
%           N   :a scalar value for iterations of a Fourier series
%           mat :a string designating which material is being used
%               : AL=aluminum, BR=brass, ST=steel
%OUTPUTS    U   : a 2D array of temperature values ordered by x-position on
%the first direction and time value on the second

%ensure t is the correct vector type for operations
if (isrow(t))==0
    t=t';
end

const=getConst();   %define constants
%[H,T0,~]=slope_exp(const.x,const.data.(mat).ss);    %find slope of experimental data

Nvec=1:N;   %array of values of N
b_n=bncoeff(Nvec,const.(mat).H.exp,const.L);  %array of values of coefficient b_n
lambda_n=(((2.*Nvec)-1).*pi)./(2*const.L);    %array of values of coefficient lambda_n

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.exp.("alpha"+string(ii))=zeros(length(const.x),length(t));
end

for jj=1:length(const.(mat).alpha2)
    
    %create summation array
        Z_sum=zeros(length(const.x),length(t));

    for j = 1:length(const.x)     %loop through x-values (positions along the bar)
    
        X=const.x(j); %pull out the current value of x
        
        for i = 1:N     %loop through N for the summation
            
            %generate a term of the summation for the current value of n
            %(across all time values)
            Z = b_n(i) .* sin(lambda_n(i) .* X).*exp(-(lambda_n(i).^2) .* const.(mat).alpha2(jj) .* t); 
            
            %add this new term to the summation
            Z_sum(j,:) = Z_sum(j,:)+Z;
    
        end
        %calculate U(x,t) for the current x-value across time
        U.exp.("alpha"+string(jj))(j,:) = const.(mat).T0 + (const.(mat).H.exp*X) + Z_sum(j,:);
    end 
end

b_n=bncoeff(Nvec,const.(mat).H.an,const.L);  %array of values of coefficient b_n

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.an.("alpha"+string(ii))=zeros(length(const.x),length(t));
end

for jj=1:length(const.(mat).alpha2)
    
    %create summation array
        Z_sum=zeros(length(const.x),length(t));

    for j = 1:length(const.x)     %loop through x-values (positions along the bar)
    
        X=const.x(j); %pull out the current value of x
        
        for i = 1:N     %loop through N for the summation
            
            %generate a term of the summation for the current value of n
            %(across all time values)
            Z = b_n(i) .* sin(lambda_n(i) .* X).*exp(-(lambda_n(i).^2) .* const.(mat).alpha2(jj) .* t); 
            
            %add this new term to the summation
            Z_sum(j,:) = Z_sum(j,:)+Z;
    
        end
        %calculate U(x,t) for the current x-value across time
        U.an.("alpha"+string(jj))(j,:) = const.(mat).T0 + (const.(mat).H.an*X) + Z_sum(j,:);
    end 
end

end
