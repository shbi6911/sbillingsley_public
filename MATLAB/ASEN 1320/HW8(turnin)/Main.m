xData = linspace(-5,30,200); %using linspace to create an xData vector of 200 evenly spaced from -5 to 30
yData = GenerateData(xData); %we call generate data to get the yData given the xData

M = zeros(1,4); %allocating memory for the M and B vectors that will store those values for the different ranges in the piecewise
B = zeros(1,4);
xranges = [-5, 10, 15, 20, 100]; %this vector allows us to us a loop to go through each section of the peicewise that M and B must be evaluated for

for ii = 1:4 %loops 4 times to get the 4 sets of values needed
    xlog = (xData >= xranges(ii)) & (xData < xranges(ii+1)); %creates a logical vector that finds the values that are between the current range
    x1 = xData(xlog); %stores the xData and yData points where the xData is within the current range
    y1 = yData(xlog);
    [m,b] = LeastSquares(x1, y1); %call the function given the x and y coords, and returns both slope and y intercept
    M(ii) = m; %stores m and b into vectors
    B(ii) = b;
end

xFit = linspace(-5,30,5000); %using linspace to create an xFit vector of length 5000 evenly spaced from -5 to 30
yFit = PiecewiseLeastFit(M,B,xFit); %getting the yfit values using the PiecewiseLeastFit function and slope, yintercept, and xvalue inputs

makePlot(xData,yData,xFit,yFit) %finally plots the data and the best fit lines using the makePlot function