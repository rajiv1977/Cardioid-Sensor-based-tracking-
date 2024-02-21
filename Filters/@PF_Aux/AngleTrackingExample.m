%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Cardioid Sensor based tracking                                  %
%                            Auxiliary  Particle Filter                                     %
%                     Copyright @2015_DRDC, version 01_02112015                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv,  and B.Balaji                                      %
%          Defence R&D Canada, 3701 Carling Avenue, Ottawa, ON, K1A 0Z4, Canada.            %
%             rajiv.sithiravel@gmail.com and Bhashyam.Balaji@drdc-rddc.gc.ca                %
%                                                                                           %
%                                   T.Kirubarajan                                           %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                                 kiruba@mcmaster.ca                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%==============================================================================
% User-Specified Parameters
%==============================================================================
nSamples              = 200;
xMin                  = -6;
xMax                  =  6;
yMin                  = -3;
yMax                  = 3;
measurementNoiseSigma = 0.4;
processNoiseSigma     = 0.3;
nRegions              = 500;
alpha                 = 0.99;

%==============================================================================
% Axes Layout Parameters
%==============================================================================
leftBuffer   = 0.10;
bottomBuffer = 0.15;
topBuffer    = 0.02;
rightBuffer  = 0.02;
rowBuffer    = 0.05;
columnBuffer = 0.10;
axesWidth    = 1-leftBuffer-rightBuffer;
axesHeight   = 1-bottomBuffer-topBuffer;

%==============================================================================
% Create Sequence of Observations from Model
%==============================================================================
x = zeros(nSamples+1,1);
y = zeros(nSamples  ,1);

x(1) = randn*processNoiseSigma;
for n=1:nSamples
    x(n+1) = alpha*x(n)     + randn*processNoiseSigma;
    y(n  ) = cos(2*pi*x(n)) + randn*measurementNoiseSigma;
end
x = x(1:nSamples);

%==============================================================================
% Plot Observed Sequence
%==============================================================================
nRows            = 2;
nColumns         = 1;
axesWidth        = (1-leftBuffer  -rightBuffer-columnBuffer*(nColumns-1))/nColumns;
axesHeight       = (1-bottomBuffer-topBuffer  -rowBuffer   *(nRows   -1))/nRows;

functionName = 'AngleTrackingExample-ObservedSequence';
hFigure = figure;
FormatFigurePosition(1,6.5,2.5);
hAxes = axes('Position',[leftBuffer bottomBuffer axesWidth axesHeight]); % Bottom
    h = plot(0:nSamples-1,x,'r','LineWidth',2);
    xlabel('$n$','Interpreter','LaTeX');
    ylabel('$x_n$','Interpreter','LaTeX');
    xlim([0 nSamples-1]);
    ylim([xMin xMax]);
    box off;
hAxes = axes('Position',[leftBuffer bottomBuffer+axesHeight+rowBuffer axesWidth axesHeight]); % Bottom
    h = plot(0:nSamples-1,y,'b','LineWidth',2);
    ylabel('$y_n$','Interpreter','LaTeX');
    xlim([0 nSamples-1]);
    ylim([yMin yMax]);
    set(hAxes,'XTickLabel','');
    box off;    
FormatFigureText(16);
set(hFigure,'Units','Inches');
position = get(hFigure,'Position');
set(hFigure,'PaperSize',[position(3:4)]);
set(hFigure,'PaperPosition',[0 0 position(3:4)]);
saveas(hFigure,['../' functionName],'pdf');

%==============================================================================
% Recursively Estimate the Marginal Posterior
%==============================================================================
xIntegration = linspace(xMin,xMax,nRegions).';
width = mean(diff(xIntegration));

posteriorFiltered  = zeros(nRegions,nSamples);
posteriorPredicted = zeros(nRegions,1);
xHatMedian         = zeros(nSamples,1);
xLower             = zeros(nSamples,1);
xUpper             = zeros(nSamples,1);

prior                  = normpdf(xIntegration,0,processNoiseSigma);
posteriorFiltered(:,1) = normpdf(y(1),cos(2*pi*xIntegration),measurementNoiseSigma).*prior;
posteriorFiltered(:,1) = posteriorFiltered(:,1)./(width*sum(posteriorFiltered(:,1)));
iLower = 1;
iUpper = 1;
for n=2:nSamples
    for cRegion=1:nRegions
        prior                       = normpdf(xIntegration(cRegion),alpha*xIntegration,processNoiseSigma);
        posteriorPredicted(cRegion) = sum(prior.*posteriorFiltered(:,n-1))*width;
    end
    likelihood = normpdf(y(n),cos(2*pi*xIntegration),measurementNoiseSigma);
    posteriorFiltered(:,n) = likelihood.*posteriorPredicted;    
    posteriorFiltered(:,n) = posteriorFiltered(:,n)./(width*sum(posteriorFiltered(:,n)));
    iMedian = 1;
    while sum(posteriorFiltered(1:iMedian,n))<sum(posteriorFiltered(iMedian+1:end,n)), iMedian = iMedian + 1; end;
    xHatMedian(n) = xIntegration(iMedian);
    
    while width*sum(posteriorFiltered(1:iLower,n))<0.025, iLower = iLower + 1; end;
    while width*sum(posteriorFiltered(1:iLower,n))>0.025, iLower = iLower - 1; end;
    while width*sum(posteriorFiltered(1:iUpper,n))<0.975, iUpper = iUpper + 1; end;
    while width*sum(posteriorFiltered(1:iUpper,n))>0.975, iUpper = iUpper - 1; end;
    
    xLower(n) = xIntegration(iLower);
    xUpper(n) = xIntegration(iUpper);
    
end

[junk,iMax] = max(posteriorFiltered);
xHatMode = xIntegration(iMax);
xHatMean = xIntegration'*posteriorFiltered.*width;

%==============================================================================
% Plot The Filtered Posterior
%==============================================================================
nRows            = 1;
nColumns         = 1;
axesWidth        = (1-leftBuffer  -rightBuffer-columnBuffer*(nColumns-1))/nColumns;
axesHeight       = (1-bottomBuffer-topBuffer  -rowBuffer   *(nRows   -1))/nRows;

for cFigure=1:2
    hFigure = figure;
    FormatFigurePosition(cFigure,6.5,4.5);
    colormap(flipud(colormap('bone')));
    hAxes = axes('Position',[leftBuffer bottomBuffer axesWidth axesHeight]); % Bottom
        h = imagesc(0:nSamples-1,xIntegration,posteriorFiltered);
        caxis([0,prctile(posteriorFiltered(:),95)]);
        hold on;
            hTrue     = plot(0:nSamples-1,x,'r.');
            hMode     = plot(0:nSamples-1,xHatMode,'g.');
            hMean     = plot(0:nSamples-1,xHatMean,'b.');
            hMedian   = plot(0:nSamples-1,xHatMedian,'m.');
            hInterval = plot(0:nSamples-1,xUpper,'y.',0:nSamples-1,xLower,'y.'); 
        hold off;
        xlabel('$n$','Interpreter','LaTeX');
        ylabel('$x_n$','Interpreter','LaTeX');
        xlim([0 nSamples-1]);
        ylim([xMin xMax]);
        box off;
        set(hAxes,'YDir','Normal');
    FormatFigureText(16);
    legend([hTrue,hMode,hMean,hMedian,hInterval(1)],'True','Mode (MAP)','Mean','Median','95% Confidence');
    set(hFigure,'Units','Inches');
    position = get(hFigure,'Position');
    set(hFigure,'PaperSize',[position(3:4)]);
    set(hFigure,'PaperPosition',[0 0 position(3:4)]);

    switch cFigure
        case 1
            functionName = 'AngleTrackingExample-FilteredPosterior';
        case 2
            functionName = 'AngleTrackingExample-FilteredPosteriorZoom';
            xlim([nSamples-20 nSamples-1]);
            ylim([-3 3]);
    end
    saveas(hFigure,['../' functionName],'pdf');
end