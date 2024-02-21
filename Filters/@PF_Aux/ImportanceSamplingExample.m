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
sigma       = 1;
nSamplesMax = 2000;

%==============================================================================
% Axes Layout Parameters
%==============================================================================
leftBuffer   = 0.10;
bottomBuffer = 0.15;
topBuffer    = 0.02;
rightBuffer  = 0.02;
axesWidth    = 1-leftBuffer-rightBuffer;
axesHeight   = 1-bottomBuffer-topBuffer;

%==============================================================================
% Processing
%==============================================================================
x             = linspace(-5*sigma,5*sigma,1000);
delta         = mean(diff(x));
meanNumerical = delta*sum(cos(x).*normpdf(x,0,sigma));

xSampled      = 5*2*sigma*(rand(nSamplesMax,1)-0.5);
importancePDF = 1/(5*2*sigma);
weights       = normpdf(xSampled,0,sigma)./importancePDF;
meanAveraged  = cumsum(cos(xSampled).*weights)./(1:nSamplesMax).';
    
yAxisMin = 0;
yAxisMax = 1.5;

%==============================================================================
% Plot the Results
%==============================================================================
functionName = 'ImportanceSamplingExample';
hFigure = figure;
FormatFigurePosition(1,6.5,2.5);
hAxes = axes('Position',[leftBuffer bottomBuffer axesWidth axesHeight]);
h = plot(1:nSamplesMax,meanAveraged,'b','LineWidth',2.0);
hold on;
    h = plot([0 nSamplesMax],meanNumerical*[1 1],'r','LineWidth',2.0);
hold off;
xlabel('$N$','Interpreter','LaTeX');
ylabel('$E[y]$','Interpreter','LaTeX');
xlim([1 nSamplesMax]);
ylim([yAxisMin yAxisMax]);
box off;
FormatFigureText(16);
set(hFigure,'Units','Inches');
position = get(hFigure,'Position');
set(hFigure,'PaperSize',[position(3:4)]);
set(hFigure,'PaperPosition',[0 0 position(3:4)]);
saveas(hFigure,['../' functionName],'pdf');

% fileIdentifier = fopen(['../' functionName '.tex'],'w');
% fprintf(fileIdentifier,'\\begin{frame}\n');
% fprintf(fileIdentifier,'\t\\frametitle{Matlab Code}\n');
% fprintf(fileIdentifier,'\t\\matlabcode{Matlab/%s.m}\n',functionName);    
% fprintf(fileIdentifier,'\\end{frame}');
% fprintf(fileIdentifier,'\n');
% fclose(fileIdentifier);