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
sigma        = 1;
nSamples     = 100;
functionName = 'GaussianPDFAsWeightedImpulses';

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
xAxisMin = -5*sigma;
xAxisMax =  5*sigma;

%xSample = 2*5*sigma*(rand(nSamples,1)-0.5);
%xSample = evrnd(0,1.5*sigma,nSamples,1);
%q       = evpdf(xSample,0,1.5*sigma);
xSample = trnd(2,nSamples,1);
q       = tpdf(xSample,2);
p       = normpdf(xSample,0,sigma);
w       = p./(q*nSamples);
w       = w/sum(w);                                        % Normalization step
x       = linspace(xAxisMin,xAxisMax,500);
y       = normpdf(x,0,sigma);

yAxisMin = 0;
yAxisMax = max(y)*1.05;

wDisplay = w*0.5*max(y)/max(w);                            % Scale the weights

%==============================================================================
% Plot the Results
%==============================================================================
hFigure = figure;
FormatFigurePosition(1,6.5,2.5);
hAxes = axes('Position',[leftBuffer bottomBuffer axesWidth axesHeight]);
plot(ones(2,1)*xSample',[zeros(1,nSamples);wDisplay'],'k','LineWidth',1.0);
hold on;
    plot(x,y,'b','LineWidth',2.0);
hold off;
xlabel('$N$','Interpreter','LaTeX');
ylabel('$E[y]$','Interpreter','LaTeX');
xlim([xAxisMin xAxisMax]);
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