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

xMin = -5;
xMax =  5;
measurementNoiseSigma = 1;
processNoiseSigma = 1;
y    = 0.2;

x = linspace(xMin,xMax,1000);
yPrior      = normpdf(x,0,processNoiseSigma);
yLikelihood = normpdf(y,cos(2*pi*x),measurementNoiseSigma);
yPosterior  = yPrior.*yLikelihood;
yPosterior  = yPosterior./(mean(diff(x))*sum(yPosterior));
yMax        = 1.05*max(yPosterior);

leftBuffer   = 0.10;
bottomBuffer = 0.15;
topBuffer    = 0.02;
rightBuffer  = 0.02;
axesWidth    = 1-leftBuffer-rightBuffer;
axesHeight   = 1-bottomBuffer-topBuffer;

functionName = 'AnglePosteriorExample';
hFigure = figure;
FormatFigurePosition(1,6.5,2.5);
hAxes = axes('Position',[leftBuffer bottomBuffer axesWidth axesHeight]);
hLikelihood = plot(x,yLikelihood,'b','LineWidth',2);
hold on;
    hPrior     = plot(x,yPrior    ,'g','LineWidth',2);    
    hPosterior = plot(x,yPosterior,'r','LineWidth',2);
hold off;
xlabel('$x$','Interpreter','LaTeX');
ylabel('$p(y|x)$','Interpreter','LaTeX');
xlim([xMin xMax]);
ylim([0 yMax]);
box off;
FormatFigureText(16);
legend([hLikelihood,hPrior,hPosterior],'Likelihood','Prior','Posterior');
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