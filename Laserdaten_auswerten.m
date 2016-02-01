function Laserdaten_auswerten(nplot)
% This script evaluates action spectroscopy data. See explanation.txt for
% details about the procedure. Change settings in the next section:
%==========================================
%plot the fit for C60(He)n

if ~exist('nplot','var')
    nplot=20;
end

% folder to import and write data
folder='C:\Users\c7441156\Desktop\PhD\Auswertung\C60He_Spektroskopie\Bereich 2\';

% filename here
A=importdata([folder,'Bereich 2_ENERGIES.txt'],'\t',1);

% output filenames here
output_diary=[folder,'diary.txt'];
output_resonances=[folder,'resonances.txt'];
output_fit=[folder,'fit_data.txt'];
output_fitparams=[folder,'fit_parameters.txt'];
output_data=[folder,'data_normalized.txt'];

% Points used to fit the linear function that determines the shift / He
% atom
linfitpoints=2:10;%setdiff([1:30],[1,5,4,23,30]);

% fit guessing parameters
sm_width = 1; % for guessing the fit starting wavelength, smooth over sm_width nm of data
%==========================================

diary(output_diary);

% Load wavelengths
xs=A.data(:,1);
xse=A.data(:,2);

% Load Traces for each C60He_n
ys=A.data(:,3:2:end);
yse=A.data(:,4:2:end);

% Determine minima and maxima
ysmax=max(ys);
ysmin=min(ys);

% Scale data
ys=(ys-repmat(ysmin,size(ys,1),1))./(repmat(ysmax-ysmin,size(ys,1),1));
% Scale errors
yse=yse./(repmat(ysmax-ysmin,size(yse,1),1));

% Import column  names
names=A.textdata(3:2:end);

x_fit=linspace(min(xs),max(xs),100);

% define a gaussian and a lorentzian with a linear baseline
fun_l=@(x,xd) x(1)+xd*x(5)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2); %lorentz
fun_g=@(x,xd) x(1)+xd*x(5)+x(4)*normpdf(xd,x(3),abs(x(2))); %gauss

% select which one to use
fun=fun_l;

% Default values for the fit guessing
x0(1)=0.8; %constant Baseline
x0(2)=0.4; %width
x0(3)=median(xs); %center
x0(4)=-0.5; %peak height
x0(5)=0; %linear part of baseline

%fit curves
peak=0;
peakerr=0;
w=0;
wkerr=0;
rsquared=0;

% number of traces to fit
l=size(ys,2);

% prepare matrices to store fit data
paramstoplot0=zeros(l,5);
paramstoplot=zeros(l,5);
fitparams=zeros(l,11);

% Loop through all traces
for i=1:l
    % This lines smoothes over sm_width nm of the data and selects the
    % lowest value. this is used as the starting point for the
    % lorentz/gauss fit
    [~,minind]=min(smooth(xs, ys(:,i),sm_width));
    
    % select the appropriate x-value
    x0(3)=xs(minind); %reset peak center start postion to the data minimum
    
    % hard coded initial height for the lorentz/gauss
    x0(4) = -0.5;
    
    % Activate this if you want different starting parameters for a certain
    % range - used for the increase for He_n < 22 in some measurements
%     if i<21
%         x0(3) = 965.4;
%         x0(4) = 0.5;
%     end
        
    % Actually fit the data - weighted with the errors
    [x,R,~,CovB] = nlinfit(xs,ys(:,i),fun,x0,'Weights',1./yse(:,i));
    
    % calculate error with 95% probability
    % https://en.wikipedia.org/wiki/Student%27s_t-distribution#Discrete_Student.27s_t-distribution
    err=1.96*sqrt(diag(CovB));
    
    % calculate R^2 for the trace
    rsquared(i)=r2(ys(:,i),fun(x,xs),1./yse(:,i));
    
    % write fit parameters to predefined arrays
    peak(i)=x(3);
    w(i)=x(2);    
    peakerr(i)=err(3);
    werr(i)=err(2);
    
    fitparams(i,1:2:end-1)=x;
    fitparams(i,2:2:end-1)=err;
    fitparams(i,end)=rsquared(i);
    
    % print out r^2
    fprintf('%i\tR² = %f\n',i,rsquared(i));

    % save parameters for plot - x0 for the green line and x for the black
    % line
    paramstoplot0(i,:)=x0;
    paramstoplot(i,:)=x;
end

% write output data
xfun=linspace(min(xs),max(xs),100);

[p,s]=wpolyfit(linfitpoints,peak(linfitpoints),1,1./peakerr(linfitpoints));
perr = sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);

fprintf('Linear Regression:\tR² = %f\n',r2(peak(linfitpoints),polyval(p,linfitpoints),1./peakerr(linfitpoints)));
fprintf('Slope: %f+-%f nm/He\n',p(1),perr(1));

subplot(2,2,1)
plot(1:l,peak,'k.',[(1:l);(1:l)],[peak-peakerr;peak+peakerr],'k-',1:0.1:l,polyval(p,1:0.1:l),'k--');
set(gca,'ylim',[min(xs),max(xs)])

diary off

%Resonances and peak widths
fid=fopen(output_resonances,'w');
fprintf(fid,'Cluster Size\tResonance Frequency\tError\tWidth\tError\n');
fclose(fid);
dlmwrite(output_resonances,[(1:l)',peak',peakerr',w',werr'],'-append','delimiter','\t','precision','%e');
fprintf('Resonances written to %s\n',output_resonances);    

%comlete list of fitting results
fid=fopen(output_fitparams,'w');
fprintf(fid,'Cluster Size\tConst Baseline\tError\tWidth\tError\tCenter\tError\tPeak Height\tError\tlin. Baseline\tError\tr squared\n');
fclose(fid);
dlmwrite(output_fitparams,[(1:l)',fitparams],'-append','delimiter','\t','precision','%e');
fprintf('Complete list of fit parameters written to %s\n',output_fitparams);  

subplot(2,2,2)
plot(1:l,w,'k.',[(1:l);(1:l)],[w-werr;w+werr],'k-');

subplot(2,2,3:4)
%[xs,ys,xse,yse]=approx_data(xdata,ydata(:,nplot),deltax);


hold off
% plot(linspace(min(xdata),max(xdata),100),fun(paramstoplot,linspace(min(xdata),max(xdata),100)),'r-',...
%     linspace(min(xdata),max(xdata),100),fun(paramstoplot0,linspace(min(xdata),max(xdata),100)),'g--',...
%     xdata,ydata(:,nplot),'k.',xs,ys,'r.');

plot(x_fit,fun(paramstoplot(nplot,:),x_fit),'k-',...
    x_fit,fun(paramstoplot0(nplot,:),x_fit),'g--',...
    xs, ys(:,nplot),'k.');
hold on

%write fit data and normalized data
fid=fopen(output_fit,'w');
fid2=fopen(output_data,'w');
fprintf(fid,'wavelength');
fprintf(fid2,'wavelength\tError');
M=zeros(length(x_fit),length(names));
for i=1:l
    fprintf(fid,'\t%s',names{i});
    fprintf(fid2,'\t%s\tError',names{i});   
    M(:,i)=fun(paramstoplot(i,:),x_fit)';
end
fprintf(fid,'\n');
fclose(fid);
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite(output_fit,[x_fit',M],'-append','delimiter','\t','precision','%e');
fprintf('Fit data written to %s\n',output_fit);

M=zeros(length(xs),2*l);
M(:,1:2:end)=ys;
M(:,2:2:end)=yse;
dlmwrite(output_data,[xs,xse,M],'-append','delimiter','\t','precision','%e');
fprintf('Normalized exp. data written to %s\n',output_data);

plotmatx=[xs-xse,xs+xse]';
plotmaty=[ys(:,nplot),ys(:,nplot)]';
plot(plotmatx,plotmaty,'k-');

plotmatx=[xs,xs]';
plotmaty=[ys(:,nplot)-yse(:,nplot),ys(:,nplot)+yse(:,nplot)]';
plot(plotmatx,plotmaty,'k-');
hold off

%plot(xdata,ydata(:,3),'k.',xfun,fun(x,xfun),'r-');
end


function out=r2(y,ycalc,w)
    % calculates r-squared
    % https://de.wikipedia.org/wiki/Bestimmtheitsma%C3%9F#Konstruktion
    out = 1 - sum(w.*(y - ycalc).^2)/sum(w.*(y - mean(y)).^2);
end
