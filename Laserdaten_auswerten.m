function Laserdaten_auswerten(nplot)
folder='F:\Laserdaten\new\';

A=importdata([folder,'Bereich 3 fein ENERGIES_test.txt'],'\t',1);

xs=A.data(:,1);
xse=A.data(:,2);

%ydata=C(1:length(xdata),2:(size(C,2)+1)/2);
ys=A.data(:,3:2:end);
yse=A.data(:,4:2:end);

ys=(ys-repmat(min(ys),size(ys,1),1))./(repmat(max(ys)-min(ys),size(ys,1),1));

names=A.textdata(3:2:end);

%==========================================
%plot the fit for C60(He)n
if ~exist('nplot','var')
    nplot=1;
end

%==========================================

output_resonances=[folder,'resonances.txt'];
output_fit=[folder,'fit_data.txt'];

%Points used to fit the linear function:
linfitpoints=setdiff([1:30],[1,5,4,23,30]);

x_fit=linspace(min(xs),max(xs),100);

fun_l=@(x,xd) x(1)+xd*x(5)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2); %lorenz
fun_g=@(x,xd) x(1)+xd*x(5)+x(4)*normpdf(xd,x(3),abs(x(2))); %gauss

fun=fun_l;

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

l=size(ys,2);

paramstoplot0=zeros(l,5);
paramstoplot=zeros(l,5);

for i=1:l
    
    [~,minind]=min(smooth(ys(:,i),10));
    x0(3)=xs(minind); %reset peak center start postion to the data minimum
        
    [x,R,~,CovB] = nlinfit(xs,ys(:,i),fun,x0,'Weights',1./yse(:,i));
        
    err=2*sqrt(diag(CovB));
    
    rsquared(i)=r2(ys(:,i),fun(x,xs),1./yse(:,i));
    peak(i)=x(3);
    w(i)=x(2);    
    peakerr(i)=err(3);
    
    werr(i)=err(2);
    
    fprintf('%i\tR² = %f\n',i,rsquared(i));

    paramstoplot0(i,:)=x0;
    paramstoplot(i,:)=x;
    
end

%write output data

xfun=linspace(min(xs),max(xs),100);

[p,s]=wpolyfit(linfitpoints,peak(linfitpoints),1,1./peakerr(linfitpoints));
perr = sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);

fprintf('Linear Regression:\tR² = %f\n',r2(peak(linfitpoints),polyval(p,linfitpoints),1./peakerr(linfitpoints)));
fprintf('Slope: %f+-%f nm/He\n',p(1),perr(1));

subplot(2,2,1)
plot(1:l,peak,'k.',[(1:l);(1:l)],[peak-peakerr;peak+peakerr],'k-',1:0.1:l,polyval(p,1:0.1:l),'k--');
set(gca,'ylim',[min(xs),max(xs)])

fid=fopen(output_resonances,'w');
fprintf(fid,'Cluster Size\tResonance Frequency\tError\tWidth\tError\n');
fclose(fid);
dlmwrite(output_resonances,[(1:l)',peak',peakerr',w',werr'],'-append','delimiter','\t','precision','%e');
fprintf('Resonances written to %s\n',output_resonances);    

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

%write fit data
fid=fopen(output_fit,'w');
fprintf(fid,'wavelength');
M=zeros(length(x_fit),length(names));
for i=1:l
    fprintf(fid,'\t%s',names{i});   
    M(:,i)=fun(paramstoplot(i,:),x_fit)';
end
fprintf(fid,'\n');
fclose(fid);
dlmwrite(output_fit,[x_fit',M],'-append','delimiter','\t','precision','%e');
fprintf('Fit data written to %s\n',output_fit);


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
  out=1 - sum(w.*(y - ycalc).^2)/sum(w.*(y - mean(y)).^2);
end

    