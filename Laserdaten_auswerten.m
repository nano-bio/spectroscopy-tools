function Laserdaten_auswerten(nplot,deltax)
A=importdata('Z:\Temp\rALSER\Bereich 3 fein ENERGIES_he100.txt','\t',1);

C=A.data;
C=(C-repmat(min(C),size(C,1),1))./(repmat(max(C)-min(C),size(C,1),1)); %normalize

xdata=A.data(:,1);
%ydata=C(1:length(xdata),2:(size(C,2)+1)/2);
ydata=C(1:length(xdata),2:(size(C,2)));

[xdata,idx]=sort(xdata);
ydata=ydata(idx,:);

%==========================================
%plot the fit for C60(He)n
if ~exist('nplot','var')
    nplot=1;
end
if ~exist('deltax','var')
    deltax=0.004;
end
%==========================================

output_file=sprintf('traces_smooth_%f_nm.txt',deltax);
output_resonances='resonances.txt';

%Points used to fit the linear function:
linfitpoints=setdiff([1:30],[1:6,22,23]);


fun_l=@(x,xd) x(1)+xd*x(5)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2); %lorenz
fun_g=@(x,xd) x(1)+xd*x(5)+x(4)*normpdf(xd,x(3),abs(x(2))); %gauss

fun=fun_l;

x0(1)=0.8; %constant Baseline
x0(2)=0.4; %width
x0(3)=median(xdata); %center
x0(4)=-0.5; %peak height
x0(5)=0; %linear part of baseline

%fit curves
peak=0;
peakerr=0;
w=0;
wkerr=0;
rsquared=0;

l=size(ydata,2);
data_out=[];
for i=1:l
    
    
%     [x,~,r,~,~,~,jac] = lsqcurvefit(fun,x0,xdata,ydata(:,i));

%     
%     err=2*sqrt(diag(1/(length(r)-4)*sum(r.^2)*inv(jac'*jac)));
%     
%     
    
    [xs,ys,xse,yse]=approx_data(xdata,ydata(:,i),deltax);
    
    
    
    
    [~,minind]=min(ys);
    x0(3)=xs(minind); %reset peak center start postion to the data minimum
    
    
    ix=yse>0;
    
    [x,R,~,CovB] = nlinfit(xs(ix),ys(ix),fun,x0,'Weights',1./yse(ix));
    
    if i==1
        data_out=zeros(length(xs(ix)),2*l+2);
        data_out(:,1)=xs(ix);
        data_out(:,2)=xse(ix);
    end
    
    data_out(:,2*(i+1)-1)=ys(ix);%-(x(1)+xs(ix)*x(5));
    data_out(:,2*(i+1))=yse(ix);
    
    err=2*sqrt(diag(CovB));
    
    rsquared(i)=r2(ys(ix),fun(x,xs(ix)),1./yse(ix));
    peak(i)=x(3);
    w(i)=x(2);    
    peakerr(i)=err(3);
    
    werr(i)=err(2);
    
    fprintf('%i\tR² = %f\n',i,rsquared(i));
    if i==nplot
        paramstoplot0=x0;
        paramstoplot=x;
    end
    
    
end

%write output data

%write title line
fid=fopen(output_file,'w');
fprintf(fid,'Wavelength\tWavelength Error');

for i=1:l
    fprintf(fid,'\t(C60)(He)%i\tError',i);
end
fprintf(fid,'\t\n');
fclose(fid);

fprintf('dlmwrite. please wait...');
    dlmwrite(output_file,data_out,'-append','delimiter','\t','precision','%e');
    fprintf(' Data written to %s.\n',output_file);


xfun=linspace(min(xdata),max(xdata),100);


[p,s]=wpolyfit(linfitpoints,peak(linfitpoints),1,1./peakerr(linfitpoints));
perr = sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);

fprintf('Linear Regression:\tR² = %f\n',r2(peak(linfitpoints),polyval(p,linfitpoints),1./peakerr(linfitpoints)));
fprintf('Slope: %f+-%f nm/He\n',p(1),perr(1));

subplot(2,2,1)
plot(1:l,peak,'k.',[(1:l);(1:l)],[peak-peakerr;peak+peakerr],'k-',1:0.1:l,polyval(p,1:0.1:l),'k--');
set(gca,'ylim',[min(xdata),max(xdata)])

fid=fopen(output_resonances,'w');
fprintf(fid,'Cluster Size\tResonance Frequency\tError\tWidth\tError\n');
fclose(fid);
dlmwrite(output_resonances,[(1:l)',peak',peakerr',w',werr'],'-append','delimiter','\t','precision','%e');
fprintf('Resonances written to %s\n',output_resonances);    

subplot(2,2,2)
plot(1:l,w,'k.',[(1:l);(1:l)],[w-werr;w+werr],'k-');

subplot(2,2,3:4)
[xs,ys,xse,yse]=approx_data(xdata,ydata(:,nplot),deltax);


hold off
% plot(linspace(min(xdata),max(xdata),100),fun(paramstoplot,linspace(min(xdata),max(xdata),100)),'r-',...
%     linspace(min(xdata),max(xdata),100),fun(paramstoplot0,linspace(min(xdata),max(xdata),100)),'g--',...
%     xdata,ydata(:,nplot),'k.',xs,ys,'r.');

plot(linspace(min(xdata),max(xdata),100),fun(paramstoplot,linspace(min(xdata),max(xdata),100)),'k-',...
    linspace(min(xdata),max(xdata),100),fun(paramstoplot0,linspace(min(xdata),max(xdata),100)),'g--',...
    xs,ys,'k.');
hold on

plotmatx=[xs-xse;xs+xse];
plotmaty=[ys;ys];
plot(plotmatx,plotmaty,'k-');

plotmatx=[xs;xs];
plotmaty=[ys-yse;ys+yse];
plot(plotmatx,plotmaty,'k-');
hold off

%plot(xdata,ydata(:,3),'k.',xfun,fun(x,xfun),'r-');
end

function [xout,yout,xerr,yerr] = approx_data(xdata,ydata,deltax)
%searches for points with similar x (<deltax) and combines them until all
%the points have an x distance >= deltax

[xdata, ix]=sort(xdata);
ydata=ydata(ix);

% groups=[];
% groups(1).ix=1;
% 
% k=1;
% for i=2:length(xdata)
%     if (xdata(i)-xdata(i-1))>deltax
%         %start a new group
%         k=k+1;
%         groups(k).ix=i;
%     else
%         %add the point to the current group
%         groups(k).ix=[groups(k).ix,i];
%     end
% end
% 
% l=length(groups);
% xout=zeros(1,l);
% yout=zeros(1,l);
% xerr=zeros(1,l);
% yerr=zeros(1,l);
% 
% for i=1:l
%     xout(i)=mean(xdata(groups(i).ix));
%     yout(i)=mean(ydata(groups(i).ix));
%     xerr(i)=std(xdata(groups(i).ix));
%     yerr(i)=std(ydata(groups(i).ix));
% end

dx=diff(xdata);
groups=setdiff(1:length(xdata),find(dx<deltax));

l=length(groups);
xout=zeros(1,l);
yout=zeros(1,l);
xerr=zeros(1,l);
yerr=zeros(1,l);

k=1;
for i=1:l
    xout(i)=mean(xdata(k:groups(i)));
    yout(i)=mean(ydata(k:groups(i)));
    xerr(i)=std(xdata(k:groups(i)));
    yerr(i)=std(ydata(k:groups(i)));
    k=groups(i)+1;
end

end

function out=r2(y,ycalc,w)
  out=1 - sum(w.*(y - ycalc).^2)/sum(w.*(y - mean(y)).^2);
end

    