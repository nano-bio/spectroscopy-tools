function [xout,yout,xerr,yerr,indices] = approx_data(xdata,ydata,deltax)
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
    indices{i}=k:groups(i);
    xout(i)=mean(xdata(k:groups(i)));
    yout(i)=mean(ydata(k:groups(i)));
    xerr(i)=std(xdata(k:groups(i)));
    yerr(i)=std(ydata(k:groups(i)));
    k=groups(i)+1;
end

end