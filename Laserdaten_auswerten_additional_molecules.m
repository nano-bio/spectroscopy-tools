function Laserdaten_auswerten(nplot)
    % This script evaluates action spectroscopy data. See explanation.txt for
    % details about the procedure. Change settings in the next section:
    %==========================================
    %plot the fit for C60(He)n

    if ~exist('nplot','var')
        nplot=20;
    end

    % folder to import and write data
    folder='Z:\Experiments\Clustof\C60 Spektroskopie Isotope Project\Extraction Test with Oasch-Molecules\';

    % filename here
    A=importdata([folder,'oaschloch_export_traces3.txt'],'\t',1);

    % amount of different additions to C60He
    n_additions = 4;

    % output filenames here
    output_diary=[folder,'diary.txt'];
    output_resonances=[folder,'resonances'];
    output_fit=[folder,'fit_data'];
    output_fitparams=[folder,'fit_parameters'];
    output_data=[folder,'data_normalized'];

    % Points used to fit the linear function that determines the shift / He
    % atom
    linfitpoints=setdiff([2:32],[4,5]);

    % fit guessing parameters
    sm_width = 1; % for guessing the fit starting wavelength, smooth over sm_width nm of data
    %==========================================

    diary(output_diary);

    % Load wavelengths
    xs=A.data(:,1);
    xse=A.data(:,2);

    % Load Traces for each C60He_n + additions
    for j=1:n_additions
        ys{j}=A.data(:,3+(j-1)*2:n_additions*2:end);
        yse{j}=A.data(:,4+(j-1)*2:n_additions*2:end);

        % Determine minima and maxima
        ysmax{j}=max(ys{j});
        ysmin{j}=min(ys{j});

        % Scale data
        ys{j}=(ys{j}-repmat(ysmin{j},size(ys{j},1),1))./(repmat(ysmax{j}-ysmin{j},size(ys{j},1),1));
        % Scale errors
        yse{j}=yse{j}./(repmat(ysmax{j}-ysmin{j},size(yse{j},1),1));

        % Import column  names
        tmp = strsplit(A.textdata{1}, '\t');
        names{j}=tmp(3+(j-1)*2:n_additions*2:end);

        x_fit=linspace(min(xs),max(xs),100);

        % define a gaussian and a lorentzian with a linear baseline
        %fun_l=@(x,xd) x(1)+xd*x(5)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2); %lorentz
        fun_l=@(x,xd) x(1)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2); %lorentz
        fun_ll=@(x,xd) x(1)+x(4)/pi*(x(2))./((xd-x(3)).^2+(x(2))^2) + x(7)/pi*(x(5))./((xd-x(6)).^2+(x(5))^2); %2 lorentz
        fun_g=@(x,xd) x(1)+xd*x(5)+x(4)*normpdf(xd,x(3),abs(x(2))); %gauss

        % select which one to use
        fun=fun_l;

        % Default values for the fit guessing
        x0{j}(1)=0.8; %constant Baseline
        x0{j}(2)=0.4; %width 1
        x0{j}(3)=median(xs); %center 1
        x0{j}(4)=-0.5; %peak height 1
        %x0{j}(5)=0.4; %width 2
        %x0{j}(6)=median(xs); %center 2
        %x0{j}(7)=0.5; %peak height 2
        %x0(5)=0; %linear part of baseline

        %fit curves
        peak{j}=0;
        peakerr{j}=0;
        w{j}=0;
        wkerr{j}=0;
        rsquared{j}=0;

        % number of traces to fit
        l{j}=size(ys{j},2);

        % prepare matrices to store fit data
        paramstoplot0{j}=zeros(l{j},4);
        paramstoplot{j}=zeros(l{j},4);
        fitparams{j}=zeros(l{j},9);

        % Loop through all traces
        for i=1:l{j}
            % This lines smoothes over sm_width nm of the data and selects the
            % lowest value. this is used as the starting point for the
            % lorentz/gauss fit
            [~,minind]=min(smooth(xs, ys{j}(:,i),sm_width));

            % select the appropriate x-value
            x0{j}(3)=xs(minind); %reset peak center start postion to the data minimum

            % hard coded initial height for the lorentz/gauss
            x0{j}(4) = -0.5;

            % Activate this if you want different starting parameters for a certain
            % range - used for the increase for He_n < 22 in some measurements
        %     if i<21
        %         x0(3) = 965.4;
        %         x0(4) = 0.5;
        %     end

            % Actually fit the data - weighted with the errors
            [x{j},R{j},~,CovB{j}] = nlinfit(xs,ys{j}(:,i),fun,x0{j},'Weights',1./yse{j}(:,i));

            % calculate error with 95% probability
            % https://en.wikipedia.org/wiki/Student%27s_t-distribution#Discrete_Student.27s_t-distribution
            err{j}=1.96*sqrt(diag(CovB{j}));

            % calculate R^2 for the trace
            rsquared{j}(i)=r2(ys{j}(:,i),fun(x{j},xs),1./yse{j}(:,i));

            % write fit parameters to predefined arrays
            peak{j}(i)=x{j}(3);
            w{j}(i)=x{j}(2);    
            peakerr{j}(i)=err{j}(3);
            werr{j}(i)=err{j}(2);

            fitparams{j}(i,1:2:end-1)=x{j};
            fitparams{j}(i,2:2:end-1)=err{j};
            fitparams{j}(i,end)=rsquared{j}(i);

            % print out r^2
            fprintf('%i\tR² = %f\n',i,rsquared{j}(i));

            % save parameters for plot - x0 for the green line and x for the black
            % line
            paramstoplot0{j}(i,:)=x0{j};
            paramstoplot{j}(i,:)=x{j};
        end

        % write output data
        xfun=linspace(min(xs),max(xs),100);

        [p{j},s{j}]=wpolyfit(linfitpoints,peak{j}(linfitpoints),1,1./peakerr{j}(linfitpoints));
        perr{j} = sqrt(diag(inv(s{j}.R)*inv(s{j}.R'))./s{j}.normr.^2./s{j}.df);

        fprintf('Linear Regression:\tR² = %f\n',r2(peak{j}(linfitpoints),polyval(p{j},linfitpoints),1./peakerr{j}(linfitpoints)));
        fprintf('Slope: %f+-%f nm/He\n',p{j}(1),perr{j}(1));
        fprintf('Offset: %f+-%f nm\n',p{j}(2),perr{j}(2));

        subplot(n_additions,3,3*(j-1)+1)
        plot(1:l{j},peak{j},'k.',[(1:l{j});(1:l{j})],[peak{j}-peakerr{j};peak{j}+peakerr{j}],'k-',1:0.1:l{j},polyval(p{j},1:0.1:l{j}),'k--');
        set(gca,'ylim',[min(xs),max(xs)])
        title(['WL/nHe: ', names{j}{1}])
    end

    diary off

    for j=1:n_additions
        %Resonances and peak widths
        fn = [output_resonances, '_', num2str(j), '.txt'];
        fid=fopen(fn,'w');
        fprintf(fid,'Cluster Size\tResonance Frequency\tError\tWidth\tError\n');
        fclose(fid);
        dlmwrite(fn,[(1:l{j})',peak{j}',peakerr{j}',w{j}',werr{j}'],'-append','delimiter','\t','precision','%e');
        fprintf('Resonances written to %s\n',fn);

        %comlete list of fitting results
        fn = [output_fitparams, '_', num2str(j), '.txt'];
        fid=fopen(fn,'w');
        fprintf(fid,'Cluster Size\tConst Baseline\tError\tWidth\tError\tCenter\tError\tPeak Height\tError\tlin. Baseline\tError\tr squared\n');
        fclose(fid);
        dlmwrite(fn,[(1:l{j})',fitparams{j}],'-append','delimiter','\t','precision','%e');
        fprintf('Complete list of fit parameters written to %s\n',fn);

        subplot(n_additions,3,3*(j-1)+2)
        plot(1:l{j},w{j},'k.',[(1:l{j});(1:l{j})],[w{j}-werr{j};w{j}+werr{j}],'k-');
        title(['Peak Width/nHe: ', names{j}{1}])

        subplot(n_additions,3,3*(j-1)+3)
        %[xs,ys,xse,yse]=approx_data(xdata,ydata(:,nplot),deltax);

        hold off

        plot(x_fit,fun(paramstoplot{j}(nplot,:),x_fit),'k-',...
            x_fit,fun(paramstoplot0{j}(nplot,:),x_fit),'g--',...
            xs, ys{j}(:,nplot),'k.');
        hold on

        %write fit data and normalized data
        fn = [output_fit, '_', num2str(j), '.txt'];
        fn2 = [output_data, '_', num2str(j), '.txt'];
        fid=fopen(fn,'w');
        fid2=fopen(fn2,'w');
        fprintf(fid,'wavelength');
        fprintf(fid2,'wavelength\tError');
        M=zeros(length(x_fit),length(names{j}));
        for i=1:l{j}
            fprintf(fid,'\t%s',names{j}{i});
            fprintf(fid2,'\t%s\tError',names{j}{i});   
            M(:,i)=fun(paramstoplot{j}(i,:),x_fit)';
        end
        fprintf(fid,'\n');
        fclose(fid);
        fprintf(fid2,'\n');
        fclose(fid2);
        dlmwrite(fn,[x_fit',M],'-append','delimiter','\t','precision','%e');
        fprintf('Fit data written to %s\n',fn);

        M=zeros(length(xs),2*l{j});
        M(:,1:2:end)=ys{j};
        M(:,2:2:end)=yse{j};
        dlmwrite(fn2,[xs,xse,M],'-append','delimiter','\t','precision','%e');
        fprintf('Normalized exp. data written to %s\n',fn2);

        plotmatx{j}=[xs-xse,xs+xse]';
        plotmaty{j}=[ys{j}(:,nplot),ys{j}(:,nplot)]';
        plot(plotmatx{j},plotmaty{j},'k-');

        plotmatx{j}=[xs,xs]';
        plotmaty{j}=[ys{j}(:,nplot)-yse{j}(:,nplot),ys{j}(:,nplot)+yse{j}(:,nplot)]';
        plot(plotmatx{j},plotmaty{j},'k-');
        title(['Signal/WL for n=',num2str(nplot),' of ' names{j}{1}])
        hold off
    end

%plot(xdata,ydata(:,3),'k.',xfun,fun(x,xfun),'r-');
end


function out=r2(y,ycalc,w)
    % calculates r-squared
    % https://de.wikipedia.org/wiki/Bestimmtheitsma%C3%9F#Konstruktion
    out = 1 - sum(w.*(y - ycalc).^2)/sum(w.*(y - mean(y)).^2);
end
