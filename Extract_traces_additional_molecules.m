%Extracts Traces with singnal normalized to the sum spectrum of IFD file 1

%IFD file calibration: LOAD H5 WITHOUT INTERNAL CALIBRATION!!!!

%==========================================================================
%========================== USER PARAMETERS ===============================
%==========================================================================
%IsotopeFit Path:
addpath('Z:\User\Josi\IsotopeFit\');

%folder
folder='Z:\Experiments\Clustof\C60 Spektroskopie Isotope Project\Final Results\';

%IFD files:
IFD_file=[];
IFD_file{1}=[folder,'H-resonance.ifd'];
%IFD_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm.ifd'];

%H5 files:
H5_file=[];
H5_file{1}=[folder,'962.2126_to_968.4814.h5'];
%H5_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm.h5'];

%Wavelength files
lambda_file=[];
lambda_file{1}=[folder,'lambda.txt'];
%lambda_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm - Wellenlaengen.txt'];

%Export filename
scan_filename=[folder,'export_traces_unscaled.txt'];

%List of molecules:
n_He=100; %number of C60-Helium traces to extract

% List of strings that should be added to C60Hen. Set to '' if nothing
% should be added
additions = {'', '[Os]', '[Ot]', '[Ou]'};

% scaling ot He series? set to 1 or 0
scaling = 0;

%-------------------- Evaluation PARAMETERS
sr=1.5; %look sr sigma above/below the current molecule
mindist=0.00001; %in nm, the minimum distance where points are combined

%==========================================================================
%======================= NOW THE MAGIC STARTS =============================
%==========================================================================

molecules=[];
i = 1;

n_additions = length(additions);

while i < n_He*n_additions
    if i==1
        for j = 1:n_additions
            molecules{i+j-1}=['[C60][He]', char(additions(j))];
        end
    else
        for j = 1:n_additions
            molecules{i+j-1}=[sprintf('[C60][He]%i',(i-1)/n_additions), char(additions(j))];
        end
    end
    i = i + n_additions;
end

% read all wavelength and sort them
% remember the files and the bin numbers!

energy_axis=[];
bufs=[];
filenum=[];
for i=1:length(lambda_file)
    temp=load(lambda_file{i});
    bufs=[bufs 1:length(temp)];
    filenum=[filenum i*ones(1,length(temp))];
    energy_axis=[energy_axis;temp];
end

[energy_axis,idx]=sort(energy_axis);
filenum=filenum(idx);
bufs=bufs(idx);

%find wavelengths which can be combined
indices={};
[eg,~,egerr,~,indices]=approx_data(energy_axis,zeros(size(energy_axis)),mindist);

%write title line to output ASCII file
fid=fopen(scan_filename,'w');
fprintf(fid,'Energy\tError');

for i=1:length(molecules)
    fprintf(fid,'\t%s\tError',molecules{i});
end

fprintf(fid,'\n');
fclose(fid);

%reference spec and IFD Data
IFD_data={};
massrange=[];


for f=1:length(IFD_file)
    data={};
    load(IFD_file{f},'-mat');
          
    data.peakdata=subtractmassoffset(data.raw_peakdata,data.calibration);
    
    IFD_data{f}=data; 
    
    %load the reference sum spectrum?
    if f==1
        %find the molecule indices
        molecule_index=find(ismember({IFD_data{1}.molecules.name},molecules));
        
        %find the part of the spectrum to look at
        minmass=min([IFD_data{1}.molecules(molecule_index).minmass])-5
        maxmass=max([IFD_data{1}.molecules(molecule_index).maxmass])+5
        
        plot_bin_length=mass2ind(data.peakdata(:,1),maxmass)-mass2ind(data.peakdata(:,1),minmass);
        ref_peakdata=data.peakdata(mass2ind(data.peakdata(:,1),minmass):mass2ind(data.peakdata(:,1),minmass)+plot_bin_length-1,:);
        
        [bgm,bgy, ~, ~]=find_bg(ref_peakdata(:,1),ref_peakdata(:,2),50,20,minmass,maxmass);
        
        %bgm=data.bgcorrectiondata.bgm;
        %bgy=data.bgcorrectiondata.bgy;
        
        %minmass=min([data.molecules(index).minmass])-sr;
        %maxmass=max([data.molecules(index).maxmass])+sr;
    end
    massrange(f,:)=mass2ind(data.peakdata(:,1),minmass):mass2ind(data.peakdata(:,1),minmass)+plot_bin_length-1;
end


%only look at the first two C60He peaks
for i=molecule_index
    IFD_data{1}.molecules(i).peakdata=IFD_data{1}.molecules(i).peakdata(1:2,:);
    IFD_data{1}.molecules(i).maxmass=IFD_data{1}.molecules(i).peakdata(2,1);
end

%we have only one write per h5 file!!
w=1;

%initialize the output matrix
output_data=zeros(length(eg),(1+length(molecule_index))*2);
output_data(:,1)=eg';
output_data(:,2)=egerr';

%extract traces for every wavelength
for i=1:length(eg); %go through all the energy groups
    fprintf('%f %%\n',100*i/length(eg));
    %sum up the signal of all bufs that belong to the current wavelength
    single_spec=zeros(1,plot_bin_length);
    for s=indices{i}
        temp=readh5buffer(H5_file{filenum(s)}, w, bufs(s));
        single_spec=single_spec+temp(massrange(filenum(s),:));
    end
              
    %design matrix for scaling: (constant baseline, linear baseline,data)
    M=[ones(plot_bin_length,1),ref_peakdata(:,1),single_spec'];
    
    %go through all the molecules
    for m=1:length(molecule_index)
        R=resolutionbycalibration(IFD_data{1}.calibration,IFD_data{1}.molecules(molecule_index(m)).com);
        
        scale_ind=findmassrange(ref_peakdata(:,1),IFD_data{1}.molecules(molecule_index(m)),R,0,sr);
        small_ind=findmassrange2(ref_peakdata(:,1),IFD_data{1}.molecules(molecule_index(m)),R,0,sr);
        
        diff_ind=setdiff(scale_ind,small_ind);
        
        %find the scaling factor
        %a=M(diff_ind,:)\ref_peakdata(diff_ind,2); %such that M*a is close to the reference spec in the msd sense
        
        % apply background correction
        %single_spec=M*a-pchip(bgm,bgy,ref_peakdata(:,1));

        if scaling==1
            output_data(i,(m+1)*2-1) = sum(single_spec(small_ind))/sum(single_spec(diff_ind));
            output_data(i,(m+1)*2)   = sqrt(1/sum(single_spec(small_ind))+1/sum(single_spec(diff_ind)))*output_data(i,(m+1)*2-1);
        else
            output_data(i,(m+1)*2-1) = sum(single_spec(small_ind));
            output_data(i,(m+1)*2)   = sqrt(sum(single_spec(small_ind)));
        end
        
        %output_data(i,(m+1)*2)   = output_data(i,(m+1)*2-1)/sqrt(sum(M(small_ind,3)));
    end
end

fprintf('dlmwrite. please wait...');
    dlmwrite(scan_filename,output_data,'-append','delimiter','\t','precision','%e');
fprintf(' done.\n');