%Extracts Traces with singnal normalized to the sum spectrum of IFD file 1

clear all;

%IFD file calibration: LOAD H5 WITHOUT INTERNAL CALIBRATION!!!!

%==========================================================================
%========================== USER PARAMETERS ===============================
%==========================================================================
%IsotopeFit Path:
addpath('Z:\Experiments\STM\matlab\IsotopeFit\');

%folder
folder='Z:\Experiments\Clustof\C60He + Spektroskopie\ID 3\';

%IFD files:
IFD_file=[];
IFD_file{1}=[folder,'964.7208_to_966.0902_gaussian.ifd'];
%IFD_file{2}=[folder,'Bereich 3 fein 1 960nm.ifd'];

%H5 files:
H5_file=[];
H5_file{1}=[folder,'964.7208_to_966.0902.h5'];
%H5_file{2}=[folder,'Bereich 3 fein 1 960nm.h5'];

%Wavelength files
lambda_file=[];
lambda_file{1}=[folder,'lambda.txt'];
%lambda_file{2}=[folder,'Bereich 3 fein 1 960nm - Wellenlaengen.txt'];

%Export filename
scan_filename_left=[folder,'traces_left.txt'];
scan_filename_right=[folder,'traces_right.txt'];

%List of molecules:
n_He=1:50; %number of C60-Helium traces to extract

%-------------------- Evaluation PARAMETERS
sr=1.3; %look sr sigma above/below the current molecule
mindist=0.01; %in nm, the minimum distance where points are combined


%==========================================================================
%======================= NOW THE MAGIC STARTS =============================
%==========================================================================

molecules=[];
for i=n_He
    if i==1
        molecules{i}='[C60][He]';
    else
        molecules{i}=sprintf('[C60][He]%i',i);
    end
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
fid=fopen(scan_filename_left,'w');
fprintf(fid,'Energy\tError');

for i=1:length(molecules)
    fprintf(fid,'\t%s\tError',molecules{i});
end

fprintf(fid,'\n');
fclose(fid);

fid=fopen(scan_filename_right,'w');
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
    [~,idx]=sort(data.peakdata(:,1));
    data.peakdata=data.peakdata(idx,:);
    
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
C60_ind=[];
for i=molecule_index
    IFD_data{1}.molecules(i).peakdata=IFD_data{1}.molecules(i).peakdata(1:2,:);
    IFD_data{1}.molecules(i).maxmass=IFD_data{1}.molecules(i).peakdata(2,1);
    
    R=resolutionbycalibration(IFD_data{1}.calibration,IFD_data{1}.molecules(i).com);
    
    C60_ind=[C60_ind, findmassrange2(ref_peakdata(:,1)',IFD_data{1}.molecules(i),R,0,sr)];
       
end
% 
% diff_ind=setdiff(500:size(ref_peakdata,1),C60_ind); %start at some offset ot avoid the c60 peak
    
%we have only one write per h5 file!!
w=1;

%initialize the output matrix
output_data_left=zeros(length(eg),(1+length(molecule_index))*2);
output_data_left(:,1)=eg';
output_data_left(:,2)=egerr';

output_data_right=output_data_left;


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
        
        %scale_ind=findmassrange(ref_peakdata(:,1),IFD_data{1}.molecules(molecule_index(m)),R,0,sr); %[+ + +]
        small_ind=findmassrange2(ref_peakdata(:,1),IFD_data{1}.molecules(molecule_index(m)),R,0,sr); %[+ - +]
        scale_ind=[small_ind(1):small_ind(end)];
        
        diff_ind=setdiff(scale_ind,small_ind); %[- + -]
        
        %use only the first c60he peak
        last_two_ind=[diff_ind(1):scale_ind(end)]; %[- + +]
        small_ind_left=setdiff(scale_ind,last_two_ind); %[+ - -]
        small_ind_right=intersect(small_ind,last_two_ind); %[- - +]
        
        %find the scaling factor
        %a=M(diff_ind,:)\ref_peakdata(diff_ind,2); %such that M*a is close to the reference spec in the msd sense
        
        % apply background correction
        %single_spec=M*a-pchip(bgm,bgy,ref_peakdata(:,1));
        
        if m==10
            plot(ref_peakdata(small_ind,1),single_spec(small_ind),'r.',ref_peakdata(diff_ind,1),single_spec(diff_ind),'k.');
            %set(gca,'ylim',[0,2.5e5]);
            title(IFD_data{1}.molecules(molecule_index(m)).name)
            pause(0.1)
            sum(single_spec(small_ind));
        end
                
        output_data_left(i,(m+1)*2-1) = sum(single_spec(small_ind_left))/sum(single_spec(diff_ind));
        output_data_left(i,(m+1)*2)   = sqrt(1/sum(single_spec(small_ind_left))+1/sum(single_spec(diff_ind)))*output_data_left(i,(m+1)*2-1);
        
        output_data_right(i,(m+1)*2-1) = sum(single_spec(small_ind_right))/sum(single_spec(diff_ind));
        output_data_right(i,(m+1)*2)   = sqrt(1/sum(single_spec(small_ind_right))+1/sum(single_spec(diff_ind)))*output_data_right(i,(m+1)*2-1);
        
        %output_data(i,(m+1)*2)   = output_data(i,(m+1)*2-1)/sqrt(sum(M(small_ind,3)));
    end
end

fprintf('dlmwrite. please wait...');
    dlmwrite(scan_filename_left,output_data_left,'-append','delimiter','\t','precision','%e');
    dlmwrite(scan_filename_right,output_data_right,'-append','delimiter','\t','precision','%e');
fprintf(' done.\n');

