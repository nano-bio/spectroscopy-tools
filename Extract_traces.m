%Extracts Traces with singnal normalized to the sum spectrum of IFD file 1

%IFD file calibration: LOAD H5 WITHOUT INTERNAL CALIBRATION!!!!

%IsotopeFit Path:
addpath('Z:\Experiments\STM\matlab\IsotopeFit\');

%folder
folder='Z:\Temp\rALSER\';

%IFD files:
IFD_file=[];
IFD_file{1}=[folder,'Bereich 3 fein 1 960nm.ifd'];
IFD_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm.ifd'];

%H5 files:
H5_file=[];
H5_file{1}=[folder,'Bereich 3 fein 1 960nm.h5'];
H5_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm.h5'];

%Wavelength files
lambda_file=[];
lambda_file{1}=[folder,'Bereich 3 fein 1 960nm - Wellenlaengen.txt'];
lambda_file{2}=[folder,'Bereich 3 Abschnitt 2 fein 960nm - Wellenlaengen.txt'];


%Export filename
scan_filename=[folder,'Bereich 3 fein ENERGIES_He100.txt'];

%List of molecules:
n_He=1:100; %number of C60-Helium traces to extract

%searchrange (for scaling)
sr=10; %look sr au above/below the current molecule


molecules=[];
he_series=[];
for i=1:length(n_He)
    if n_He(i)==1
        molecules{i}='[C60][He]';
    else
        molecules{i}=sprintf('[C60][He]%i',n_He(i));
    end
    he_series{i}=sprintf('[He]%i',n_He(i)+180);
end

if length(molecules)~=length(he_series)
    fprintf('PANIC!!! Not all molecules found! I abort!\n');
    break
end

energy_axis=[];
n_bufs=[];
sort_idx={};
for i=1:length(lambda_file)
    temp=load(lambda_file{i});
    [temp,sort_idx{i}]=sort(temp);
    n_bufs(i)=length(temp);
    energy_axis=[energy_axis;temp];
end

ES_mat=zeros(length(energy_axis),length(molecules));

position=1;

%write title line
fid=fopen(scan_filename,'w');
fprintf(fid,'Energy\t');

for i=1:length(molecules)
    fprintf(fid,'%s\t',molecules{i});
end

fprintf(fid,'\n');
fclose(fid);

for f=1:length(IFD_file)
    data={};
    load(IFD_file{f},'-mat');
    
    %find the indices of the molecules of interest
    index=find(ismember({data.molecules.name},molecules));
    index_he=find(ismember({data.molecules.name},he_series));
    
    %n_writes=getnumberofinstancesinh5(H5_file{f},'writes')
    %n_bufs=getnumberofinstancesinh5(H5_file{f},'buffers')
           
    data.peakdata=subtractmassoffset(data.raw_peakdata,data.calibration);
    
    %load the reference sum spectrum?
    if f==1
         ref_peakdata=data.peakdata;
         ref_cal=data.calibration;
         %minmass=min([data.molecules(index).minmass])-sr;
         %maxmass=max([data.molecules(index).maxmass])+sr;
         minmass=750;
         maxmass=1200;
         
    end
    
    %scale the buffers: which indices?
    ref_ind_buf=mass2ind(data.peakdata(:,1),minmass):mass2ind(data.peakdata(:,1),maxmass);
    
    %where is the data in the reference spec?
    ref_ind=mass2ind(ref_peakdata(:,1),minmass):mass2ind(ref_peakdata(:,1),minmass)+length(ref_ind_buf)-1; %make sure the two ranges have the same length
        
    %scaling: dont look at the peaks which are affected by the laser
    for i=1:length(index)
        R=resolutionbycalibration(ref_cal,data.molecules(index(i)).com);
        [ref_ind_buf,ix]=setdiff(ref_ind_buf,findmassrange2(data.peakdata(:,1)',data.molecules(index(i)),R,0,1));
        ref_ind=ref_ind(ix);
    end
    
    %for w=1:n_writes
    w=1;
    plot_spec=[];
        for b=sort_idx{f}'
            single_spec=readh5buffer(H5_file{f}, w, b);
            if mod(position,10)==0, fprintf('%f %%\n',100*position/length(energy_axis)); end;
            
            %find the scaling factor
            %design matrix for scaling: (constant baseline, linear baseline,
            %data)
            M=[ones(length(ref_ind_buf),1),data.peakdata(ref_ind_buf,1),single_spec(ref_ind_buf)'];
            %find the scaling factor
            a=M\ref_peakdata(ref_ind,2); %such that M*a is close to the reference spec in the msd sense
            
            %scale the data
            single_spec=[ones(length(single_spec),1),data.peakdata(:,1),single_spec']*a;
            
            for i=1:length(index)
                %a=single_spec(ref_ind_buf)'\ref_peakdata(ref_ind,2);
                
                R=resolutionbycalibration(ref_cal,data.molecules(index(i)).com);
                                
                %ES_mat(position,i)=sum(a*single_spec(small_ind));
                ES_mat(position,i)=sum(single_spec(findmassrange2(data.peakdata(:,1)',data.molecules(index(i)),R,0,0.5)));
            end
            position=position+1;
        end
    %end
   
end
 %write ASCII file
    fprintf('dlmwrite. please wait...');
    dlmwrite(scan_filename,[energy_axis,ES_mat],'-append','delimiter','\t','precision','%e');
    fprintf(' done.\n');