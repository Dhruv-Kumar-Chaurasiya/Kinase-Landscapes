%% Computes Weighted Contact Map, Pairwise Electrostatic Energy between Charged Atoms and Defined Block Approximation based on the Protein Size
clear; clc; tic;
tline='5UY6_sp_Other_CAMKK2';% PDB filename
%% Input Parameters
pdb = tline;
cdata=cat(2,deblank(pdb(1,:)),'.pdb');
pdbr=fopen(cdata,'rt');
pdbc=fopen('buffer1.txt','wt');
pH = 7; % pH value
srcutoff=5.0; % Cutoff distance for VdW interaction (in Angstrom)
BlockSize = 3; % can be manually set to a specific value. A heuristic for blocksize is also provided in lines 132-139
ecutoff=1000; % Cutoff dist for electrostatic interaction (in Angstrom)
seqsepcutoff = 3; % cut-off for local contacts (here, j-i<3 for local contacts)
% In linux uncomment 2 lines below 
%command = cat(2,'./stride ',cdata,' -fstruct',pdblist,'.txt');
%system(command); 
structfile = cat(2,'struct',pdb,'.txt'); % file with sec. str. information
%% Remove Hydrogens and hetero atoms including nucleotides
line=fgetl(pdbr);
while feof(pdbr)~=1
    if(length(line)>4)
        if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0)
            fprintf(pdbc,'%s\n',line);
        end
        if(strcmp(line(1:3),'TER')==1)
            break;
        end
    end
    line=fgetl(pdbr);
end
fclose(pdbr);
fclose(pdbc);
    
%% Get atmoic details with their corresponding residue number
% atomn: Atom name, resno: Residue number, x,y,z: coordinates
i=1; j=1; prevres=0; nnegres=0; nposres=0; nhphilres=0; iin=1; iim=1;
aares={'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER','ASP','ASN','THR','GLU','GLN','HIS','LYS','ARG','PRO','CYS'};
aacode={'G','A','V','L','I','M','F','Y','W','S','D','N','T','E','Q','H','K','R','P','C'};
charres={'HIS','LYS','ARG','GLU','ASP'};
posres=[16 17 18]; % [HIS LYS ARG]
negres=[11 14]; % [ASP GLU]
hphilres=[8 10 12 13 15 20]; % [TYR SER ASN THR GLN CYS]
atomc={'NE ','NH1','NH2','NZ ','OD1','OD2','OE1','OE2','ND1','NE2'};
atombb={'C','N','CA'};
charmag7=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0 0]'; % pH 7
charmag5=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0.5 0.5]'; % pH 5
charmag3p5=[0.33 0.33 0.33 1 -0.25 -0.25 -0.25 -0.25 0.5 0.5]'; % pH 3.5
charmag2=[0.33 0.33 0.33 1 0 0 0 0 0.5 0.5]'; % pH 2   
pdbc=fopen('buffer1.txt','rt');
    
while feof(pdbc)~=1
    line=fgetl(pdbc);
    if(length(line)>4)
        atomn(i,:)=line(14:16);
        resno(i,1)=str2double(line(23:26));
        x(i,1)=str2double(line(32:38));
        y(i,1)=str2double(line(40:46));
        z(i,1)=str2double(line(48:54));
        
        % Assign charge magnitude at specified pH
        charmag(i,1)=0;
        if sum(strcmp(charres,line(18:20)))
            [~,x1] = find(strcmp(atomc,atomn(i,:)));
            if isempty(x1)==0
                if pH == 7
                    charmag(i,1)=charmag7(x1);
                elseif pH == 5
                    charmag(i,1)=charmag5(x1);
                elseif pH == 3.5
                    charmag(i,1)=charmag3p5(x1);
                elseif pH == 2
                    charmag(i,1)=charmag2(x1);
                end
            end
        end
        
        % Calculate Seq Composition
        if (resno(i)-prevres~=0)
            % Get sequence from PDB
            [~,x2]=find(strcmp(line(18:20),aares));
            protseq(1,j)=aacode(x2);
            prevres=resno(i);
            j=j+1;
        end
    end
    i=i+1;
end
fclose(pdbc);
%% Compute Short and Long-Range interaction map
% Short_Range: no. of atomic contacts and h-bonds b/w residues
% Long_Range: details for electrostatic interactions
%-----Initializarion---------------------------------------------------
atomt=length(x(:,1));
nres=resno(atomt,1)-resno(1,1)+1;
srcont=zeros(nres);
loccont=zeros(nres); locconte=zeros(nres);
nloccont=zeros(nres); nlocconte=zeros(nres);
%----------------------------------------------------------------------
startres=resno(1,1);
resno(:,1)=resno(:,1)-startres+1; %Changing res. no. according to the PDB file
k=1;
for i=1:(atomt-1)
    for j=(i+1):atomt
        dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
        % Calculate Short-Range Interactions
        if((resno(j)-resno(i))>=1 && (charmag(i)==0 || charmag(j)==0))% for j>i in case of nearest neighbour
            dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
            if(dist<=srcutoff)
                srcont(resno(i),resno(j))=srcont(resno(i),resno(j))+1;                
            end
        end
        % Calculate Long-Range Interactions and associated Energy
        if ((resno(j)-resno(i))>=1 && (charmag(i)~=0 && charmag(j)~=0))
            diste=sqrt(((x(i)-x(j))^2)+((y(i)-y(j))^2)+((z(i)-z(j))^2));
            if(diste<=ecutoff)
                intene=332*4.184*charmag(i)*charmag(j)/(29*diste); %Elec. Contibution with dielectric 29
                % Elec. Info in matrix form
                ElecMat(k,:)=[resno(i) resno(j) diste abs(resno(i)-resno(j)) intene];
                k=k+1;
            end
        end
    end
end
% Heuristic to choose the block size % given in JMB (2021) doi;- https://doi.org/10.1016/j.jmb.2021.167325
% if nres<=300
%     BlockSize = ceil(nres/100);
% elseif nres<600
%     BlockSize = 5;
% else
%     BlockSize = 6;
% end   
%% Get Secondary Structure(SecStr) Info
sfile=fopen(structfile,'rt'); i=1; j=1; disr=[]; ppos=[]; %structure file
lsfile=fgetl(sfile);
while lsfile > 0
    if length(lsfile)> 25 && strcmp(lsfile(1:3),'ASG')
        h(i,1) = 0;
        %H,G,E are alpha,3-10 and extended conformations, used here
        if (strcmp(lsfile(25),'H') || strcmp(lsfile(25),'E') || strcmp(lsfile(25),'G'))
            h(i,1) = 1;
        end
        i=i+1; % increasing helix row vector
        if strcmp(lsfile(6:8),'GLY') || (strcmp(lsfile(25),'H')==0 && strcmp(lsfile(25),'E')==0 && strcmp(lsfile(25),'G')==0)
            disr=[disr str2double(lsfile(18:20))];% non-helical, non-strand and glycine residue position               
        end
        if strcmp(lsfile(6:8),'PRO')
            ppos = [ppos str2double(lsfile(18:20))]; % proline residue position
        end
    end
    lsfile=fgetl(sfile);
end
fclose(sfile);   
    
h1=[[0;h] [h;0]];
hbeg = find(h1(:,1)==0 & h1(:,2)==1);
hend = find(h1(:,1)==1 & h1(:,2)==0);
StructBlock = startres; hbeg = hbeg+startres-1; hend=hend+startres-1;
for i=1:length(hbeg)
    StructBlock = [StructBlock;hbeg(i);hend(i)];
end
StructBlock = [StructBlock;resno(end)+startres];
residual = rem(StructBlock(2:end)-StructBlock(1:end-1),BlockSize);
nBlocks = floor((StructBlock(2:end)-StructBlock(1:end-1))/BlockSize);
k=1; BlockID = 0; uresno = unique(resno);
BlockMat = zeros(length(uresno),2); flag=0;
for i=2:length(StructBlock)
    for j=1:nBlocks(i-1)
        BlockID = BlockID+1;
        BlockMat(k:k+BlockSize-1,:) = [uresno(k:k+BlockSize-1) BlockID*ones(BlockSize,1)];
        k=k+BlockSize; flag=1;
    end
    if residual(i-1) == 1
        if flag==0; BlockID = BlockID+1; end
        BlockMat(k,:) = [uresno(k) BlockID];
        k=k+1;
    elseif residual(i-1) > 1
        BlockID = BlockID+1;
        BlockMat(k:k+residual(i-1)-1,:) = [uresno(k:k+residual(i-1)-1) BlockID*ones(residual(i-1),1)];
        k=k+residual(i-1);
    end
end
    
contactMapB = zeros(BlockID);
for i=1:nres
    for j=i+1:nres
        contactMapB(BlockMat(i,2),BlockMat(j,2)) = contactMapB(BlockMat(i,2),BlockMat(j,2)) + srcont(i,j);
    end
end
contdistMapB = zeros(length(ElecMat(:,1)),5);
for i=1:length(ElecMat(:,1))
    contdistMapB(i,:) = [BlockMat(ElecMat(i,1),2) BlockMat(ElecMat(i,2),2) ElecMat(i,3:end)];
end
%% Printing Regular Contact Map files
eval(['save contactmapmatElec',pdb,'.dat srcont -ascii']);
eval(['save contdistElec',pdb,'.dat ElecMat -ascii']);
eval(['save contdistElecB',pdb,'.dat contdistMapB -ascii']);
eval(['save contactmapmatElecB',pdb,'.dat contactMapB -ascii']);
eval(['save BlockDet',pdb,'.dat BlockMat -ascii']);
eval(['save BlockSize',pdb,'.dat BlockSize -ascii']);
eval(['save disr',pdb,'.dat disr -ascii']);
eval(['save ppos',pdb,'.dat ppos -ascii']);
toc;
