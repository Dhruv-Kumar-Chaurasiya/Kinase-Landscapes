%% plot Coupling data of Alanine-scanning mutagenesis 
clear;clc;close all;
tic;
load AlaScanData.mat;
%% Input
% only PRKACA, PIM1, CSNK1A1, MAPK1, NEK1, ULK3, MAP4K1, RIPK2, ABL1 as name
% only 3OVV, 1XWS, 6GZD, 3SA0, 4APC, 6FDY, 7M0M, 6S1F, 6XR6 as PDB
% only AGC, CAMK, CK1, CMGC, NEK, STE, TKL, TYR, Other as family
pdb='3OVV'; 
name='';
fam='';
mutpos=1; % residue position undergoing alanine mutation
%% verification
idx=[];
if ~isempty(pdb) & isempty(name) & isempty(fam)
    idx=find(strcmpi(data(:,3),pdb));
    fam=data(idx,1);% family
    name=data(idx,2);% name
    flag=1;
    if isempty(idx),disp('invalid or unused kinase PDB id');end
elseif ~isempty(name) & isempty(pdb) & isempty(fam)
    idx=find(strcmpi(data(:,2),name));
    fam=data(idx,1);% family
    pdb=data(idx,3);% PDB id
    flag=1;
    if isempty(idx),disp('Invalid or unused kinase name');flag=0;end
elseif ~isempty(fam) & isempty(pdb) & isempty(name)
    idx=find(strcmpi(data(:,1),fam));
    name=data(idx,2);% Name
    pdb=data(idx,3);% PDB id
    flag=1;
    if isempty(idx),disp('Invalid or unused kinase family');flag=0;end
elseif ~(isempty(pdb) & isempty(name)) & isempty(fam)
    idx1=find(strcmpi(data(:,3),pdb));
    if isempty(idx1),idx1=0;end % defining any random fraction
    idx=idx1;flag=1;
    idx2=find(strcmpi(data(:,2),name));
    if isempty(idx2),idx2=0.005;end % defining any random fraction
    if sum(idx1==idx2)==0,disp('Invalid kinase name or PDB id');flag=0;end
    if flag==1
        fam=data(idx,1);% family
    end
elseif ~(isempty(pdb) & isempty(fam)) & isempty(name)
    idx1=find(strcmpi(data(:,3),pdb));
    if isempty(idx1),idx1=0;end % defining any random fraction
    idx=idx1;flag=1;
    idx2=find(strcmpi(data(:,1),fam));
    if isempty(idx2),idx2=0.005;end % defining any random fraction
    if sum(idx1==idx2)==0,disp('Invalid kinase family or PDB id');flag=0;end
    if flag==1
        name=data(idx,2);% name
    end
elseif ~(isempty(name) & isempty(fam)) & isempty(pdb)
    idx1=find(strcmpi(data(:,2),name));
    if isempty(idx1),idx1=0;end % defining any random fraction
    idx=idx1;flag=1;
    idx2=find(strcmpi(data(:,1),fam));
    if isempty(idx2),idx2=0.005;end % defining any random fraction
    if sum(idx1==idx2)==0,disp('Invalid kinase name or family');flag=0;end
    if flag==1
        pdb=data(idx,3);% PDB
    end
elseif ~(isempty(name) & isempty(fam) & isempty(pdb))
    idx1=find(strcmpi(data(:,3),pdb));
    if isempty(idx1),idx1=0;end % defining any random fraction
    idx=idx1;flag=1;
    idx2=find(strcmpi(data(:,2),name));
    if isempty(idx2),idx2=0.005;end % defining any random fraction
    idx3=find(strcmpi(data(:,1),fam));
    if isempty(idx3),idx2=0.001;end % defining any random fraction
    if isequal(idx1,idx2,idx3)==0,disp('Invalid kinase family, name or pdb id');flag=0;end    
else
    disp('Provide either kinase name, family or PDB id.');idx=0;flag=0;
end
%% plots
if flag~=0
    mutfrom=string(data{idx,4}); %mutation
    mutnum=data{idx,5}; % mutated residue number
    CA_D=data{idx,6};% C_a distance 
    feswt=data{idx,7};% wild-type fes data
    fesA=data{idx,8};% all mutated fes data
    dGpwt=data{idx,9};% wild-type dG+
    dGpA=data{idx,10};% all mutated dG+
    DCI=data{idx,11};% DCI obtained using (DCM) ddG+=dG+,mut-dG+,wt
    xfit=data{idx,12};% fit parameters
    ci=data{idx,13};% fit parameter confidence interval
    yfit=data{idx,14};% y-fit values
    xm=data{idx,16};% mean of CA_dist at each 2 Angstrom 
    xs=data{idx,15};% std. dev. of CA_dist at each 2 Angstrom 
    ym=data{idx,17};% mean DCI at each 2 Angstrom CA_dist  
    ys=data{idx,18};% std. dev. DCI at each 2 Angstrom CA_dist  
    idmut=find(mutpos==mutnum);
    if ~isempty(idmut)
        mut=strcat(mutfrom(idmut),num2str(mutnum(idmut)),"A");
        % Free energy
        fx=linspace(1/size(feswt,1),1,size(feswt,1))';
        fesmut=fesA(:,idmut);
        figure;
        plot(fx,feswt,'b',fx,fesmut,'k--','linewidth',1);xlim([0,1]);
        xlabel('Fraction of structured blocks');ylabel('FE (kJ mol^{-1})');
        xt = linspace(0,fx(end),6); xt1 = cell(1,length(xt));
        for i=1:length(xt); xt1{i}=num2str(xt(i)/fx(end)); end
        set(gca, 'XTick', xt, 'XTickLabel', xt1);
        legend("WT",mut);
        title(strcat(upper(fam),": ",upper(name)," (",mut,") free energy profile"));
        % MR vs Dist
        MR=DCI(:,idmut);
        Dx=CA_D(:,idmut);
        xft=round(xfit(idmut,:),2);
        Dxm=xm{idmut};Dxs=xs{idmut};
        MRfit=yfit{idmut};MRm=ym{idmut};MRs=ys{idmut};
        ub=round(max(max(MRfit),max(MRm+MRs)),1);
        figure;
        plot(Dx,MR,'bo','linewidth',1); 
        hold on; errorbar(Dxm,MRm,Dxs,'ko','horizontal');errorbar(Dxm,MRm,MRs,'ro');
        plot(Dxm,MRfit,'r-','linewidth',1); hold off;
        axis([0 50 0 (ub+0.5)]);
        xlabel(['Distance from mutated site (',char(197),')']);
        ylabel('|\Delta\DeltaG_+| (kJmol^{-1})');
        xt = linspace(0,50,6);
        xt1 = cell(1,length(xt));for i=1:length(xt);xt1{i}=num2str(xt(i));end
        set(gca,'XTick', xt, 'XTickLabel', xt1);
        txt=['y = (' num2str(xft(1)) ') e^{-x/(' num2str(xft(2)) ')} + ' num2str(xft(3))];
        text(45,(ub+0.2),txt,'HorizontalAlignment','right','VerticalAlignment','top','fontsize',12)
        title(strcat(upper(fam),": ",upper(name)," (",mut,") mutational response"));
    else
        disp('Invalid mutation positions (refrain from selecting Alanine, Glycine, Proline residue positions).');
    end
end
toc;