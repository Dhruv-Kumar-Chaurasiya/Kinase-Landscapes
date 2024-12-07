%% plot data
clear;clc;close all;
tic;
load FullVarData.mat;
%% Input
pdb='7M0K'; % Input PDB id (like 6XR6) or
name='MAP4K1'; % Input Kinase Name (like ABL1)
%% verification
idx=[];
if ~isempty(pdb) & isempty(name)
    idx=find(strcmpi(varsave(:,3),pdb));
    fam=varsave(idx,1);% family
    name=varsave(idx,2);% name
    state=varsave(idx,4);% state
    flag=1;
    if isempty(idx),disp('invalid or unused kinase PDB id');flag=0;end
elseif ~isempty(name) & isempty(pdb)
    idx2=find(strcmpi(varsave(:,2),name));
    if length(idx2)==2
        prp='Conformational State? A/I [A]:'; % A-active, I-inactive
        txt=input(prp,"s");
        if isempty(txt),txt="A";end
        if strcmpi(txt,"A") | strcmpi(txt,"Active")
            idx=idx2(1);
        elseif strcmpi(txt,"I") | strcmpi(txt,"Inactive")
            idx=idx2(2);
        end
    elseif length(idx2)==3
        prp='Conformational State of ABL1? A/I1/I2 [A]:';
        txt=input(prp,"s");
        if isempty(txt),txt="A";end
        if strcmpi(txt,"A") | strcmpi(txt,"Active")
            idx=idx2(1); % for 6XR6
        elseif strcmpi(txt,"I1") | strcmpi(txt,"Inactive1")
            idx=idx2(3); % for 6XR7 (Inactive1)
        elseif strcmpi(txt,"I2") | strcmpi(txt,"Inactive2")
            idx=idx2(2); % for 6XRG (Inactive2)
        end    
        if isempty(idx),disp('invalid or unused kinase PDB id');flag=0;end
    else
        idx=idx2;
    end
    fam=varsave(idx,1);% family
    pdb=varsave(idx,3);% PDB id
    state=varsave(idx,4);% state
    flag=1;
    if isempty(idx),disp('Invalid or unused kinase name');flag=0;end
elseif ~(isempty(pdb) & isempty(name))
    idx1=find(strcmpi(varsave(:,3),pdb));
    if isempty(idx1),idx1=0;end % defining any random fraction
    idx=idx1;flag=1;
    idx2=find(strcmpi(varsave(:,2),name));
    if isempty(idx2),idx2=0.005;end % defining any random fraction
    if sum(idx1==idx2)==0,disp('Invalid kinase name or PDB id.');flag=0;end
    if flag==1
        fam=varsave(idx,1);% family
        state=varsave(idx,4);% state
    end
else
    disp('Provide kinase name or PDB id.');flag=0;
end
%% plots
if flag~=0
    ene=varsave(idx,5);% ene-value
    block=varsave{idx,6};% blocks
    fes=varsave{idx,7};% fes data
    ResProb=varsave{idx,8};% ResProb
    Fpath=varsave{idx,9};% Fpath
    fes2D=varsave{idx,10};% fes2D
    dGc=varsave{idx,14};% dGc
    BlockDet=varsave{idx,15};% residue in blocks
    nmic=varsave{idx,16};% microstates sampled
    nres=BlockDet(end,1);
    P_r=ResProb(BlockDet(:,2));
    R=8.314/1000;T=310;% Universal gas constant and temprature used
    Gs_r=(-R*T)*log(P_r./(1-P_r));% residue folding stability
    dGc_r=dGc(BlockDet(:,2),BlockDet(:,2));
    mdGc_r=mean(dGc_r,'omitnan');% <dGc>
    % Free energy
    figure;
    plot(block,fes,'b','linewidth',1); axis([0 block(end) 0 ceil(max(fes))]);
    xlabel('Fraction of structured blocks');ylabel('FE (kJ mol^{-1})');
    xt = linspace(0,block(end),6); xt1 = cell(1,length(xt));
    for i=1:length(xt); xt1{i}=num2str(xt(i)/block(end)); end
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    title(strcat(state,": ",upper(name)," free energy profile"));
    ftxt=strcat("# microstates: ",num2str(nmic));
    ft=text(ceil(max(fes))/2,(ceil(min(fes)))/2,ftxt);
    ft.VerticalAlignment='bottom';ft.HorizontalAlignment='center';
    %residue folding probability across reaction cordinates
    figure;
    pcolor(Fpath); colormap jet; shading interp;
    xlabel('Fraction of structured blocks');
    ylabel('Residue Index');
    set(gca, 'XTick', xt, 'XTickLabel', xt1);box off;
    title(strcat(state,": ",upper(name)," folding probability"));
    % folding stability
    figure;
    plot(1:nres,Gs_r,'b',1:nres,zeros(1,nres),'k--','linewidth',1);
    axis([0 nres floor(min(Gs_r)) ceil(max(Gs_r))]);
    xlabel('Residue Index');ylabel('\DeltaG_s (kJ mol^{-1})');
    xt = round(linspace(0,nres,6)); xt1 = cell(1,length(xt));
    for i=1:length(xt); xt1{i}=num2str(xt(i)); end    
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    title(strcat(state,": ",upper(name)," residue folding stability"));
    % 2D free energy landscape
    figure;
    hv = round(block(end)/2);
    surf((0:block(end)-hv),(0:hv),fes2D); shading interp; colormap jet;
    xlabel('# of structured blocks in C-terminal');
    ylabel('# of structured blocks in N-terminal');
    zlabel('FE (kJ mol^{-1})');colorbar;view([138,19]);
    title(strcat(state,': ',upper(name),' free energy landscape'));
    % effective thermodynamic coupling (dGc)
    figure; 
    pcolor(dGc_r);shading interp;colormap jet;colorbar;
    xlabel('Residue index');ylabel('Residue index');
    title(strcat(state,": ",upper(name)," effective thermodynamic coupling"))
    % mean dGc
    figure;
    plot(1:nres,mdGc_r,'b',1:nres,zeros(1,nres),'k--','linewidth',1); 
    axis([0 nres floor(min(mdGc_r)) ceil(max(mdGc_r))]);
    xlabel('Residue Index');ylabel('<\DeltaG_c > (kJ mol^{-1})');
    xt = round(linspace(0,nres,6)); xt1 = cell(1,length(xt));
    for i=1:length(xt); xt1{i}=num2str(xt(i)); end    
    set(gca, 'XTick', xt, 'XTickLabel', xt1);    
    title(strcat(state,": ",upper(name)," <\DeltaG_c >"));
end
toc;