%% The Code Computes 1D free energy profile, 2D free energy surface and Residue Folding Probabilities as the function of RC
clear; clc; tic;
% Run after running cmapCalcElecBlock.m
%% Input Parameters
pdb=char();
pdb = '1BLX_sp_CMGC_CDK6'; % Input PDB file name
ene1 = -78.6/1000;% vdW interaction energy (kJ/mol) per native contact
DS = -14.5/1000; % Entropic cost (kJ/mol.K) per residue
DCp = -0.3579/1000; % Heat capacity change (kJ/mol.K) per native contact
T = 310; % Temperature in K
IS = 0.1; % Ionic strength in Molar units

%% Contact Maps obtained through PDB
aa=pdb(1,:);
eval(['load contactmapmatElecB',aa,'.dat;']);
eval(['cmapmask=contactmapmatElecB',aa,';']);
eval(['load contdistElecB',aa,'.dat;']);
eval(['contdist=contdistElecB',aa,';']);
eval(['load BlockSize',aa,'.dat;']);
eval(['bs=BlockSize',aa,';']);
eval(['load BlockDet',aa,'.dat;']);
eval(['BlockDet=BlockDet',aa,';']);
nres=length(cmapmask);
nres1 = length(unique(BlockDet(:,1)));
Mw=nres1*110;
wi=10;% for number of calculating windows in case of 'NaN's in the free energy projections
eval(['load disr',aa,'.dat;']);
eval(['disr=disr',aa,';']);% non-helical/non-strand/non-310 helix + Glycine residues positions
eval(['load ppos',aa,'.dat;']);
eval(['ppos=ppos',aa,';']);% proline residues positions

C=0;
%% Constants
R=0.008314;
Tref = 385;
zval=exp(DS./R);
zvalc=exp((DS-(6.0606/1000))./R); %% Excess entropic cost: dDS = -6.0606 J/mol.K, doi:- 10.1021/acs.jpcb.6b00658
zjj=zval.*ones(nres1,1);
zjj(disr)=zvalc;
zjj(ppos)=1;
zvec = ones(BlockDet(end,2),1);
for i=1:length(BlockDet)
    zvec(BlockDet(i,2)) = zvec(BlockDet(i,2)) * zjj(BlockDet(i,1));
end

%% calculate Zfin
ResProb = zeros(nres,length(T));
Zfin=zeros(length(T));
for kk=1:length(T)
    for ll=1:length(C)
        %% Elec. energy involving IS contribution
        ISfac=5.66*sqrt(IS/T(kk))*sqrt(80/29);%dielectric constant=29
        emapmask=zeros(nres);
        for i=1:nres
            for iin=i:nres
                x1=find(contdist(:,1)==i & contdist(:,2)==iin);
                emapmask(i,iin)=sum(contdist(x1,5).*exp(-ISfac.*contdist(x1,3)));
            end
        end
        pepval=zeros(nchoosek(nres+1,2)+2*nchoosek(nres+1,4),7);
        k=1;
        %% Generating Combinations for SSA
        for i=1:nres
            for iin=1:nres-i+1
                sw=0;
                stabEtemp=0; eneEtemp=0; nconttemp=0;
                wii=floor(iin/wi);
                swgen=zeros(wi,1);
                wistart=i;
                for wiii=1:wi-1
                    nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                    eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                    swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                    wistart=i+wii*wiii+1;
                end
                nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                sw=prod(swgen);
                
                pepval(k,:)=[(iin) sw i iin 0 0 1];
                k=k+1;
            end
            disp(i);
        end
        %% Generating Combinations for DSA
        for i=1:nres
            for iin=1:nres-i+1
                for j=i+iin+1:nres
                    for jin=1:nres-j+1    
                        %for island 1
                        sw1=0;
                        stabEtemp=0; eneEtemp=0; nconttemp=0;
                        wii=floor(iin/wi);
                        swgen=zeros(wi,1);
                        wistart=i;
                        for wiii=1:wi-1
                            nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                            swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                            wistart=i+wii*wiii+1;
                            stabEtemp=0; eneEtemp=0; nconttemp=0;
                        end
                        nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                        eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                        swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                        sw1=prod(swgen);
                        %for island 2
                        sw2=0;
                        stabEtemp=0; eneEtemp=0; nconttemp=0;
                        wii=floor(jin/wi);
                        swgen=zeros(wi,1);
                        wistart=j;
                        for wiii=1:wi-1
                            nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,j:j+jin-1)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,j:j+jin-1)));
                            swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                            wistart=j+wii*wiii+1;
                        end
                        nconttemp=sum(sum(cmapmask(wistart:j+jin-1,j:j+jin-1)));
                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                        eneEtemp=sum(sum(emapmask(wistart:j+jin-1,j:j+jin-1)));
                        swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                        sw2=prod(swgen);
                        
                        sw=sw1.*sw2;
                        
                        pepval(k,:)=[(iin+jin) sw i iin j jin 2];
                        k=k+1;
                    end
                end
            end
            disp([i j]);
        end
        %% Generating Combinations for DSAw/L
        for i=1:nres
            for iin=1:nres-i+1
                for j=i+iin+1:nres
                    for jin=1:nres-j+1
                        stabE=0; eneE=0; sw=0;
                        vv = [(i:i+iin-1) (j:j+jin-1)];
                        ncont=sum(sum(cmapmask(vv,vv)));
                        stabE=(ncont*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                        eneE=sum(sum(emapmask(vv,vv)));
                        if sum(sum(cmapmask(i:i+iin-1,j:j+jin-1)))~=0 || sum(sum(emapmask(i:i+iin-1,j:j+jin-1)))~=0
                            %from island 1
                            sw1=0;
                            stabEtemp=0; eneEtemp=0; nconttemp=0;
                            wii=floor(iin/wi);
                            swgen=zeros(wi,1);
                            wistart=i;
                            for wiii=1:wi-1
                                nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,vv)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,vv)));
                                swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                                wistart=i+wii*wiii+1;
                            end
                            nconttemp=sum(sum(cmapmask(wistart:i+iin-1,vv)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:i+iin-1,vv)));
                            swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                            sw1=prod(swgen);
                            
                            %from island 2
                            sw2=0;
                            stabEtemp=0; eneEtemp=0; nconttemp=0;
                            wii=floor(jin/wi);
                            swgen=zeros(wi,1);
                            wistart=j;
                            for wiii=1:wi-1
                                nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,vv)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,vv)));
                                swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                                wistart=j+wii*wiii+1;
                            end
                            nconttemp=sum(sum(cmapmask(wistart:j+jin-1,vv)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:j+jin-1,vv)));
                            swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                            sw2=prod(swgen);
                            
                            sw=sw1*sw2*zvalc^(j-(i+iin));
                        else
                            sw = 0;
                        end
                        pepval(k,:)=[(iin+jin) sw i iin j jin 3];
                        k=k+1;
                    end
                end
            end
            disp([i j]);
        end
        Zfin(kk) = sum(pepval(:,2));
        pepval = pepval(pepval(:,2)~=0,:);
        %% 1D free energy profile
        fes = zeros(nres,1);
        for i=1:nres
            fes(i) = sum(pepval(pepval(:,1)==i,2));
        end
        fes = fes./sum(fes);
        fes = -R*T(kk)*log(fes);
        
        subplot(2,1,1);
        plot((1:nres),fes,'b'); axis([0 inf 0 inf]);
        xlabel('# of structured blocks');
        ylabel('FE (kJ mol^{-1})');

        %% Residue folding probability as the function of RC and 2D free-energy profile
        hv = round(nres/2);
        Fpath = zeros(nres1,nres);
        FpathZ = zeros(nres,1);
        conv = zeros(nres,2);
        fes2D = zeros(hv+1,nres-hv+1);
        fes2DResProb = zeros(hv+1,nres-hv+1,nres);
        for i=1:nres
            f = find(BlockDet(:,2)==i);
            conv(i,:) = [BlockDet(f(1),1) BlockDet(f(end),1)];
        end
        for i=1:length(pepval)
            pept = zeros(nres1,1);
            pept(conv(pepval(i,3),1):conv(pepval(i,3)+pepval(i,4)-1,2))=1;
            if pepval(i,end)>1; pept(conv(pepval(i,5),1):conv(pepval(i,5)+pepval(i,6)-1,2))=1; end
            Fpath(:,pepval(i,1)) = Fpath(:,pepval(i,1)) + pepval(i,2)*pept;
            FpathZ(pepval(i,1)) = FpathZ(pepval(i,1)) + pepval(i,2);
            pept = zeros(nres,1);
            pept(pepval(i,3):pepval(i,3)+pepval(i,4)-1)=1;
            if pepval(i,end)>1; pept(pepval(i,5):pepval(i,5)+pepval(i,6)-1)=1; end
            ResProb(:,kk) = ResProb(:,kk) + pepval(i,2)*pept; % To calculate residue folding probability
            fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1) = fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1)+pepval(i,2);
            fes2DResProb(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1,:) = fes2DResProb(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1,:)+reshape(pepval(i,2)*pept,1,1,length(pept));
        end
        ResProb(:,kk) = ResProb(:,kk)./Zfin(kk);
        
        for i=1:nres
            Fpath(:,i) = Fpath(:,i)./FpathZ(i);
        end
            
        subplot(2,1,2);
        pcolor(Fpath); colormap jet; shading interp;
        xlabel('# of structured blocks');
        ylabel('Residue Index');

        %% 2D free energy profile and residue folding probability at each point
        for i=1:length(fes2D(:,1))
            for j=1:length(fes2D(1,:))
                fes2DResProb(i,j,:) = fes2DResProb(i,j,:)./fes2D(i,j);
            end
        end
        % To visualize residue folding probability at (n(C),n(N)) :
        % plot(reshape(fes2DResProb(n(N)+1,n(C)+1,:),nres,1));
        fes2D = fes2D./sum(sum(fes2D));
        fes2D = -R*T(kk)*log(fes2D);
        figure; surf((0:nres-hv),(0:hv),fes2D); shading interp; colormap jet;
        xlabel('# of structured blocks in C-terminal');
        ylabel('# of structured blocks in N-terminal');
        zlabel('FE (kJ mol^{-1})')
    end
end
%% coupling
pepvalx=pepval; clear pepval;
npert=1;
x=squeeze(pepvalx(:,:,1));
nmic=length(squeeze(pepvalx(:,1,1)));
pv=zeros(nmic,npert);
respa=zeros(nres,npert);

respn=zeros(nres,npert);
pjfkf=zeros(nres,nres,npert);
pjfku=zeros(nres,nres,npert);
pjukf=zeros(nres,nres,npert);
pjuku=zeros(nres,nres,npert);
for kk=1:npert
    Z(kk,1)=sum(squeeze(pepvalx(:,2,kk)));
    pv(:,kk)=squeeze(pepvalx(:,2,kk))./Z(kk);        
    for j=1:nres           
        % index for all states in which res j is structured
        a=find( (x(:,4)+x(:,3)-1)>=j & x(:,3)<=j | ((x(:,6)+x(:,5)-1)>=j & x(:,5)<=j) );
        % probability of states where res j is structured
        respa(j,kk)=sum(squeeze(pepvalx(a,2,kk)))./Z(kk); 
        
        n=setdiff(1:length(x),a,'stable'); % index for all states in which j in unstructured
        % probability of states where res j is structured
        respn(j,kk)=sum(squeeze(pepvalx(n,2,kk)))./Z(kk); 
        
        for k=1:nres
            % index for all states in which res k is structured
            a1=find( (x(:,4)+x(:,3)-1)>=k & x(:,3)<=k | ((x(:,6)+x(:,5)-1)>=k & x(:,5)<=k));
            % probability for all states in which res j and k are structured
            pjfkf(j,k,kk)=sum(pv(intersect(a,a1,'stable'),kk));
            
            % index for all states in which res k is unstructured
            n1=setdiff(1:length(x),a1,'stable');
            % index of all states in which res j is structured but k is not
            a2=intersect(a,n1,'stable');
            % probability for all states in which res j is structured but k is not
            pjfku(j,k,kk)=sum(pv(a2,kk));
            
            % index for all states in which res j is unstructured but k is structured
            a3=intersect(n,a1,'stable');
            % probability for all states in which res j is unstructured but k is structured
            pjukf(j,k,kk)=sum(pv(a3,kk));
            
            % index for all states in which res j and k are unstructured
            a4=intersect(n,n1,'stable');
            % probability for all states which res j and k are unstructured
            pjuku(j,k,kk)=sum(pv(a4,kk));
            clear a1 a2 a3 a4 n1
            [j k]
        end
    end
end
dGcp=R*T*log( (pjuku(:,:,end).*pjfkf(:,:,end))./ (pjfku(:,:,end).*pjukf(:,:,end)));%dGc
chiplus=R*T*log(squeeze(pjfkf(:,:,end))./(pjukf(:,:,end))); %dGplus
chiminus=R*T*log(squeeze(pjfku(:,:,end))./(pjuku(:,:,end))); %dGminus
% dGplus or symetric chi_plus
Gpsym=(chiplus+chiplus')./2;
%dGminus or symetric chi_minus
Gmsym=(chiminus+chiminus')./2;
a3=R*T*log((pjfkf(:,:,end)+pjuku(:,:,end))./(pjfku(:,:,end)+pjukf(:,:,end)));
%infinity replacement by NaN
Gpsym(isinf(Gpsym))=NaN; Gmsym(isinf(Gmsym))=NaN; a3(isinf(a3))=NaN;
dGcp(isinf(dGcp))=NaN;
figure; pcolor(dGcp);shading interp;colormap jet;colorbar;
xlabel('Residue index');
ylabel('Residue index');
save([aa,'FECoupdata.mat']);
toc;