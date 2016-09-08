% Polarimetric Biological Radar Simulator
% Version 4
%
% Updated May 16, 2012:
%       - Calculate velocity by difference in position
%       - Handle errors when scatterers are beyond radar range
% 
% Updated December 21, 2012:
%       - Added 4/3 Earth refraction
%       - Added conditional polarization statements
%
% Updated June 26, 2013:
%       - Changed book keeping to store data as long array and only 
%       - construct sweeps after volume scan is complete. 
%
% Updated July 7, 2014:
%       - Added calculation of orientation angles. 
%
% Updated November 21, 2014:
%       - Added system phase offset and depolarized contributions.        


% clear all
% close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%  Define Paths and Files  %%%%%%%%%%%%%%%%%%%%%%%%

% biological behavior model input file name
%bioFileName = '/Users/phillipstepanian/backupRoot/research/bioSim/data/ieeeMig/chilPDF/100perSqkmCo/100perSqkmCo0001.mat';
bioFileName = 'science1.mat';
% scatter model look-up table file name
sctFileName = 'batScatter.mat';

% output IQ file save name
iqSaveName  = 'polBat4IQ_30';

% to disable user I/O interaction, set to '0'
interactiveIO = 1;
  

%%%%%%%%%%%%%%%%%%%%%%%  Define Physical Constants  %%%%%%%%%%%%%%%%%%%%%%%

cc        = 3e8;            % speed of light         [m/s]
ii        = sqrt(-1);         % imaginary unit
eta0        = 377;            % impedance of free space  [ohms]
k_boltz     = 1.38e-23;         % Boltzmann's constant     [J/degree]
r_earth     = 6378.1e3;         % Radius of Earth        [m]
k_e         = 4/3;            % Effective Earth radius factor (D-Z p.21)


%%%%%%%%%%%%%%%%%%%%%%  Set Radar System Attributes  %%%%%%%%%%%%%%%%%%%%%%

lam         = 0.107;          % radar wavelength        [m]
bmWid       = .96;            % 3-dB beam width         [degrees]
PRT         = 1e-3;           % pulse repetition time     [s]
pulsWid     = 1.57e-6;        % pulse width             [s]
sampTm      = 1.57e-6;        % time between rx samples   [s]
gain        = 45;             % antenna gain            [dBi]
txPow       = 1000e3;         % peak transmit power       [W]
rxBndWid    = 0.63e6;         % receiver bandwidth        [Hz]
nsTemp      = 290;            % system temperature        [K]
rxNsFig     = 3;              % receiver noise figure     [dB]
sysDP       = 30;              % system diff phase (v ahead of h) [deg]


%%%%%%%%%%%%%%%%%%%%%%%%%%  Set Scan Attributes  %%%%%%%%%%%%%%%%%%%%%%%%%%

polType  = 'STSR';            % polarization type        [HORZ VERT STSR]
scanRate = 3;                 % pedestal scan rate       [degrees/sec]
elSet    = [0.5];             % elevation angles         [degrees]
dAz      = 1;                 % azimuthal steps        [degrees]
azSet    = 1:dAz:360;         % azimuthal angles         [degrees]
nVols    = 5;                 % number of volume scans   [#]


%%%%%%%%%%%%%%%%  Define Biological Model Reference Frame  %%%%%%%%%%%%%%%%

% The biological model is set in a reference frame with respect to (0,0,0).
% We will assume the radar location is at (0,0,0), so to have biology
% located away from the radar, the position of the biological reference
% frame origin will be moved to a position (bOrX,bOrY,bOrZ) with respect to 
% the radar reference frame origin.

bOrX = 10000;                    % biological x origin w.r.t radar     [m] 
bOrY = 10000;                    % biological y origin w.r.t radar     [m] 
bOrZ = 0;                        % biological z origin w.r.t radar     [m] 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End User Input  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%  Calculate Dependent Scan Parameters  %%%%%%%%%%%%%%%%%%

freq   = cc/lam;                       % radar frequency        [Hz]
Va     = lam/4/PRT;                    % aliasing velocity      [m/s]
ra     = PRT*cc/2;                     % max unambiguous range  [m]
rWid6  = cc*pulsWid/2;                 % 6-dB range width       [m]
azRes  = azSet(2)-azSet(1);              % azimuthal resolution   [degrees]
nGates = floor((PRT-pulsWid)/sampTm);      % number of range gates  [#]
nAzs   = length(azSet);                % number of azimuths     [#]
nEls   = length(elSet);                % number of elevations   [#]
rngSet = cc*(pulsWid+sampTm*(0:nGates-1))/2; % range gates          [m] 


% find pulses per dwell [#] and modify scan rate to keep integer pulses
if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
    dwellNum    = floor(1/scanRate*azRes/PRT); 
    scanRate    = 1/(dwellNum*PRT/azRes);
elseif strcmp(polType,'ATAR')
    dwellNum    = floor(1/scanRate*azRes/PRT/2);
    scanRate    = 1/(dwellNum*PRT/azRes*2);
end

% volume scan update time [s]
volTime     = 1/scanRate*(azSet(end)-azSet(1)+dAz)*length(elSet);


% save run header data
header.wavelength  = lam;
header.bmWid       = bmWid;      
header.PRT         = PRT;       
header.pulsWid     = pulsWid;   
header.sampTm      = sampTm;  
header.gain        = gain;     
header.txPow       = txPow;      
header.rxBndWid    = rxBndWid;   
header.nsTemp      = nsTemp;    
header.rxNsFig     = rxNsFig;          
header.offSet      = [bOrX bOrY bOrZ];
header.polType     = polType;
header.scanRate    = scanRate;


%%%%%%%%%%%%%%%%%%%%%%%%%  Display Run Attributes  %%%%%%%%%%%%%%%%%%%%%%%%


if interactiveIO
    fprintf('behavior model input filename: %s\n',bioFileName)
    fprintf('scatter model input filename:  %s\n',sctFileName)
    fprintf('output IQ save filename: %s\n',iqSaveName)
    fprintf('wavelength \t\t= %g m\nfrequency \t\t= %g Hz\nPRT \t\t\t= %g s\n',lam,freq,PRT)
    fprintf('6-dB range resolution \t= %g m\n3-dB beam width \t= %g degrees\n',rWid6,bmWid)
    fprintf('azimuthal sampling \t= %g degrees \naliasing velocity \t= %g m/s\n',azRes,Va)
    fprintf('max unambiguous range \t= %g m\npulses per dwell \t= %g\n',ra,dwellNum)
    fprintf('scan rate \t\t= %g degrees/s\nvolume update time \t= %g s\n',scanRate,volTime)
    fprintf('dual-polarization type \t= %s\n',polType)
    fprintf('bio reference location \t= (%g,%g,%g)\n',bOrX,bOrY,bOrZ)
    
    response = input('\nContinue with simulation? (n/y): ','s');
    if ~strcmp(response,'y')
        return
    end
    clear response
end


%%%%%%%%%%%%%%%%%%%%%%%%%%  Begin Simulation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the scatter definitions
load(sctFileName)

% find azimuth for each pulse in volume scan
azTm = (azSet(1)-bmWid/2)+scanRate*PRT*(1:(dwellNum*nAzs));
azTm(azTm>=360) = azTm(azTm>=360)-360;
azTm = single(repmat(azTm,[1 nEls]));

% find elevation for each pulse in volume scan
elTm = single(kron(elSet,ones(1,dwellNum*nAzs)));


% loop over and save each volume scan
for volCtr = 1:nVols
    
    
%     preallocate space for time series data for this volume
    if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
        Vh = single(nan(nAzs*dwellNum*nEls,nGates));
    end
    if strcmp(polType,'STSR')||strcmp(polType,'VERT')
        Vv = single(nan(nAzs*dwellNum*nEls,nGates));
    end
    
    
%     load biological behavior model input for this volume scan
    if interactiveIO,fprintf('Loading bio model input: %s\n',bioFileName), end
    load(bioFileName); % [m]
    
    agntPos = batLoc;
    
    agntVel = diff(agntPos,1,3);
    agntVel = cat(3,agntVel,agntVel(:,:,end));
    idx            = str2double(bioFileName(find(bioFileName=='.')-1));
    if ceil(volCtr*volTime)+1 > idx*size(agntPos,3)*dt
        %% at end of bio model file, check for others
        if interactiveIO,disp('Passed end of current bio model file. Looking for extensions.'), end
        idx            = find(bioFileName=='.')-1;
        bioFileName(idx) = num2str(str2double(bioFileName(idx))+1);
        idx            = dir(bioFileName);
        if isempty(idx)
          if interactiveIO, disp('No additional file found. Ending program.'), end
          break
        else
          if interactiveIO,fprintf('Additional file found.\nLoading bio model input: %s\n',bioFileName), end
          pad    = mod(ceil(volCtr*volTime)+1 - size(agntPos,3),size(agntPos,3));
          bioPos = single(agntPos(:,:,floor((volCtr-1)*volTime)+1:end));
          bioVel = single(agntVel(:,:,floor((volCtr-1)*volTime)+1:end));
          clear agntPos agntVel dt idx
          load(bioFileName,'agntPos','dt')
          agntVel = diff(agntPos,1,3);
          agntVel = cat(3,agntVel,agntVel(:,:,end));
          bioPos = cat(3,bioPos,single(agntPos(:,:,1:pad)));
          bioVel = cat(3,bioVel,single(agntVel(:,:,1:pad)));
          clear angtPos agntVel
          bioPos = bioPos(bioPos(:,3,end)>0,:,:);
          bioVel = bioVel(bioPos(:,3,end)>0,:,:);
          nAgent = size(bioPos,1);
        end
    else
        bioPos = single(agntPos(:,:,mod(floor((volCtr-1)*volTime)+1,size(agntPos,3)):mod(ceil(volCtr*volTime)+1,size(agntPos,3))));
        bioPos = bioPos(bioPos(:,3,end)>-5,:,:);
        bioVel = single(agntVel(:,:,mod(floor((volCtr-1)*volTime)+1,size(agntPos,3)):mod(ceil(volCtr*volTime)+1,size(agntPos,3))));
        bioVel = bioVel(bioPos(:,3,end)>-5,:,:);
        nAgent = size(bioPos,1);
        clear agntPos agntVel
    end
    if interactiveIO, disp('Bio model input successfully loaded.'), end
    

    
    %% number of PRTs per biologial model time step
    PRTperMod = ceil(dt/PRT);


    %% shift bio coordinate system to radar-relative cartesian coordinates
    bioPos(:,1,:) = bioPos(:,1,:) + bOrX;
    bioPos(:,2,:) = bioPos(:,2,:) + bOrY;
    bioPos(:,3,:) = bioPos(:,3,:) + bOrZ;
    
    
    %% initialize counters
    pulseCtr     = 1;
    
    %% loop over each bio position at resolution of biological model (dt)
    for bioModTmCtr = 1:size(bioPos,3)
        
        disp(bioModTmCtr)
         
        %% If this is the last iteration for this volume scan, make sure
        %%data ends at end of scan
        if (pulseCtr+PRTperMod) > length(azTm)
          PRTperMod = length(azTm)-pulseCtr;
          if PRTperMod == 0
              %% arrange pulses into sweeps
              if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
                tmp = reshape(Vh,dwellNum,nAzs,nEls,nGates);
                Vh  = permute(tmp,[2 1 4 3]);
              end
              if strcmp(polType,'STSR')||strcmp(polType,'VERT')
                tmp = reshape(Vv,dwellNum,nAzs,nEls,nGates);
                Vv = permute(tmp,[2 1 4 3]);
              end

              %% save time series data for completed volume scan
              saveName = [iqSaveName num2str(ceil(volCtr*volTime)) '.mat'];
              if interactiveIO,fprintf('Saving: %s/\n',saveName), end
              if strcmp(polType,'STSR')
                save('-v7.3',saveName,'Vh','Vv','rngSet','elSet','azSet','volTime','header')
              elseif strcmp(polType,'HORZ')
                save('-v7.3',saveName,'Vh','rngSet','elSet','azSet','volTime','header')
              elseif strcmp(polType,'VERT')
                save('-v7.3',saveName,'Vv','rngSet','elSet','azSet','volTime','header')
              end
              if interactiveIO,fprintf('Volume scan %g of %g completed.\n',volCtr,nVols), end
             clear bioPos bioVel saveName Vv Vh
              break
          end
        end
        
        
        %% convert agents to radar-relative spherical coordinates (az,el,rng)
        bioPosSph(:,1) = 180/pi*atan2(bioPos(:,1,bioModTmCtr),bioPos(:,2,bioModTmCtr));
        bioPosSph(:,3) = sqrt(bioPos(:,1,bioModTmCtr).^2+bioPos(:,2,bioModTmCtr).^2+bioPos(:,3,bioModTmCtr).^2);
        bioPosSph(:,2) = asind(bioPos(:,3,bioModTmCtr)./bioPosSph(:,3));
        bioPosSph(bioPosSph(:,1)<0,1) = bioPosSph(bioPosSph(:,1)<0,1)+360;
  
         
        %% account for standard refraction using 4/3 Earth model
        refracH = (sqrt(bioPosSph(:,3).^2+(k_e*r_earth)^2+...
          2.*sind(bioPosSph(:,2))*k_e*r_earth.*bioPosSph(:,3))-k_e*r_earth);
        refracD = k_e*r_earth.*asin(bioPosSph(:,3).*cosd(bioPosSph(:,2))./(k_e*r_earth+refracH));
        
        %% modify agent locations to account for refraction offset
        bioPosSph(:,3) = sqrt(refracH.^2+refracD.^2);
        bioPosSph(:,2) = atand(refracH./refracD);
        clear refracH refracD
  
        
        %% check bounds on radar resolution volumes during this time step
        azBds = [azTm(pulseCtr)-4*bmWid azTm(pulseCtr+PRTperMod)+4*bmWid];   
        elBds = [elTm(pulseCtr)-4*bmWid elTm(pulseCtr+PRTperMod)+4*bmWid];
        
        
        %% find agents within resolution volumes 
        % (Find indexes of scatterers in the volume covered in this time
        % step
        if azBds(1)<azBds(2)
          bioId = find(bioPosSph(:,1)>=azBds(1) & bioPosSph(:,1)<=azBds(2)...
              & bioPosSph(:,2)>=elBds(1) & bioPosSph(:,2)<=elBds(2)...
              & bioPosSph(:,3)<=(rngSet(end)) & bioPosSph(:,3)>500 ...
              & bioPos(:,3,bioModTmCtr)>-1);
        else
          bioId = find(((bioPosSph(:,1)>=azBds(1)&bioPosSph(:,1)<=360)|...
              (bioPosSph(:,1)<=azBds(2)&bioPosSph(:,1)>=0))...
              & bioPosSph(:,2)>=elBds(1) & bioPosSph(:,2)<=elBds(2)...
              & bioPosSph(:,3)<=(rngSet(end))& bioPosSph(:,3)>500 ...
              & bioPos(:,3,bioModTmCtr)>-1);
        end
        clear azBds elBds
        
        
        %% if the length of bioId is 1, we have problems with squeezing out
        %% a dimension. So we'll just add another agent.
        if length(bioId) <= 1
          bioId = vertcat(bioId,[1 2]');
        end
        bioId = 1:nAgent;disp('L-322')
        
        %% calculate radial velocity for agents during this time step
        radVel(:,1) = squeeze(bioVel(bioId,1,bioModTmCtr).*cosd(bioPosSph(bioId,1)-90).*cosd(bioPosSph(bioId,2)));
        radVel(:,2) = squeeze(bioVel(bioId,2,bioModTmCtr).*cosd(bioPosSph(bioId,1)).*cosd(bioPosSph(bioId,2)));
        radVel(:,3) = squeeze(bioVel(bioId,3,bioModTmCtr).*cosd(90-bioPosSph(bioId,2)));
        RadVelSign = -1*(acosd(dot(bioPos(bioId,:,bioModTmCtr),radVel,2)./(sqrt(sum(bioPos(bioId,:,bioModTmCtr).^2,2)).*sqrt(sum(radVel.^2,2))))>90);
        RadVelSign(RadVelSign==0) = 1;
        radVel      = RadVelSign .* sqrt(...
          (bioVel(bioId,1,bioModTmCtr).*cosd(bioPosSph(bioId,1)-90).*cosd(bioPosSph(bioId,2))).^2 ...
          +(bioVel(bioId,2,bioModTmCtr).*cosd(bioPosSph(bioId,1)).*cosd(bioPosSph(bioId,2))).^2 ...
          +(bioVel(bioId,3,bioModTmCtr).*cosd(90-bioPosSph(bioId,2))).^2 );
        clear RadVelSign
        
      
        %% interpolate positions to PRT
        initialPos  = repmat(bioPos(bioId,:,bioModTmCtr),[1,1,PRTperMod]);
        initialVel  = repmat(bioVel(bioId,:,bioModTmCtr),[1,1,PRTperMod]);
        temp        = repmat(0:PRT:dt-PRT,[length(bioId),1])';
        temp        = temp(1:PRTperMod,:);  
        timeMatrix  = cat(3,temp,temp,temp);
        clear temp
        timeMatrix  = shiftdim(timeMatrix,1);
        bioPosIntrp = initialPos+initialVel.*timeMatrix;
        clear initialPos initialVel timeMatrix
        
        
        %% use velocity vectors to find orientation wrt radar 
        a_orien = sign(bioVel(bioId,1,bioModTmCtr).*cosd(bioPosSph(bioId,1))+...
          bioVel(bioId,2,bioModTmCtr).*sind(bioPosSph(bioId,1))).*...
          acosd((bioPos(bioId,1,bioModTmCtr).*bioVel(bioId,1,bioModTmCtr)+...
          bioPos(bioId,2,bioModTmCtr).*bioVel(bioId,2,bioModTmCtr))./...
          (sqrt(bioPos(bioId,1,bioModTmCtr).^2+bioPos(bioId,2,bioModTmCtr).^2).*...
          sqrt(bioVel(bioId,1,bioModTmCtr).^2+bioVel(bioId,2,bioModTmCtr).^2)));
        b_orien = -bioPosSph(bioId,2);
        clear bioPosSph
        
        
        
        %% use orientation to find RCS
        
        for ct = 1:length(a_orien)
        aId = find(azBat==round(a_orien(ct)),1,'first');
        bId = find(elBat==round(b_orien(ct)),1,'first');
        
        RCSH(ct,1) = rcsHH(bId,aId);
        RCSV(ct,1) = rcsVV(bId,aId);
        RCSHV(ct,1)= rcsHV(bId,aId);
        phsH(ct,1) = phHH(bId,aId);
        phsV(ct,1) = phVV(bId,aId);
        phsHV(ct,1)= phHV(bId,aId);
        end

        
        %% make boresight unit vector for each time step
        [cartBoreX cartBoreY cartBoreZ] = sph2cart(...
          repmat(-pi*(azTm(pulseCtr:pulseCtr+PRTperMod-1)-90)./180,[length(bioId) 1]),...
          repmat(pi*elTm(pulseCtr:pulseCtr+PRTperMod-1)./180,[length(bioId) 1]),...
          ones(length(bioId),PRTperMod));
        
        
        %% find angle of each agent with respect to boresight for each time step
        theta = real(acosd((cartBoreX.*squeeze(bioPosIntrp(:,1,:))+...
          cartBoreY.*squeeze(bioPosIntrp(:,2,:))+...
          cartBoreZ.*squeeze(bioPosIntrp(:,3,:)))./...
          (sqrt(cartBoreX.^2+cartBoreY.^2+cartBoreZ.^2).*...
          sqrt(squeeze(bioPosIntrp(:,1,:)).^2+...
          squeeze(bioPosIntrp(:,2,:)).^2+squeeze(bioPosIntrp(:,3,:)).^2))));
        clear cartBoreX cartBoreY cartBoreZ bioId
        
        
        %% find range to each agent for each time step
        rng = squeeze(sqrt(bioPosIntrp(:,1,:).^2+bioPosIntrp(:,2,:).^2+bioPosIntrp(:,3,:).^2));
        clear bioPosIntrp
        
        
        %% find beam weight for each agent for each time step (via D-Z p.34)
        dishDim   = 1.27*lam/(pi*bmWid/180);
        term1     = pi*dishDim*sind(theta)/lam;
        term2     = besseli(2,term1);
        bmWt      = -10*log10(((8*term2./term1.^2).^2))+gain;
        bmWt      = sqrt(10.^((bmWt)/10)./10.^(gain/10));
        clear dishDim term1 term2 theta
        
        
        %% find range gate for each agent for each time step
        gateNum = ceil((rng-rngSet(1))/rWid6);
        
        
        %% handle errors when agent location is NaN
        gateNum(isnan(gateNum)) = 1;
        
        
        %% handle errors when agents are beyond the radar range
        gateNum(gateNum>length(rngSet)) = 1;
        rng(gateNum>length(rngSet))     = nan;
        
        
        %% handle errors when agents are within the minimum radar range
        gateNum(gateNum==0) = 1;
        rng(gateNum==0)     = nan;
        
        
        %% find distance from center of range gate
        rngDist = rng-rngSet(gateNum);
        
        
        %% find range weight for each agent for each time step (via D-Z p.79)
        a       = pi/(2*sqrt(log(2)));
        b       = rxBndWid*pi/(4*sqrt(log(2)));
        x       = 2*a*rxBndWid/cc*rngDist;
        rngWt   = (erf(x+b)-erf(x-b))/2;
        clear a b x rngDist

        RCSH  = repmat(RCSH,[1 PRTperMod]);
        RCSV  = repmat(RCSV,[1 PRTperMod]);
        RCSHV  = repmat(RCSHV,[1 PRTperMod]);
        %% calculate echo power from each agent (D-Z p. 46, eq. 3.24)
        if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
          rxPowH = txPow/2*(10^(gain/10)).^2*lam^2*RCSH.*rngWt.*bmWt.^4./((4*pi).^3*rng.^4);
          rxPowHV= txPow/2*(10^(gain/10)).^2*lam^2*RCSHV.*rngWt.*bmWt.^4./((4*pi).^3*rng.^4);
        end
        if strcmp(polType,'STSR')||strcmp(polType,'VERT')
          rxPowV = txPow/2*(10^(gain/10)).^2*lam^2*RCSV.*rngWt.*bmWt.^4./((4*pi).^3*rng.^4);
          rxPowVH= txPow/2*(10^(gain/10)).^2*lam^2*RCSHV.*rngWt.*bmWt.^4./((4*pi).^3*rng.^4);
        end
        clear bmWt rngWt RCSH RCSV RCSHV
         
        
        %% calculate echo phase from each agent
        radVel  = repmat(radVel,[1 PRTperMod]);
        phsH  = repmat(phsH,[1 PRTperMod]);
        phsV  = repmat(phsV,[1 PRTperMod]);
        phsHV  = repmat(phsHV,[1 PRTperMod]);
        
        
        phaseH = 4/lam*pi*rng+4/lam*pi*PRT*radVel+phsH;
        phaseV = 4/lam*pi*rng+4/lam*pi*PRT*radVel+phsV+sysDP;
        phaseHV = 4/lam*pi*rng+4/lam*pi*PRT*radVel+phsHV;
        phaseVH = 4/lam*pi*rng+4/lam*pi*PRT*radVel+phsHV+sysDP;
        clear radVel rng phsH phsV phsHV %timeSet
        
        
        %% calculate complex echo voltage for each agent and pulse
        if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
          IhAgnt = ((rxPowH/2).^(.5).*cos(phaseH))+((rxPowVH/2).^(.5).*cos(phaseVH));
          QhAgnt = (-(rxPowH/2).^(.5).*sin(phaseH))+(-(rxPowVH/2).^(.5).*sin(phaseVH));
          VhAgnt = IhAgnt - ii.*QhAgnt;
        end
        if strcmp(polType,'STSR')||strcmp(polType,'VERT')
          IvAgnt = ((rxPowV/2).^(.5).*cos(phaseV))+((rxPowHV/2).^(.5).*cos(phaseHV));
          QvAgnt = (-(rxPowV/2).^(.5).*sin(phaseV))+(-(rxPowHV/2).^(.5).*sin(phaseHV));
          VvAgnt = IvAgnt - ii.*QvAgnt;
        end
        clear IhAgnt QhAgnt IvAgnt QvAgnt rxPowH rxPowV rxPowHV rxPowVH...
          phaseH phaseV phaseHV phaseVH
        
        
        %% sum echoes in each range gate to get total signal
        if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
          sumH = zeros(nGates,PRTperMod);
        end
        if strcmp(polType,'STSR')||strcmp(polType,'VERT')
          sumV = zeros(nGates,PRTperMod);
        end
        for gateCtr = 1:nGates
          gateMask = gateNum==gateCtr;
          if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
              tempH  = nansum(VhAgnt.*gateMask);
              sumH(gateCtr,:) = tempH;
          end
          if strcmp(polType,'STSR')||strcmp(polType,'VERT')
              tempV  = nansum(VvAgnt.*gateMask);
              sumV(gateCtr,:) = tempV;
          end
        end 
        clear gateCtr gateMask gateNum tempH tempV        
        
       
        %% concatenate pulses
        if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
          Vh(pulseCtr:pulseCtr+PRTperMod-1,:) = sumH';
        end
        if strcmp(polType,'STSR')||strcmp(polType,'VERT')
          Vv(pulseCtr:pulseCtr+PRTperMod-1,:) = sumV';
        end
        clear sumH sumV
        
        
        %% iterate pulse counter
        pulseCtr = pulseCtr + PRTperMod;
        
        
    end % end bio model time step
    
    
end % end radar volume scans

%% Reshape stuff for plotting
if strcmp(polType,'STSR')||strcmp(polType,'HORZ')
    tmp = reshape(Vh,dwellNum,nAzs,nEls,nGates);
    Vh  = permute(tmp,[2 1 4 3]);
end
if strcmp(polType,'STSR')||strcmp(polType,'VERT')
    tmp = reshape(Vv,dwellNum,nAzs,nEls,nGates);
    Vv = permute(tmp,[2 1 4 3]);
end

 %% plot stuff


elId = 1;
azSet_tmp = [azSet azSet(1)];
[AZ DD] = meshgrid(azSet_tmp,rngSet);
EL = elSet(elId)*ones(size(AZ));

yy = DD.*cosd(AZ).*cosd(EL);
xx = DD.*sind(AZ).*cosd(EL);
powH = squeeze(nanmean(abs(double(Vh(:,:,:,elId))),2))';
powH = [powH powH(:,1)];
powH(10*log10(powH)<-120)=nan;

powV = squeeze(nanmean(abs(double(Vv(:,:,:,elId))),2))';
size(powV)
powV = [powV powV(:,1)];
powV(10*log10(powV)<-120)=nan;

figure
pcolor(xx, yy, 10*log10(powH)),colorbar
shading flat
axis square
hold on
spurRng=[0 140000];
spurAz = [0:45:360];
[spurAz, spurRng] = meshgrid(spurAz,spurRng);
spurXX = spurRng.*sind(spurAz);
spurYY = spurRng.*cosd(spurAz);
ringAz = 0:.25:360;
ringRng = [25000:25000:200000];
[ringRng, ringAz] = meshgrid(ringRng,ringAz);
ringXX = ringRng.*sind(ringAz);
ringYY = ringRng.*cosd(ringAz);
plot(ringXX,ringYY,'k','linewidth',1)
plot(spurXX,spurYY,'k','linewidth',1)

plot(bioPos(:,1,21),bioPos(:,2,21),'k.')
plot(0,0,'ko')
axis([-50000 50000 -50000 50000])
