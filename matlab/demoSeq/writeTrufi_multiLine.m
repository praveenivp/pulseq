% we define here a really crude True-FISP a.k.a. bSSFP sequence
% there is no user control for TR/TE, you just specify the ADC time and the
% RF parameters and the rest is calculated to find the fastest posible
% timing. The sequence intensively uses the extended trapezoid
% functionality to achieve near-optimal timing. Due to the requirement for
% splitting the sequence into blocks the TR is increased by approximately
% 40-60 us (rfDeadTime+adcDeadTime) in comparison to the true minimum TR

% set system limits
% had to slow down ramps and increase adc_duration to avoid stimulation
sys = mr.opts('MaxGrad',50,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 00e-6, 'rfDeadTime', 0e-6, ...
    'adcDeadTime', 20e-6,'B0',9.4);

seq=mr.Sequence(sys);              % Create a new sequence object
fov=220e-3; Nx=floor(fov/0.8e-3); Ny=floor(fov/0.8e-3);     % Define FOV and resolution


Ry=4;
TR_desired=10e-3;


% ADC duration (controls TR/TE)
dwell=3.5; %us
adc_dur=dwell*Nx; %us

% RF parameters
alpha=15; % deg
thick=0.8*16; %mm
rf_dur=1000; % us
rf_apo=0.5;
rf_bwt=4;
 slrFac=1;
% Create 'alpha' degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,'system',sys,'maxslew',slrFac*sys.maxSlew);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys,'maxslew',slrFac*sys.maxSlew);
gx_Neg = mr.makeTrapezoid('x','FlatArea',-Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys,'maxslew',slrFac*sys.maxSlew);
PreX=mr.makeTrapezoid('x','Area',0.5*gx.area,'maxslew',slrFac*sys.maxSlew,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
% gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;


prephaser_max =mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'system',sys);
prephaser_dur_max=mr.calcDuration(Prephaser_max);


gy_blips=mr.makeTrapezoid('y','Area',deltak*Ry,'delay',gx.riseTime+gx.flatTime,'system',sys,'maxslew',slrFac*sys.maxSlew);
ADC_time=mr.calcDuration(gx);
ReadOut=TR_desired- (mr.calcDuration(gz)+mr.calcDuration(PreX)*2);

Nlines=floor(ReadOut/ADC_time);
Nseg=ceil(Ny/(Nlines*Ry));

% Loop over phase encodes and define sequence blocks
%X is read
Lines_list=1:Ry:Ny;
Lines_list=Lines_list(1:Nlines:end);
phaseAreas = ((0:Ry:Ny-1)-Ny/2)*deltak;
phaseAreas_pre=phaseAreas(1:Nlines:end);
phaseAreas_post=phaseAreas_pre+Nlines*deltak*Ry;
for i=1:(Nseg)
    rf.phaseOffset=pi*mod(i,2);
    adc.phaseOffset=pi*mod(i,2);
    
    seq.addBlock(rf,gz);
    %1:Rxy:Ny lines needs to be measured
   
    %prePhasers
    
    PreY=mr.makeTrapezoid('y','Area',phaseAreas_pre(i),'maxslew',slrFac*sys.maxSlew,'system',sys);
    %Rephasers
    if(mod(Nlines,2)==0)
        PostX=mr.makeTrapezoid('x','Area',-0.5*gx.area,'maxslew',slrFac*sys.maxSlew,'system',sys);
    else
        PostX=mr.makeTrapezoid('x','Area',0.5*gx.area,'maxslew',slrFac*sys.maxSlew,'system',sys);
    end
    PostY=mr.makeTrapezoid('y','Area',-phaseAreas_post(i),'maxslew',slrFac*sys.maxSlew,'system',sys);
    seq.addBlock(PreX,PreY,gzReph);
    for ml=1:Nlines
        if(mod(ml,2)==0)
            seq.addBlock(gx,adc,gy_blips);
        else
            seq.addBlock(gx_Neg,adc,gy_blips);
        end
        
    end
    seq.addBlock(PostX,PostY,gzReph);
    
    if(length(seq.blockDurations)<15)
        TR_actual=sum(seq.blockDurations);
    end
    
     seq.addBlock(mr.makeDelay(TR_desired-TR_actual))
    
end
fprintf('Sequence ready\n');

seq.plot('timeDisp','ms')
 xlim([0 TR_desired*1e3])
% fprintf('TR=%03f ms  TE=%03f ms\n', TR*1e3, TE*1e3);

fprintf('Echo Spacing: %2.3f ms \n',1e3*(gx.flatTime+2*gx.riseTime))
fprintf('read BW/pxl: %2.3f ms \n',1/(gx.flatTime))
fprintf('%d shots and %d lines per short\n',Nseg,Nlines)
fprintf('toatalreadour %.2f ms \n',1e3*Nlines*gx.flatTime)


%% get trajectory in
[gw]=seq.waveforms_and_times();
% tq=0:10e-6:max([gw{1}(1,end),gw{2}(1,end),gw{3}(1,end)]);
tq=0:10e-6:(Nseg)*TR_desired-10e-6;
G_RPS=zeros(3,length(tq));
for i=1:3
G_RPS(i,:)=interp1(gw{i}(1,:),gw{i}(2,:),tq,'linear');
end
G_RPS(isnan(G_RPS))=0;
G_RPS=reshape(G_RPS,3,[],Nseg);
% figure,plot(G_RPS(:,:)')

k_epi=cumsum(G_RPS(:,:)./42.567e3,2)*sp.Gyro*10e-6;
k_epi=reshape(k_epi,size(G_RPS));

% t_epi_adc=abs(abs(G_RPS(1,:,1))-gx.amplitude)<1;
t_epi_adc=abs(G_RPS(1,:,1)-gx.amplitude)<1 | circshift(abs(G_RPS(1,:,1)-gx_Neg.amplitude)<1,-1,2);
k_epi(:,t_epi_adc==0,1:Nseg)=nan;
figure(19),clf,plotk(k_epi(:,:,:))

%% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces

 figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
    ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
title('2D k-space');


fprintf('Echo Spacing: %2.3f ms \n',1e3*(gx.flatTime+2*gx.riseTime))
fprintf('read BW/pxl: %2.3f ms \n',1/(gx.flatTime))
fprintf('%d shots and %d lines per short\n',Nseg,Nlines)

filename=sprintf('bSSFP2D_ML%d_TR%d_R%d_DW%d_FA%d.seq',Nlines,round(TR_desired*1e3),Ry,dwell,alpha);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare export

seq.setDefinition('FOV', [fov fov thick*1e-3]);
seq.setDefinition('Name', 'trufi');
seq.write(fullfile('\\mrz10\MRZ9T\Upload9T\USERS\Praveen\VE12U_builds\pulseq',filename))       % Write to pulseq file

% seq.install('siemens');

%% plot entire interpolated wave forms -- good for debugging of gaps, jumps,
% etc, but is relatively slow
%gw=seq.gradient_waveforms();
%figure; plot(gw'); % plot the entire gradient waveform
gw_data=seq.waveforms_and_times();
figure; plot(gw_data{1}(1,:),gw_data{1}(2,:),gw_data{2}(1,:),gw_data{2}(2,:),gw_data{3}(1,:),gw_data{3}(2,:)); % plot the entire gradient waveform
title('gradient waveforms');


%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits

rep = seq.testReport;
fprintf([rep{:}]);

