clc;
clearvars -except TRAP_ind EVAP_ind Sim_Temp_ind Sim_Theta_ind;
close all;
run Constants
if Evaptranp_Cal==1
    run EvapTransp_Cal;
end
% function MainLoop
global i tS KT Delt_t TEND TIME MN NN NL ML ND hOLD TOLD h hh T TT P_gOLD P_g P_gg Delt_t0
global KIT NIT TimeStep Processing
global SUMTIME hhh TTT P_ggg Theta_LLL DSTOR Thmrlefc CHK Theta_LL Theta_L
global NBCh AVAIL Evap DSTOR0 EXCESS QMT RS BCh hN hSAVE NBChh DSTMAX Soilairefc Trap sumTRAP_dir sumEVAP_dir
global TSAVE IRPT1 IRPT2 AVAIL0 TIMEOLD TIMELAST SRT ALPHA BX alpha_h bx Srt
global QL QL_h QL_T QV Qa KL_h Chh ChT Khh KhT TRAP_ind


run StartInit;   % Initialize Temperature, Matric potential and soil air pressure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TIMEOLD=0;
TIMELAST=0;
for i=1:tS+1                          % Notice here: In this code, the 'i' is strictly used for timestep loop and the arraies index of meteorological forcing data.
    KT=KT+1                          % Counting Number of timesteps
    if KT>1 && Delt_t>(TEND-TIME)
        Delt_t=TEND-TIME;           % If Delt_t is changed due to excessive change of state variables, the judgement of the last time step is excuted.
    end
    TIME=TIME+Delt_t;               % The time elapsed since start of simulation
    TimeStep(KT,1)=Delt_t;
    SUMTIME(KT,1)=TIME;
    Processing=TIME/TEND
    %%%%%% Updating the state variables. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if TIME>=489*3600 && TIME<=492*3600    %7-13 9-10 p=52mm
        NBChh=1;
    elseif TIME>=952*3600 && TIME<=964*3600  %8-1 16-18 p=60mm
        NBChh=1;
    elseif TIME>=1288*3600 && TIME<=1297*3600  %8-15 16-17 p=67mm
        NBChh=1;
    elseif TIME>=1862*3600 && TIME<=1870*3600  %9-8 14-18 p=93.11mm
        NBChh=1;
    else
        NBChh=2;
    end
    if IRPT1==0 && IRPT2==0
        for MN=1:NN
            hOLD(MN)=h(MN);
            h(MN)=hh(MN);
            hhh(MN,KT)=hh(MN);
            KL_h(MN,KT)=KL_h(MN,2);
            Chh(MN,KT)=Chh(MN,2);
            ChT(MN,KT)=ChT(MN,2);
            Khh(MN,KT)=Khh(MN,2);
            KhT(MN,KT)=KhT(MN,2);
            
            if Thmrlefc==1
                TOLD(MN)=T(MN);
                T(MN)=TT(MN);
                TTT(MN,KT)=TT(MN);
            end
            if Soilairefc==1
                P_gOLD(MN)=P_g(MN);
                P_g(MN)=P_gg(MN);
                P_ggg(MN,KT)=P_gg(MN);
            end
            if rwuef==1
                SRT(MN,KT)=Srt(MN,1);
                ALPHA(MN,KT)=alpha_h(MN,1);
                BX(MN,KT)=bx(MN,1);
            end
        end
        DSTOR0=DSTOR;
        if KT>1
            run SOIL1
        end
    end
    
    if Delt_t~=Delt_t0
        for MN=1:NN
            hh(MN)=h(MN)+(h(MN)-hOLD(MN))*Delt_t/Delt_t0;
            TT(MN)=T(MN)+(T(MN)-TOLD(MN))*Delt_t/Delt_t0;
        end
    end
    hSAVE=hh(NN);
    TSAVE=TT(NN);
    if NBCh==1
        hN=BCh;
        hh(NN)=hN;
        hSAVE=hN;
    elseif NBCh==2
        if NBChh~=2
            if BCh<0
                hN=DSTOR0;
                hh(NN)=hN;
                hSAVE=hN;
            else
                hN=-1e6;
                hh(NN)=hN;
                hSAVE=hN;
            end
        end
    else
        if NBChh~=2
            hN=DSTOR0;
            hh(NN)=hN;
            hSAVE=hN;
        end
    end
    run Forcing_PARM
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for KIT=1:NIT   % Start the iteration procedure in a time step.
        
        run SOIL2;
        run CondL_T;
        run Density_V;
        run CondL_Tdisp;
        run Latent;
        run Density_DA;
        run CondT_coeff;
        run Condg_k_g;
        run CondV_DE;
        run CondV_DVg;
        
        run h_sub;
        
        if NBCh==1
            DSTOR=0;
            RS=0;
        elseif NBCh==2
            AVAIL=-BCh;
            EXCESS=(AVAIL+QMT(KT))*Delt_t;
            if abs(EXCESS/Delt_t)<=1e-10,EXCESS=0;end
            DSTOR=min(EXCESS,DSTMAX);
            RS=(EXCESS-DSTOR)/Delt_t;
        else
            AVAIL=AVAIL0-Evap(KT);
            EXCESS=(AVAIL+QMT(KT))*Delt_t;
            if abs(EXCESS/Delt_t)<=1e-10,EXCESS=0;end
            DSTOR=0;
            RS=0;
        end
        
        if Soilairefc==1
            run Air_sub;
        end
        
        if Thmrlefc==1
            run Enrgy_sub;
        end
        
        if max(CHK)<0.001
            break
        end
        hSAVE=hh(NN);
        TSAVE=TT(NN);
    end
    TIMEOLD=KT;
    KIT
    KIT=0;
    run SOIL2;
    % run TimestepCHK
    
    if IRPT1==0 && IRPT2==0
        if KT        % In case last time step is not convergent and needs to be repeated.
            MN=0;
            for ML=1:NL
                for ND=1:2
                    MN=ML+ND-1;
                    Theta_LLL(ML,ND,KT)=Theta_LL(ML,ND);
                    Theta_L(ML,ND)=Theta_LL(ML,ND);
                    
                end
            end
            run ObservationPoints
        end
        if (TEND-TIME)<1E-3
            for MN=1:NN
                hOLD(MN)=h(MN);
                h(MN)=hh(MN);
                hhh(MN,KT)=hh(MN);
                if Thmrlefc==1
                    TOLD(MN)=T(MN);
                    T(MN)=TT(MN);
                    TTT(MN,KT)=TT(MN);
                end
                if Soilairefc==1
                    P_gOLD(MN)=P_g(MN);
                    P_g(MN)=P_gg(MN);
                    P_ggg(MN,KT)=P_gg(MN);
                end
            end
            break
        end
    end
    for MN=1:NN
        QL(MN,KT)=QL(MN);
        QL_h(MN,KT)=QL_h(MN);
        QL_T(MN,KT)=QL_T(MN);
        Qa(MN,KT)=Qa(MN);
        QV(MN,KT)=QV(MN);
    end
end
% run PlotResults
if Evaptranp_Cal==1   % save the variables for ETind scenario
    Sim_Theta_ind=Sim_Theta;
    Sim_Temp_ind=Sim_Temp;
    TRAP=36000.*Trap;
    TRAP_ind=TRAP';
    EVAP=36000.*Evap;
    EVAP_ind=EVAP';
    disp ('Convergence Achieved for ETind scenario. Please switch to ETdir scenario and run again.')
else
    TRAP=36000.*Trap;
    TRAP_dir=TRAP';
    EVAP=36000.*Evap;
    EVAP_dir=EVAP';
    for i=1:KT/24
%         sumTRAP_ind(i)=0;
%         sumEVAP_ind(i)=0;
        sumTRAP_dir(i)=0;
        sumEVAP_dir(i)=0;
        for j=(i-1)*24+1:i*24
%             sumTRAP_ind(i)=TRAP_ind(j)+sumTRAP_ind(i);
%             sumEVAP_ind(i)=EVAP_ind(j)+sumEVAP_ind(i);
            sumTRAP_dir(i)=TRAP(j)+sumTRAP_dir(i);
            sumEVAP_dir(i)=EVAP(j)+sumEVAP_dir(i);
        end
    end
    run PlotResults1
end
