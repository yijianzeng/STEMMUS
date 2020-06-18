function h_BC
global RHS KL_h Precip Evap NN C4 Trap
global NBCh NBChB BCh BChB hN KT Delt_t DSTOR0 NBChh TIME h_SUR AVAIL0 C4_a
%%%%%%%%%% Apply the bottom boundary condition called for by NBChB %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBChB==1            %-----> Specify matric head at bottom to be ---BChB;
    RHS(1)=BChB;  
    C4(1,1)=1;
    RHS(2)=RHS(2)-C4(1,2)*RHS(1);
    C4(1,2)=0; 
    C4_a(1)=0;
elseif NBChB==2        %-----> Specify flux at bottom to be ---BChB (Positive upwards);
    RHS(1)=RHS(1)+BChB;
elseif NBChB==3        %-----> NBChB=3,Gravity drainage at bottom--specify flux= hydraulic conductivity;
    RHS(1)=RHS(1)-KL_h(1,1);   
end

%%%%%%%%%% Apply the surface boundary condition called for by NBCh  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBCh==1             %-----> Specified matric head at surface---equal to hN;
    RHS(NN)=h_SUR(KT); 
    C4(NN,1)=1;
    RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
    C4(NN-1,2)=0; 
    C4_a(NN-1)=0;
elseif NBCh==2
    if NBChh==1
        RHS(NN)=hN; 
        C4(NN,1)=1;
        RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
        C4(NN-1,2)=0; 
    else
        RHS(NN)=RHS(NN)-BCh;   %> and a specified matric head (saturation or dryness)was applied; 
    end    
else     
    Evap_Cal;
    if TIME>=80*3600 && TIME<=82*3600     %6-26 8-10 p=13.1mm
        Precip(KT)=0.25/3600;
    elseif TIME>=223*3600 && TIME<=227*3600   %7-2 7-11 p=9mm
        Precip(KT)=0.18/3600;
     elseif TIME>=489*3600 && TIME<=490*3600    %7-13 9-14 p=52mm
        Precip(KT)=5/2/3600;
        % elseif TIME>=583*3600 && TIME<=586*3600    %7-17 7-10 p=2mm
        %     Precip(KT)=0.05/3600;
    elseif TIME>=703*3600 && TIME<=706*3600   %7-22 7-8 p=5.5mm
        Precip(KT)=0.05/3600;
    elseif TIME>=952*3600 && TIME<=954*3600  %8-1 16-18 p=60mm
        Precip(KT)=6/3/3600;
        %         elseif TIME>80*3600 && TIME<82*3600  %8-11 7-8 p=0.8mm
        %         Precip(KT)=0.66/3600;
    elseif TIME>=1288*3600 && TIME<=1290*3600  %8-15 16-17 p=67mm
        Precip(KT)=6.7/3/3600;
    elseif TIME>=1862*3600 && TIME<=1866*3600  %9-8 14-18 p=93.11mm
        Precip(KT)=5.3/5/3600;
    else
        Precip(KT)=0;
    end
    
    AVAIL0=Precip(KT)+DSTOR0/Delt_t;
    if NBChh==1
        RHS(NN)=hN; 
        C4(NN,1)=1;
        RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
        C4(NN-1,2)=0; 
        C4_a(NN-1)=0;
    else
        RHS(NN)=RHS(NN)+AVAIL0-Evap(KT); 
    end
end

