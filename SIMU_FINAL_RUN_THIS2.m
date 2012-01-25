
%
%
%
%   RADAR SIMULATION
%
%   *** RUN THIS FILE ***
%
%
%


clc;
clear all;
close all;

f=3*10^9;           % frequency
c=3*10^8;           % velocity
lambda=c/f;         % wavelength

k=2*pi/lambda;

len_dipole=lambda/10; % dipole length

rot_l=-45;          % rotation angle lower limit
rot_h=45;           % rotation angle upper limit
rot_step=1;         % rotation angle step 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plate-1


p1_l=5;                     % lengthmax(E_plate1_mag)
p1_w=3;                     % width
p1_angle=30;                % angle
p1_dist=5-2.5*sin(pi*30/180);% plate-1 distance from reference plate


% Calculating RCS for plate-1

E_plate1=RCS_PLATE(p1_l,p1_angle,p1_w,rot_l,rot_h,rot_step,len_dipole,k,p1_dist);

len = size(E_plate1);
    
for n=1:len(2)
    E_plate1_mag(n)=abs(E_plate1(n));   % Magnitude
   
end


maxx=max(E_plate1_mag);
x= rot_l:rot_step:rot_h;    % angles




% Plate-2


p2_l=3;
p2_w=3;
p2_angle=0;
p2_dist=0;

% Calculating RCS for plate-1

E_plate2=RCS_PLATE(p2_l,p2_angle,p2_w,rot_l,rot_h,rot_step,len_dipole,k,p2_dist);

len = size(E_plate2);
    
for n=1:len(2)
    E_plate2_mag(n)=abs(E_plate2(n));   % Magnitude
end







% Plate-3

p3_l=5;
p3_w=3;
p3_angle=-30;
p3_dist=5-2.5*sin(pi*30/180);


% Calculating RCS for plate-1

E_plate3=RCS_PLATE(p3_l,p3_angle,p3_w,rot_l,rot_h,rot_step,len_dipole,k, p3_dist);

len = size(E_plate3);
    
for n=1:len(2)
    E_plate3_mag(n)=abs(E_plate3(n));   % Magnitude
end




E_plates= E_plate1_mag + E_plate2_mag + E_plate3_mag;




% For Scatterer

ref_dist=5;

rcs_scat1=1;
ang_scat1=30;
dist_scat1=3;

e_scat1 = RCS_SCATT(rcs_scat1,ref_dist, dist_scat1, ang_scat1, rot_l, rot_h, rot_step,k,len_dipole);


rcs_scat2=1;
ang_scat2=0;
dist_scat2=3;

e_scat2 = RCS_SCATT(rcs_scat2,ref_dist, dist_scat2, ang_scat2, rot_l, rot_h, rot_step,k,len_dipole);



rcs_scat3=1;
ang_scat3=-30;
dist_scat3=3;

e_scat3 = RCS_SCATT(rcs_scat3,ref_dist, dist_scat3, ang_scat3, rot_l, rot_h, rot_step,k,len_dipole);



rcs_scat4=1;
ang_scat4=-130;
dist_scat4=3;

e_scat4 = RCS_SCATT(rcs_scat4,ref_dist, dist_scat4, ang_scat4, rot_l, rot_h, rot_step,k,len_dipole);


rcs_scat5=1;
ang_scat5=-180;
dist_scat5=3;

e_scat5 = RCS_SCATT(rcs_scat5,ref_dist, dist_scat5, ang_scat5, rot_l, rot_h, rot_step,k,len_dipole);


rcs_scat6=1;
ang_scat6=160;
dist_scat6=3;

e_scat6 = RCS_SCATT(rcs_scat6,ref_dist, dist_scat6, ang_scat6, rot_l, rot_h, rot_step,k,len_dipole);





E_scats=e_scat1 + e_scat2 + e_scat3 + e_scat4 + e_scat5 + e_scat6;





l=size(E_scats);

for n=1:l(2)
    E_scats(n)=abs(E_scats(n));
end;




% TOTAL

l=size(E_scats);
for n=1:l(2)
    E_final(n)=E_plates(n)+E_scats(n);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part-2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumptions:
%
%   1. Plane height is not taken in account. All objects are at same
%      height.
%   2. Radar is having gaine 1 in all directions.








%   Plane Motion

% Radar specs

Rrep=1*10^2;     % Hz
Tp=10^-6;           % Sec- pulse duration
Trep=1/Rrep;
R_ht=10;        % Radar height



% Plane specs

P_dist= 50*10^3;        % distance of plane
P_ht= 50;               % hight of plane
mark=350;               % one mark velocity in meters/sec
P_velo=5*mark;          % plane speed
P_angl=30;



% motion simulation specs

Tsim= 4;      % simulation time in seconds
delta_angl=1;  % change in angle w.r.t. radar possible in  pulse-to-pulse in degrees


         

% clutter definitions

C_dist=[35000, 35500, 36000, 34000,41000, 43000, 40000, 45000];    %distance of trees in meters

C_num=8;        % number of trees

     
C_RCS=[1.5, 1.5, 1.6, 1.7, 1.6, 1.1, 1, 1.6]; % in square meters



NxtTx=0.0;
NxtRx=Tsim;
RealTm=0.0;
LstEvnt=0.0;


C=3*10^8;

Ievnt=1;    % increases as event happens
Irfl=1;     % pointing to empty location
Ip=1;       % for plane refln only



dist=P_dist;

Gr_p=1;    %gain in direction of plane


% Plane Tragectories:

% NOTE: For changing Trajectory Un-comment apropriate following lines.

V_traj=TrajRay(Tsim,Trep,P_velo,P_dist);	% for  Rayleigh Trajectory

%V_traj=TrajSin(Tsim,Trep,P_velo,P_dist);	% for Sine Trajectory
%V_traj=TrajExp(Tsim,Trep,P_velo,P_dist);	% for Exponential Trajectory

Tx_num=0;
 


 while(RealTm<(Tsim-Trep))
            
            if(NxtTx < NxtRx)           % if nxt event is Tx
                
             
                RealTm=NxtTx;           % update time
                
                NxtTx=RealTm + Trep;    % Update NxtTx
                
                Tx_num=Tx_num+1;
                
                % for plane angle
                
                P_angl=P_angl+delta_angl*RANDOM();
                
                P_ang_eff=ceil(P_angl);
                
                
                if (abs(P_angl)>=(rot_h-1))
                    
                    P_angl=p1_angle*RANDOM();
                end
                
                    
                
                % events happened              
                               
                Tmeet=RealTm+dist/(C+V_traj(Tx_num));   % t at which wave and plane meet
                Dmeet=C*(Tmeet-RealTm);         % distance from radar where wave and plane meet
                Trec=Tmeet+Dmeet/C;             % time at which wave receves bk to radar
                Drec=Dmeet-(Trec-Tmeet)*V_traj(Tx_num);
                
                            
                
                
                               
                % Getting RCS of Plane depending on its angle.                               
                
                [min_diff, min_idx]=min(abs(x-P_ang_eff));  % min_idx contains index of which rcs to be taken
                
                P_RCS=E_final(min_idx); % Plane RCS seen by radar
                
                PathLoss=1/(Dmeet^4);    % Path loss
                
                RCSrec=P_RCS*PathLoss;
                
                
                
                
                
                dist=dist-Trep*V_traj(Tx_num);  % update distance
                
                % Reflection 2D Vector --> RFL;
                
                RFL(1,Irfl)=RCSrec;
                RFL(2,Irfl)=Trec;
               
                Irfl=Irfl+1;
                
                % Adding Clutter Reflections...
                
                
                for i=1:C_num
                     
                    C_Tmeet=RealTm+C_dist(i)/C;
                    C_Dmeet=C_dist(i);
                    C_Trec=C_Tmeet+C_Dmeet/C;
                   
                    C_PathLoss=1/(C_Dmeet^4);    % Path loss
             
                    C_RCSrec=C_RCS(i)*C_PathLoss;
                    
                    
                    RFL(1,Irfl)=C_RCSrec;
                    RFL(2,Irfl)=C_Trec;
                                 
                    Irfl=Irfl+1;
                              
                                        
                end
                
                
                
                RESULT(1,Ievnt)=1;              % Transmitted power
                RESULT(2,Ievnt)=0;              % Received Power
                RESULT(3,Ievnt)=RealTm;         % Real time
                
                Ievnt=Ievnt+1;                  % increase event

                
                NxtRx=min(RFL(2,:));
               
                
               
               
                   
            else                % if nxt event is Rx     
                
                RealTm=NxtRx;
               
                
                
                % remove entry that is being received
                
                if(isempty(RFL)==0)
                    [v_remove,remove]=min(abs(RFL(2,:)-RealTm));
                    
                    Thres=C*(Trep-(NxtTx-RealTm))/2;
                    
                
                    RESULT(1,Ievnt)=0;              % Tx
                    RESULT(2,Ievnt)=RFL(1,remove)*(Thres^4);  % Rx
                    RESULT(3,Ievnt)=RealTm;         % Time
                    
                    Ievnt=Ievnt+1;                  % event increase
                
                                                    
                    RFL(:,remove)=[];       % Remove entry
                    
                    Irfl=Irfl-1;            % update size of RFL
                    
                end
                
                
                
                
                % Next reception event
                
                if(isempty(RFL))
                    NxtRx=Tsim;
                else
                    NxtRx=min(RFL(2,:));
                end
                
                
                
                                    
                
                                
             end
             
 end
            
     
        
 
 
 
        Tm= RESULT(3,:);    % Time
        RES=RESULT;         % Rx values
             
        RES(3,:)=[];
        
      
       % just scaling for better viewing of plot(ZOOMing)
        
        mul=max(RES(2,:))*2;
        RES(1,:)=mul*RES(1,:);
        
        
        
        
        
        figure;
        stem(Tm, RES');
        title('Sent and Received pulses (After 1/R^2 thresholding)');
        legend('Tx Pulses','Rx Pulses',2);
        xlabel('Time(in seconds)');
        ylabel('Amplitude(RCS/R^4)');
        
          
       
       
            

        
        figure;
       
        stem(Tm,RESULT(2,:));
        title('Received Signal (After 1/R^2 Thresholding)');
        legend('Rx Pulses',1);
        xlabel('Time(Sec.)');
        ylabel('Amplitude(RCS/R^4)');
        
       
 %%%%%%%%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%%       
        
        
 %%%%%%%%%%%%%% 3/5 (Max) %%%%%%%%%%%%%%%%%%%%%%
 
        
        
        % Removing last ZEROS
        
        s=size(RESULT);
        T_RESULT=fliplr(RESULT);
        
        [m, k]=min(abs(T_RESULT - 1));
        
        for i=1:(k)
            T_RESULT(:,i)=[];
        end
        
        F_RESULT=fliplr(T_RESULT);
        
        RESULT=F_RESULT;
        
        s=size(RESULT);
        
        
        % collecting signal received in one pulse duration
        
        
        k=1;
        for i=1:(C_num+2):s(2)
            for j=0:C_num+1
                CHOP(k,(j+1))=RESULT(2,(i+j));
                TMM(k,(j+1))=RESULT(3,(i+j));
                
            end
            k=k+1;
        end
        
        
        ss=size(CHOP);
        
        

        % Finding Max peak:

        for i=1:ss(1)
    
            [CHOP_val(i), INX]=max(CHOP(i,:));
            CHOP_Tm(i)=TMM(i,INX);
            
    
        end
        
         
        Tx_Tm=0:Trep:(Tsim-Trep);
        
        
        figure;
        stem(CHOP_Tm,CHOP_val);
        title('Filtered Received signal(3/5 or max)');
        xlabel('Time(Sec.)');
        ylabel('RCS');
        
        
        
        
        
        %%%%%%%%%%%%3/5%%%%%%%%%%%%%%%%%%%%%
        
        Th=max(C_RCS)+0.6;
        
        s=size(CHOP_val)
        
        sz_tf=floor(s(2)/5)
        count=0;
        for i=1:(sz_tf-1)
            count=0;
            for j=0:4
                
                if(CHOP_val(i+j)>Th)
                    count=count+1;
                end
                
                temp(3,(j+1))=CHOP_val(i+j);
                temp(2,(j+1))=CHOP_Tm(i+j);
                temp(1,(j+1))=Tx_Tm(i+j);
            end
            
            if(count>=3)
                [sampl_val, sampl_inx]=max(temp(1,j));
                TF(1,i)=temp(1,sampl_inx);
                TF(2,i)=temp(2,sampl_inx);
                TF(3,i)=temp(3,sampl_inx);
            else
                [sampl_val, sampl_inx]=max(temp(1,j));
                TF(1,i)=temp(1,sampl_inx);
                TF(2,i)=temp(2,sampl_inx);
                TF(3,i)=0;
            end
        end
        
        figure;
        stem(TF(3,:));
                
        
        
        CHOP_val=TF(3,:);
        CHPO_Tm=TF(2,:);
        
        Tx_Tm=TF(1,:);
        
        
        % calculating Distance and velocity of plane:
        
        
        % distance -Actual detected
        
        s=size(CHOP_val);
        
        for i=1:s(2)
            
            Est_dist(i)=(CHOP_Tm(i)-Tx_Tm(i))*C/2;
            
        end
        
        
        
        
        
        % shift in time -actual
        
        for i=1:s(2)
            
            Est_sh(i)=(CHOP_Tm(i)-Tx_Tm(i));
        end
        
        Vmax=max(V_traj);
        
        delta_dist=Vmax*Trep;
        
        delta_shift=delta_dist*2/C;
        
        
        % correction in estimated time(WINDOWING)
        
        
        Est_delta_sh=diff(Est_sh);  % shift bteween Rx pulses 
        
        s=size(Est_delta_sh);
        
        for i=1:s(2)
            
            if(abs(Est_delta_sh(i))>delta_shift)
                Est_delta_sh(i)=0;
            end
        end
        
        
        Est_sh_corr(1)=Est_sh(1);
                
        
        for i=1:s(2)
            
            Est_sh_corr(i+1)=Est_sh_corr(i)+Est_delta_sh(i);
        end
        
        
        
        % corrected dist;
        
        for i=1:(s(2)+1)
            
            Est_dist_corr(i)=Est_sh_corr(i)*C/2;
        end
        
        
        
        
        % inst_velo
        
            
        Est_inst_velo = -1*diff(Est_dist_corr)/Trep;
        
        Temp_Tm=CHOP_Tm;
        s=size(Temp_Tm);
        Temp_Tm(s(2))=[];
        
        figure;
        subplot(1,2,1);
        plot(Est_inst_velo,'rx');
        title('Instantaneous Velocity Traj(3/5 or max)');
        ylabel('Velocity(m/s)');
        xlabel('Time(Sec.)');
        
       
        subplot(1,2,2);
        plot(Est_dist,'rx');
        title('Instantanepus Distance Traj(3/5 or max)');
        ylabel('Distance(m)');
        xlabel('Time(Sec.)');
       
        
        
        % coreecting inst_velo
        
        s=size(Est_inst_velo);

        
        
        
        for i=2:s(2)
            
            if(Est_inst_velo(i) == 0)
                Est_inst_velo(i)=Est_inst_velo(i-1);
            else
               Est_inst_velo(i)=Est_inst_velo(i);
            end
            
        end
         
        
        
        
        
        figure;
        subplot(1,2,1)
        plot( Est_inst_velo);
        Title('Estimated Velocity Traj(3/5 or max)');
        xlabel('Time(Sec.)');
        ylabel('Velocity(m/s)');
        
        
       
        subplot(1,2,2)
        plot(Est_dist_corr);
        Title('Estimated Distance Traj(3/5 or max)');
        xlabel('Time(Sec.)');
        ylabel('distance(m)');
        
        
%         %%%%%%%%%%% Threshold detection %%%%%%%%%%%%
%         
%         Th=max(C_RCS)+0.6;
%         
%         s=size(RESULT(2,:));
%         
%         for i=1:s(2)
%             if(RESULT(2,i)<=Th)
%                 RES_th(i)=0;
%             else
%                 RES_th(i)=RESULT(2,i);
%             end
%         end
%         
%         %figure;
%         %stem(RESULT(3,:),RES_th);
%           
%         % Finding Maax;
%         
%         s=size(RES_th);
%         k=1;
%         for i=1:(C_num+2):s(2)
%             for j=0:C_num+1
%                 CHOP(k,(j+1))=RES_th((i+j));
%                 TMM(k,(j+1))=RESULT(3,(i+j));
%                 
%             end
%             k=k+1;
%         end
%         
%         
%         ss=size(CHOP);
% 
% 
%         for i=1:ss(1)
%     
%             [CHOP_val_th(i), INX]=max(CHOP(i,:));
%             CHOP_Tm_th(i)=TMM(i,INX);
%             
%     
%         end
%         
%         Tx_Tm=0:Trep:(Tsim-Trep);
%         
%         
%         figure;
%         stem(CHOP_Tm_th,CHOP_val_th);
%         title('Filtered Received signal(Threshold)');
%         xlabel('Time(Sec.)');
%         ylabel('RCS');
%         
%         
%         % distance -Actual detected
%         
%         s=size(CHOP_val_th);
%         
%         for i=1:s(2)
%             
%             Est_dist(i)=(CHOP_Tm_th(i)-Tx_Tm(i))*C/2;
%             
%         end
%         
%         % shift in time -actual
%         
%         for i=1:s(2)
%             
%             Est_sh(i)=(CHOP_Tm_th(i)-Tx_Tm(i));
%         end
%         
%         Vmax=max(V_traj);
%         
%         delta_dist=Vmax*Trep;
%         
%         delta_shift=delta_dist*2/C;
%         
%         
%         % correction in estimated time(WINDOWING)
%         
%         
%         Est_delta_sh=diff(Est_sh);  % shift bteween Rx pulses 
%         
%         s=size(Est_delta_sh);
%         
%         for i=1:s(2)
%             
%             if(abs(Est_delta_sh(i))>delta_shift)
%                 Est_delta_sh(i)=0;
%             end
%         end
%         
%         
%          
%         Est_sh_corr(1)=Est_sh(1);
%                 
%         
%         for i=1:s(2)
%             
%             Est_sh_corr(i+1)=Est_sh_corr(i)+Est_delta_sh(i);
%         end
%         
%         
%         
%         % corrected dist;
%         
%         for i=1:(s(2)+1)
%             
%             Est_dist_corr(i)=Est_sh_corr(i)*C/2;
%         end
%         
%          % inst_velo
%         
%             
%         Est_inst_velo = -1*diff(Est_dist_corr)/Trep;
%         
%         Temp_Tm=CHOP_Tm;
%         s=size(Temp_Tm);
%         Temp_Tm(s(2))=[];
%         
%         figure;
%         subplot(1,2,1);
%         plot(Temp_Tm,Est_inst_velo,'rx');
%         title('Instantaneous Velocity Traj(Threshold)');
%         ylabel('Velocity(m/s)');
%         xlabel('Time(Sec.)');
%         
%        
%         subplot(1,2,2);
%         plot(CHOP_Tm_th,Est_dist,'rx');
%         title('Instantanepus Distance Traj(Threshold)');
%         ylabel('Distance(m)');
%         xlabel('Time(Sec.)');
%        
%         
%         % coreecting inst_velo
%         
%         s=size(Est_inst_velo);
% 
%         
%         
%         
%         for i=3:s(2)
%             
%             if(Est_inst_velo(i) == 0)
%                 Est_inst_velo(i)=(2*Est_inst_velo(i-1)+Est_inst_velo(i-2))/3;
%             else
%                Est_inst_velo(i)=Est_inst_velo(i);
%             end
%             
%         end
%          
%         
%         
%         
%         
%         
%         
%         figure;
%         subplot(1,2,1)
%         plot(Temp_Tm, Est_inst_velo);
%         Title('Estimated Velocity Traj(Threshold)');
%         xlabel('Time(Sec.)');
%         ylabel('Velocity(m/s)');
%         
%         
%        
%         subplot(1,2,2)
%         plot(CHOP_Tm_th,Est_dist_corr);
%         Title('Estimated Distance Traj(Threshold)');
%         xlabel('Time(Sec.)');
%         ylabel('distance(m)');
%         
%         
%         %%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%