%%
clc; clear all;
close all;

global u1 u2 u3
global val_weights_1 val_weights_2 val_weights_3 f_1 f_2 f_3 g_1 g_2 g_3 theta_2_star theta_3_star
r=20;

h=0.01;

monte_check=3;
init_val_weights_1 = 0.01*rand(1,6);
init_val_weights_2 = 0.01*rand(1,6);
init_val_weights_3 = 0.01*rand(1,6);
for monte=1:monte_check
    x_init = 0.1*rand(1,1)+0.5*[(40*pi)/180 0 (-3*pi)/180 0 (3*pi)/180 0];
    %      x_init = 0.4*[1.0468    0.3487   -0.2963    0.3487    0.4011    0.3487];
    Reward_1(1)=0;
    Reward_2(1)=0;
    Reward_3(1)=0;
    x_store_1=zeros(100,6);x_store_2=zeros(100,6);x_store_3=zeros(100,6);
    
    timer_val1=uint64(0);timer_val2=uint64(0);timer_val3=uint64(0);
    last_event1 = 0;last_event2 = 0;last_event3 = 0;event2=0;event1=0;event3=0;
    time_base1 = [];time_base2 = [];time_base3 = [];
    current_event1 = 0;current_event2 = 0;current_event3 = 0;
    inter_event_time1 = [];inter_event_time2 = [];inter_event_time3 = [];
    
    x0 = x_init;
    val_weights_1 = init_val_weights_1;
    val_weights_2 = init_val_weights_2;
    val_weights_3 = init_val_weights_3;
    
    x_store_1=zeros(100,6);x_store_2=zeros(100,6);x_store_3=zeros(100,6);
    event_threshold_1(1) = zeros;event_counter_1=0;
    event_indicating_function_1(1) = [x0(1) x0(2)]*[x0(1);x0(2)];
    event_threshold_2(1) = zeros;event_counter_2=0;
    event_indicating_function_2(1) = [x0(3) x0(4)]*[x0(3);x0(4)];
    event_threshold_3(1) = zeros;event_counter_3=0;
    event_indicating_function_3(1) = [x0(5) x0(6)]*[x0(5);x0(6)];
    
    monte
    
    if (monte==1)
               alpha_1 = 2.5;beta_1 = 0.15;
                alpha_2 = 5;beta_2 = 0.25;
                alpha_3 = 5;beta_3 = 0.25;
                Ru_1=0.5*0.08;Ru_11=-0.05*0.05;Qx_1=20*0.0001;
                Ru_2=0.5*0.02;Ru_22=-0.050*0.02;Qx_2=15*0.005;
                Ru_3=0.5*0.02;Ru_33=-0.01*0.02;Qx_3=15*0.005;
                gamma_1=1.086;gamma_2=1.086;gamma_3=1.086;
            elseif (monte==3)
        alpha_1 = 2.5;beta_1 = 0.15;
        alpha_2 = 5;beta_2 = 0.25;
        alpha_3 = 5;beta_3 = 0.25;
        Ru_1=0.5*0.2;Ru_11=-0.05*0.05;Qx_1=30*0.0001;
        Ru_2=0.5*0.02;Ru_22=-0.05*0.02;Qx_2=15*0.005;
        Ru_3=0.5*0.02;Ru_33=-0.01*0.02;Qx_3=15*0.005;
        gamma_1=1.086;gamma_2=1.086;gamma_3=1.086;
    elseif (monte==2)
        alpha_1 = 2.5;beta_1 = 0.15;
        alpha_2 = 5;beta_2 = 0.25;
        alpha_3 = 5;beta_3 = 0.25;
        Ru_1=0.5*0.2;Ru_11=-0.05*0.05;Qx_1=30*0.0001;
        Ru_2=0.5*0.02;Ru_22=-0.05*0.02;Qx_2=15*0.005;
        Ru_3=0.5*0.02;Ru_33=-0.01*0.02;Qx_3=15*0.005;
        gamma_1=0*1.086;gamma_2=0*1.086;gamma_3=0*1.086;
    end
    u1=0;u2=0;u3=0;
    theta_2_star = 0;-0.227;
    theta_3_star = 0;0.559;
    p2 = 6070;
    p3 = 192;
    X(1,:)=x0;
    x1=0; x2=0; x3=0; x4=0;x5=0; x6=0;
    temp_flag_1=1;temp_flag_2=1;temp_flag_3=1;
    for k=1:r/h
        event_threshold_1(k) = gamma_1*(2-exp(toc(timer_val1)))*[x1 x2]*[x1;x2];
        event_indicating_function_1(k) = [x0(1) x0(2)]*[x0(1) x0(2)]';
        if k==1 || event_indicating_function_1(k) >= event_threshold_1(k)
            last_event1=current_event1;
            flag_event_1=1;
            event_counter_1 = event_counter_1+1;tic;
            timer_val1=tic;
            u_h1=u1;
            event1=event1+1;
            current_event1=k;
            time_base1=[time_base1,current_event1];
            inter_event_time1 = [inter_event_time1,current_event1-last_event1];
        else
              flag_event_1=0;
            i_hyb_count_1 = 0;
        end
        event_threshold_2(k) = gamma_2*(2-exp(toc(timer_val2)))*[(x3-theta_2_star) x4]*[(x3-theta_2_star);x4];
        event_indicating_function_2(k) = [x0(3)-theta_2_star x0(4)]*[x0(3)-theta_2_star x0(4)]';
        if k==1 || event_indicating_function_2(k) >= event_threshold_2(k)
            last_event2=current_event2;
            flag_event_2=1;
            event_counter_2 = event_counter_2+1;tic;
            timer_val2=tic;
            u_h2=u2;
            event2=event2+1;
            current_event2=k;
            time_base2=[time_base2,current_event2];
            inter_event_time2 = [inter_event_time2,current_event2-last_event2];
        else
            flag_event_2=0;
            i_hyb_count_2 = 0;
        end
        event_threshold_3(k) = gamma_3*(2-exp(toc(timer_val3)))*[x5-theta_3_star x6]*[x5-theta_3_star;x6];
        event_indicating_function_3(k) = [x0(5)-theta_3_star x0(6)]*[x0(5)-theta_3_star x0(6)]';
        if k==1 || event_indicating_function_3(k) >= event_threshold_3(k)
            last_event3=current_event3;
            flag_event_3=1;
            event_counter_3 = event_counter_3+1;tic;
            timer_val3=tic;
            u_h3=u3;
            event3=event3+1;
            current_event3=k;
            time_base3=[time_base3,current_event3];
            inter_event_time3 = [inter_event_time3,current_event3-last_event3];
        else
            flag_event_3=0;
            i_hyb_count_3 = 0;
        end
        
        [t,y] = ode23('e_closed_dyn_ex_biped',[0 h],x0);
        X(k+1,:) = y(end,:);
        x0=X(k+1,:);
        
        if flag_event_1==1
            x1=x0(1);x2=x0(2);
            Hybrid_evolve_1(k)=0;
        elseif flag_event_1==0
            x1=X(k,1);x2=X(k,2);
            Hybrid_count_1 = 0;
            x_store_1(temp_flag_1,:) = [x0(1);x0(2);x0(3);x0(4);x0(5);x0(6)];
            u_store_1(temp_flag_1)=u1;
            temp_flag_1=temp_flag_1+1;
        end
        if flag_event_2==1
            x3=x0(3);x4=x0(4);
            Hybrid_evolve_2(k)=0;
        elseif flag_event_2==0
            x3=X(k,3);x4=X(k,4);
            Hybrid_count_2 = 0;
            x_store_2(temp_flag_2,:) = [x0(1);x0(2);x0(3);x0(4);x0(5);x0(6)];
            u_store_2(temp_flag_2)=u2;
            temp_flag_2=temp_flag_2+1;
        end
        if flag_event_3==1
            x5=x0(5);x6=x0(6);
            Hybrid_evolve_3(k)=0;
        elseif flag_event_3==0
            x5=X(k,5);x6=X(k,6);
            Hybrid_count_3 = 0;
            x_store_3(temp_flag_3,:) = [x0(1);x0(2);x0(3);x0(4);x0(5);x0(6)];
            u_store_3(temp_flag_3)=u3;
            temp_flag_3=temp_flag_3+1;
        end
        
        
        
        
        Q_1(k)=Qx_1*x1*x1+10*Qx_1*x2*x2+(Qx_2/30)*(x3-theta_2_star)*(x3-theta_2_star)+10*(Qx_2/30)*x4*x4+(Qx_3/30)*(x5-theta_3_star)*(x5-theta_3_star)+10*(Qx_3/30)*x6*x6;
        R_1(k)= u1*Ru_1*u1;
        R_1_h(k)=u_h1*Ru_11*u_h1;
        Q_2(k)=Qx_2*(x3-theta_2_star)*(x3-theta_2_star)+10*Qx_2*x4*x4+(Qx_1/30)*x1*x1+10*(Qx_1/30)*x2*x2+(Qx_3/30)*(x5-theta_3_star)*(x5-theta_3_star)+10*(Qx_3/30)*x6*x6;
        R_2(k)= u2*Ru_2*u2;
        R_2_h(k)=u_h2*Ru_22*u_h2;
        Q_3(k)=Qx_3*(x5-theta_3_star)*(x5-theta_3_star)+10*Qx_3*x6*x6+(Qx_1/30)*x1*x1+10*(Qx_1/30)*x2*x2+(Qx_2/30)*(x3-theta_2_star)*(x3-theta_2_star)+10*(Qx_2/30)*x4*x4;
        R_3(k)= u3*Ru_3*u3;
        R_3_h(k)=u_h3*Ru_33*u_h3;
        V_1(k)=Q_1(k)+R_1(k)+R_2_h(k)+R_3_h(k);
        Reward_1(k+1)=Reward_1(k)+V_1(k);
        V_2(k)=Q_2(k)+R_2(k)+R_1_h(k)+R_3_h(k);
        Reward_2(k+1)=Reward_2(k)+V_2(k);
        V_3(k)=Q_3(k)+R_3(k)+R_1_h(k)+R_2_h(k);
        Reward_3(k+1)=Reward_3(k)+V_3(k);
        
        
        
        %         x1=x0(1);x2=x0(2);x3=x0(3);x4=x0(4);x5=x0(5);x6=x0(6);
        
        
        % theta_2_star = -0.559;
        % theta_3_star = 0.226;
        % p2 = 226;
        % p3 = 5240;
        f1= x2;f3=x4;f5=x6;
        f2=0.1*(1-5.25*(x1)^2)*x2-x1+0.01*x6*(x5-theta_3_star)+0.01*x4*(x3-theta_2_star);
        f4=(0.01*(1-p2*(x3-theta_2_star)^2)*x4)-(4*(x3-theta_2_star))+(0.057*x1*x2)+(0.1*(x4-x6));
        f6=(0.01*(1-p3*(x5-theta_3_star)^2)*x6)-(4*(x5-theta_3_star))+(0.057*x1*x2)+(0.1*(x6-x4));
        f_1=[f1;f2];g_1=[0;1];
        f_2=[f3;f4];g_2=[0;1];
        f_3=[f5;f6];g_3=[0;1];
        f = [f1;f2;f3;f4;f5;f6];
        g = [0;1;0;1;0;1];
        
        phi_x_1 = ([x1^2 x2^2 (x1*x3*x4) (x2*x3*x4) (x1*x5*x6) (x2*x5*x6)]);
        phi_x_2 = ([x3^2 x4^2 (x2*x1*x3) (x2*x1*x4) (x6*x4) x6*x3]);
        phi_x_3 = ([x5^2 x6^2 (x1*x2*x5) (x1*x2*x6) (x4*x5) (x4*x6)]);
        del_phi_x_1 =[diag([x1*2 x2*2]) diag([x3*x4 x3*x4]) diag([x5*x6 x5*x6])]; % phi_x_1*(eye(length(phi_x_1))-phi_x_1)%dtansig(phi_x_1,([x1^2 x2^2 (x1*x3*x4) (x2*x3*x4) (x1*x5*x6) (x2*x5*x6)]));
        del_phi_x_2 = [diag([x3*2 x4*2]) diag([x2*x1 x2*x1]) diag([x6 x6])]; %phi_x_2*(eye(length(phi_x_1))-phi_x_2);%dtansig(phi_x_2,([x3^2 x4^2 (x2*x1*x3) (x2*x1*x4) (x6*x4) x6*x3]));
        del_phi_x_3 = [diag([x5*2 x6*2]) diag([x2*x1 x2*x1]) diag([x4 x4])];% phi_x_3*(eye(length(phi_x_1))-phi_x_3);%dtansig(phi_x_3,([x5^2 x6^2 (x1*x2*x5) (x1*x2*x6) (x4*x5) (x4*x6)]));

        V_estimate_1(k) = val_weights_1*phi_x_1';
        V_estimate_2(k) = val_weights_2*phi_x_2';
        V_estimate_3(k) = val_weights_3*phi_x_3';
        del_V_estimate_1 = (val_weights_1*del_phi_x_1');
        del_V_estimate_2 = (val_weights_2*del_phi_x_2');
        del_V_estimate_3 = (val_weights_3*del_phi_x_3');
        
        
        Hamiltonian_1(k) = Q_1(k)+R_2_h(k)+R_3_h(k)+ del_V_estimate_1*f_1-0.25*del_V_estimate_1*g_1*Ru_1^-1*g_1'*del_V_estimate_1';
        Hamiltonian_2(k) = Q_2(k)+R_1_h(k)+R_3_h(k)+ del_V_estimate_2*f_2-0.25*del_V_estimate_2*g_2*Ru_2^-1*g_2'*del_V_estimate_2';
        Hamiltonian_3(k) = Q_3(k)+R_2_h(k)+R_1_h(k)+ del_V_estimate_3*f_3-0.25*del_V_estimate_3*g_3*Ru_3^-1*g_3'*del_V_estimate_3';
        del_Hamiltonian_1 = del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*val_weights_1';
        del_Hamiltonian_2 = del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*val_weights_2';
        del_Hamiltonian_3 = del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*val_weights_3';
        
        
        if flag_event_1==1
            %             alpha_1=25;
            a_val_weights_1 = h*(-((alpha_1*del_Hamiltonian_1*Hamiltonian_1(k))/(del_Hamiltonian_1'*del_Hamiltonian_1+1)^2)+0.5*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            b_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')*(Q_1(k)+ ((val_weights_1+a_val_weights_1/2)*del_phi_x_1')*f_1-0.25*((val_weights_1+a_val_weights_1/2)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+a_val_weights_1/2)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')+1)^2)+0.5*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            c_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')*(Q_1(k)+ ((val_weights_1+b_val_weights_1/2)*del_phi_x_1')*f_1-0.25*((val_weights_1+b_val_weights_1/2)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+b_val_weights_1/2)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')+1)^2)+0.5*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            d_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')*(Q_1(k)+ ((val_weights_1+c_val_weights_1)*del_phi_x_1')*f_1-0.25*((val_weights_1+c_val_weights_1)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+c_val_weights_1)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')+1)^2)+0.5*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            
            val_weights_1 = val_weights_1+a_val_weights_1/6+b_val_weights_1/3+c_val_weights_1/3+d_val_weights_1/6;
            plot_val_weights_1(k,:) = (val_weights_1);
            u1 = -inv(Ru_1)*g_1'*del_V_estimate_1'+randn*0.00;
        elseif flag_event_1==0 && i_hyb_count_1 <30
            %             alpha_1=5;
            a_val_weights_1 = h*(-((alpha_1*del_Hamiltonian_1*Hamiltonian_1(k))/(del_Hamiltonian_1'*del_Hamiltonian_1+1)^2)+0.0*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            b_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')*(Q_1(k)+ ((val_weights_1+a_val_weights_1/2)*del_phi_x_1')*f_1-0.25*((val_weights_1+a_val_weights_1/2)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+a_val_weights_1/2)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+a_val_weights_1/2)')+1)^2)+0.0*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            c_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')*(Q_1(k)+ ((val_weights_1+b_val_weights_1/2)*del_phi_x_1')*f_1-0.25*((val_weights_1+b_val_weights_1/2)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+b_val_weights_1/2)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+b_val_weights_1/2)')+1)^2)+0.0*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            d_val_weights_1 = h*(-((alpha_1*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')*(Q_1(k)+ ((val_weights_1+c_val_weights_1)*del_phi_x_1')*f_1-0.25*((val_weights_1+c_val_weights_1)*del_phi_x_1')*g_1*Ru_1^-1*g_1'*((val_weights_1+c_val_weights_1)*del_phi_x_1')'))/((del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')'*(del_phi_x_1'*f_1-0.5*del_phi_x_1'*g_1*Ru_1^-1*g_1'*del_phi_x_1*(val_weights_1+c_val_weights_1)')+1)^2)+0.0*beta_1*del_phi_x_1'*g_1*Ru_1^-1*g_1'*[x1;x2])';
            
            val_weights_1 = val_weights_1+a_val_weights_1/6+b_val_weights_1/3+c_val_weights_1/3+d_val_weights_1/6;
            i_hyb_count_1=i_hyb_count_1+1;
        end
        if flag_event_2==1
            a_val_weights_2 = h*(-((alpha_2*del_Hamiltonian_2*Hamiltonian_2(k))/(del_Hamiltonian_2'*del_Hamiltonian_2+1)^2)+0.5*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            b_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')*(Q_2(k)+ ((val_weights_2+a_val_weights_2/2)*del_phi_x_2')*f_2-0.25*((val_weights_2+a_val_weights_2/2)*del_phi_x_2')*g_2*Ru_2^-1*g_2'*((val_weights_2+a_val_weights_2/2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')+1)^2)+0.5*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            c_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')*(Q_2(k)+ ((val_weights_2+b_val_weights_2/2)*del_phi_x_2')*f_2-0.25*((val_weights_2+b_val_weights_2/2)*del_phi_x_2')*g_2*Ru_2^-1*g_2'*((val_weights_2+b_val_weights_2/2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')+1)^2)+0.5*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            d_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')*(Q_2(k)+   ((val_weights_2+c_val_weights_2)*  del_phi_x_2')*f_2-0.25*((val_weights_2+c_val_weights_2)*del_phi_x_2')  *g_2*Ru_2^-1*g_2'*((val_weights_2+c_val_weights_2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')+1)^2)+0.5*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            
            val_weights_2 = val_weights_2+a_val_weights_2/6+b_val_weights_2/3+c_val_weights_2/3+d_val_weights_2/6;
            plot_val_weights_2(k,:) = (val_weights_2);
            u2 = -inv(Ru_2)*g_2'*del_V_estimate_2'+randn*0.00;
        elseif flag_event_2==0 && i_hyb_count_2 <25
            a_val_weights_2 = h*(-((alpha_2*del_Hamiltonian_2*Hamiltonian_2(k))/(del_Hamiltonian_2'*del_Hamiltonian_2+1)^2)+0.0*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            b_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')*(Q_2(k)+ ((val_weights_2+a_val_weights_2/2)*del_phi_x_2')*f_2-0.25*((val_weights_2+a_val_weights_2/2)*del_phi_x_2')*g_2*Ru_2^-1*g_2'*((val_weights_2+a_val_weights_2/2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+a_val_weights_2/2)')+1)^2)+0.0*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            c_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')*(Q_2(k)+ ((val_weights_2+b_val_weights_2/2)*del_phi_x_2')*f_2-0.25*((val_weights_2+b_val_weights_2/2)*del_phi_x_2')*g_2*Ru_2^-1*g_2'*((val_weights_2+b_val_weights_2/2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+b_val_weights_2/2)')+1)^2)+0.0*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            d_val_weights_2 = h*(-((alpha_2*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')*(Q_2(k)+   ((val_weights_2+c_val_weights_2)*  del_phi_x_2')*f_2-0.25*((val_weights_2+c_val_weights_2)*del_phi_x_2')  *g_2*Ru_2^-1*g_2'*((val_weights_2+c_val_weights_2)*del_phi_x_2')'))/((del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')'*(del_phi_x_2'*f_2-0.5*del_phi_x_2'*g_2*Ru_2^-1*g_2'*del_phi_x_2*(val_weights_2+c_val_weights_2)')+1)^2)+0.0*beta_2*del_phi_x_2'*g_2*Ru_2^-1*g_2'*[x3;x4])';
            
            val_weights_2 = val_weights_2+a_val_weights_2/6+b_val_weights_2/3+c_val_weights_2/3+d_val_weights_2/6;
            i_hyb_count_2=i_hyb_count_2+1;
        end
        if flag_event_3==1
            a_val_weights_3 = h*(-((alpha_3*del_Hamiltonian_3*Hamiltonian_3(k))/(del_Hamiltonian_3'*del_Hamiltonian_3+1)^2)+0.5*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            b_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')*(Q_3(k)+ ((val_weights_3+a_val_weights_3/2)*del_phi_x_3')*f_3-0.25*((val_weights_3+a_val_weights_3/2)*del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+a_val_weights_3/2)*del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')'*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')+1)^2)+0.5*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            c_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')*(Q_3(k)+ ((val_weights_3+b_val_weights_3/2)*del_phi_x_3')*f_3-0.25*((val_weights_3+b_val_weights_3/2)*del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+b_val_weights_3/2)*del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')'*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')+1)^2)+0.5*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            d_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')*(Q_3(k)+   ((val_weights_3+c_val_weights_3)*del_phi_x_3')*  f_3-0.25*((val_weights_3+c_val_weights_3)*  del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+c_val_weights_3)*  del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')'  *(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')+1)^2)+  0.5*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            
            val_weights_3 = val_weights_3+a_val_weights_3/6+b_val_weights_3/3+c_val_weights_3/3+d_val_weights_3/6;
            plot_val_weights_3(k,:) = (val_weights_3);
            u3 = -inv(Ru_3)*g_3'*del_V_estimate_3'+randn*0.00;
        elseif flag_event_3==0 && i_hyb_count_3 <30
            a_val_weights_3 = h*(-((alpha_3*del_Hamiltonian_3*Hamiltonian_3(k))/(del_Hamiltonian_3'*del_Hamiltonian_3+1)^2)+0.0*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            b_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')*(Q_3(k)+ ((val_weights_3+a_val_weights_3/2)*del_phi_x_3')*f_3-0.25*((val_weights_3+a_val_weights_3/2)*del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+a_val_weights_3/2)*del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')'*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+a_val_weights_3/2)')+1)^2)+0.0*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            c_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')*(Q_3(k)+ ((val_weights_3+b_val_weights_3/2)*del_phi_x_3')*f_3-0.25*((val_weights_3+b_val_weights_3/2)*del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+b_val_weights_3/2)*del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')'*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+b_val_weights_3/2)')+1)^2)+0.0*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            d_val_weights_3 = h*(-((alpha_3*(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')*(Q_3(k)+   ((val_weights_3+c_val_weights_3)*del_phi_x_3')*  f_3-0.25*((val_weights_3+c_val_weights_3)*  del_phi_x_3')*g_3*Ru_3^-1*g_3'*((val_weights_3+c_val_weights_3)*  del_phi_x_3')'))/((del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')'  *(del_phi_x_3'*f_3-0.5*del_phi_x_3'*g_3*Ru_3^-1*g_3'*del_phi_x_3*(val_weights_3+c_val_weights_3)')+1)^2)+  0.0*beta_3*del_phi_x_3'*g_3*Ru_3^-1*g_3'*[x5;x6])';
            
            val_weights_3 = val_weights_3+a_val_weights_3/6+b_val_weights_3/3+c_val_weights_3/3+d_val_weights_3/6;
            i_hyb_count_3=i_hyb_count_3+1;
            
        end
        control_1(k)=u1;
        control_2(k)=u2;
        control_3(k)=u3;
    end
    
    
    
    %%
    if monte ==1
        t = 0:h:r;
        figure(1);
        subplot(2,1, 1)
        plot(t, X(:,1),'--b','LineWidth',1);hold on;
        plot(t, X(:,3)-theta_2_star,'--r','LineWidth',1);hold on;
        plot(t, X(:,5)-theta_3_star,'--m','LineWidth',1);hold on;
        xlabel('time','FontWeight','b','FontSize',12);ylabel('states','FontWeight','b','FontSize',12);
        title('State Trajectory','FontWeight','b','FontSize',12);grid on
        set( gca, 'FontWeight', 'b','FontSize', 12 );
        legend('x1','x2','x3')
        subplot(2,1, 2)
        plot(t(1:end-1), (control_1),'--b','LineWidth', 1);hold on;
        plot(t(1:end-1), (control_2),'--r','LineWidth', 1);hold on;
        plot(t(1:end-1),(control_3),'--m','LineWidth', 1);grid on;box on;
        title('Control input','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Force applied on the cart','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        
        figure(4);
        subplot(2,1,1)
        plot(t(1:end-1), 100*(abs(Hamiltonian_1)+abs(Hamiltonian_2)+abs(Hamiltonian_3)),'--b','LineWidth', 1);hold on;
        title('Estimated Hamiltonian','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Bellman error','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        subplot(2,1,2)
        plot(t, 100*(Reward_1+Reward_2+Reward_3),'--b','LineWidth', 1);hold on;
        title('Cumulative cost','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Cost','FontWeight', 'b', 'FontSize', 12)
        legend('ss1','ss2','ss3')
        hold on;set( gca, 'FontWeight', 'b', 'FontSize', 12 );
    elseif monte==2
        t = 0:h:r;
        figure(1);
        subplot(2,1, 1)
        plot(t, X(:,1),'b','LineWidth',1);hold on;
        plot(t, X(:,3)-theta_2_star,'r','LineWidth',1);hold on;
        plot(t, X(:,5)-theta_3_star,'m','LineWidth',1);hold on;
        xlabel('time','FontWeight','b','FontSize',12);ylabel('states','FontWeight','b','FontSize',12);
        title('State Trajectory','FontWeight','b','FontSize',12);grid on
        set( gca, 'FontWeight', 'b','FontSize', 12 );
        legend('x1','x2','x3')
        subplot(2,1, 2)
        plot(t(1:end-1), (control_1),'b','LineWidth', 1);hold on;
        plot(t(1:end-1), (control_2),'r','LineWidth', 1);hold on;
        plot(t(1:end-1),(control_3),'m','LineWidth', 1);grid on;box on;
        title('Control input','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Force applied on the cart','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        
        figure(4);
        subplot(2,1,1)
        plot(t(1:end-1), 100*(abs(Hamiltonian_1)+abs(Hamiltonian_2)+abs(Hamiltonian_3)),'b','LineWidth', 1);hold on;
        title('Estimated Hamiltonian','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Bellman error','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        subplot(2,1,2)
        plot(t, 100*(Reward_1+Reward_2+Reward_3),'b','LineWidth', 1);hold on;
        title('Cumulative cost','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Cost','FontWeight', 'b', 'FontSize', 12)
        legend('ss1','ss2','ss3')
        hold on;set( gca, 'FontWeight', 'b', 'FontSize', 12 );
    else
        t = 0:h:r;
        figure(1);
        subplot(2,1, 1)
        plot(t, X(:,1:2),'--b','LineWidth',1);hold on;
        plot(t, X(:,3:4)-theta_2_star,'--r','LineWidth',1);hold on;
        plot(t, X(:,5:6)-theta_3_star,'--m','LineWidth',1);hold on;
        xlabel('time','FontWeight','b','FontSize',12);ylabel('states','FontWeight','b','FontSize',12);
        title('State Trajectory','FontWeight','b','FontSize',12);grid on
        set( gca, 'FontWeight', 'b','FontSize', 12 );
        legend('x1','x2','x3')
        subplot(2,1, 2)
        plot(t(1:end-1), (control_1),'--b','LineWidth', 1);hold on;
        plot(t(1:end-1), (control_2),'--r','LineWidth', 1);hold on;
        plot(t(1:end-1),(control_3),'--m','LineWidth', 1);grid on;box on;
        title('Control input','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Force applied on the cart','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        
        figure(4);
        subplot(2,1,1)
        plot(t(1:end-1), 100*(abs(Hamiltonian_1)+abs(Hamiltonian_2)+abs(Hamiltonian_3)),'--b','LineWidth', 1);hold on;
        title('Estimated Hamiltonian','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Bellman error','FontWeight', 'b', 'FontSize', 12)
        hold on;legend('ss1','ss2','ss3')
        set( gca, 'FontWeight', 'b', 'FontSize', 12 );
        subplot(2,1,2)
        plot(t, 100*(Reward_1+Reward_2+Reward_3),'--b','LineWidth', 1);hold on;
        title('Cumulative cost','FontWeight', 'b', 'FontSize', 12)
        xlabel('time in sec','FontWeight', 'b', 'FontSize', 12)
        ylabel('Cost','FontWeight', 'b', 'FontSize', 12)
        legend('ss1','ss2','ss3')
        hold on;set( gca, 'FontWeight', 'b', 'FontSize', 12 );
    end
    tot_reward(:,monte)=(Reward_1+Reward_2+Reward_3);
    tot_hamil(:,monte)=abs(Hamiltonian_1)+abs(Hamiltonian_2)+abs(Hamiltonian_3);
end
% ratio_diff(monte)=(Reward_1(end)+Reward_2(end)+Reward_3(end))/save_normalizing_reward;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% figure(7);
% subplot(2,2,1);
% stem(t(1:end-1), event_indicating_function_1,'b');hold on;
% plot(t(1:end-1), event_threshold_1,'b','LineWidth',2);
% subplot(2,2,2);stem(t(1:end-1), event_indicating_function_2,'g');
% hold on;plot(t(1:end-1), event_threshold_2,'g','LineWidth',2);
% subplot(2,2,3);stem(t(1:end-1), event_indicating_function_3,'r');
% hold on;plot(t(1:end-1), event_threshold_3,'r','LineWidth',2);
% set( gca, 'FontWeight', 'b','FontSize', 12 );
% grid on;
% box on;
% xlim([0 10])
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% figure(7);subplot(224)
% plot(time_base1*h, inter_event_time1*h,'*b');hold on;
% plot(time_base2*h, inter_event_time2*h,'*g');hold on;
% plot(time_base3*h, inter_event_time3*h,'*r')
% set( gca, 'FontWeight', 'b','FontSize', 12 );
% grid on;box on;
% xlabel('Time','FontWeight', 'b','FontSize', 12);
% ylabel('Inter-event time','FontWeight', 'b','FontSize', 12);
% title('Inter-event time','FontWeight', 'b','FontSize', 12)
% legend('ss1','ss2','ss3')
% tuned_hybrid_ss1=val_weights_1;
% tuned_hybrid_ss2=val_weights_2;
% tuned_hybrid_ss3=val_weights_3;
% tuned_hybrid = [tuned_hybrid_ss1 tuned_hybrid_ss2 tuned_hybrid_ss3];
