%Function callling for Fuzzy Climate Chane Control Scheme

function xdot = e_closed_dyn_ex_biped(t,x)

global u1 u2 u3
theta_2_star = 0;-0.227;
theta_3_star = 0;0.559;
p2 = 6070;
p3 = 192;

% theta_2_star = -0.559;
% theta_3_star = 0.226;
% p2 = 226;
% p3 = 5240;
% f1= x(2);
% f3=x(4);
% f5=x(6);
% f2=0.1*(1-5.25*(x(1))^2)*x(2)-x(1)+0.01*x(6)*(x(5)-theta_3_star)+0.01*x(4)*(x(3)-theta_2_star);
% f4=(0.01*(1-p2*(x(3)-theta_2_star)^2)*x(4))-(4*(x(3)-theta_2_star))+(0.057*x(1)*x(2))+(0.1*(x(4)-x(6)));
% f6=(0.01*(1-p3*(x(5)-theta_3_star)^2)*x(6))-(4*(x(5)-theta_3_star))+(0.057*x(1)*x(2))+(0.1*(x(6)-x(4)));
% f = [f1;f2;f3;f4;f5;f6];
% g = [0;1;0;1;0;1];

xdot(1)=x(2);
xdot(2)= 0.1*(1-5.25*(x(1))^2)*x(2)-x(1)+0.01*x(6)*(x(5)-theta_3_star)+0.01*x(4)*(x(3)-theta_2_star)+u1;
xdot(3)=x(4);
xdot(4)=(0.01*(1-p2*(x(3)-theta_2_star)^2)*x(4))-(4*(x(3)-theta_2_star))+(0.057*x(1)*x(2))+(0.1*(x(4)-x(6)))+u2;
xdot(5)=x(6);
xdot(6)=(0.01*(1-p3*(x(5)-theta_3_star)^2)*x(6))-(4*(x(5)-theta_3_star))+(0.057*x(1)*x(2))+(0.1*(x(6)-x(4)))+u3;
xdot=xdot';
