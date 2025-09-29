clear all
close all
clc
j = sqrt(-1);
% Enter data for forming admittance matrix
data = [1 2 0.02 0.06 
    1 3 0.08 0.24 
    2 3 0.06 0.25
    2 4 0.06 0.18
    2 5 0.04 0.12
    3 4 0.01 0.03
    4 5 0.08 0.24];
nl = data(:,1);                % From which bus
nr = data(:,2);                % To which bus
R = data(:,3);                 % Line resistance
X = data(:,4);                 % Line reactance
nbr = length(data(:,1));       % Number of branch
nbus = max(max(nl), max(nr));   % Number of bus
Z = R + j*X;
y = ones(nbr,1)./Z;             % ./ performs element wise division
Ybus = zeros(nbus, nbus);
% Off Diagonal element
for k = 1:nbr
    if  not(nl(k) == nr(k)) 
        Ybus(nl(k),nr(k)) = Ybus(nl(k), nr(k)) - y(k);
        Ybus(nr(k),nl(k)) = Ybus(nl(k),nr(k));
    end
end
% Diagonal element element
for n= 1 : nbus
    for k = 1 : nbr
        if nl(k) == n || nr(k) ==n
            Ybus(n,n) = Ybus(n,n) + y(k);
        
        end
    end
end
Ybus 


% For Load Bus Calculation

zdt = [1 0 0 0 0 1.06+j*0
    2 .2 .10 .40 .30 1+j*0
    3 .45 .15 0 0 1+j*0
    4 .40 .05 0 0 1+j*0
    5 .60 .10 0 0 1+j*0];
bn = zdt(:,1);                      % bus number
Pl = zdt(:,2);                      % MW Load
Ql = zdt(:,3);                      % MVAR Load
Pg = zdt(:,4);                      % MW Generator
Qg = zdt(:,5);                      % MVAR Generator
V = zdt(:,6);                       % Bus Voltage


for k=1:4
    S(k) = (zdt(k+1,4) - zdt(k+1,2)) + (1j)*(zdt(k+1,5) - zdt(k+1,3));
    
end
S;
S_conj = conj(S)

V_conj = conj(V);
iter = 0;
    max_iter = 500;
while iter < max_iter
    
for i = 2:5
    sum = 0;
        for k = 1:nbus
            if k ~= i                          
                sum = sum + Ybus(i,k)*V(k);
                sum;
            end
            sum;
        end
        V(i) = (1/Ybus(i,i))*((S_conj(i-1)/V_conj(i)) -sum);
        V_conj(i) = conj(V(i));
        
end
    iter = iter +1;
end
V
V_conj = conj(V);
[theta, Voltage] = cart2pol(real(V), imag(V));
Voltage
Phase_Angle = rad2deg(theta);
Phase_Angle

%Voltage_conj = conj(Voltage);
s = 0;

for i = 1:1:5
    for k = 1:1:5
        if i~=k
            s = s + (Ybus(i,k)*V(k));
        end
    end
    S(i) = (V_conj(i)*((s) + (V(i)*Ybus(i,i))));
end
S;
       
Power_Slack_Bus = (real(S(1))*100)
Q_Slack_Bus = (imag(S(1))*100)
z = conj(Z(1));

sq = ((V(1)-V(2))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_1 = conj(Z(1));
p_loss_Line_One_Two1 = ((mag_1/z_1)*100)
p_loss_Line_One_Two = real(p_loss_Line_One_Two1)
Q_loss_Line_One_Two = imag(p_loss_Line_One_Two1)

sq = ((V(1)-V(3))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_2 = conj(Z(2));
p_loss_Line_One_THREE1 = ((mag_1/z_2)*100);
p_loss_Line_One_THREE = real (p_loss_Line_One_THREE1)
Q_loss_Line_One_THREE = imag(p_loss_Line_One_THREE1)

sq = ((V(2)-V(3))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_3 = conj(Z(3));
p_loss_Line_tWO_THREE1 = ((mag_1/z_3)*100);
p_loss_Line_tWO_THREE = real(p_loss_Line_tWO_THREE1)
Q_loss_Line_tWO_THREE = imag(p_loss_Line_tWO_THREE1)

sq = ((V(2)-V(4))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_4 = conj(Z(4));
p_loss_Line_tWO_fOUR1 = ((mag_1/z_4)*100);
p_loss_Line_tWO_fOUR = real(p_loss_Line_tWO_fOUR1)
Q_loss_Line_tWO_fOUR = imag(p_loss_Line_tWO_fOUR1)

sq = ((V(2)-V(5))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_5 = conj(Z(5));
p_loss_Line_TWO_fIVE1 = ((mag_1/z_5)*100)
p_loss_Line_TWO_fIVE = real(p_loss_Line_TWO_fIVE1)
Q_loss_Line_TWO_fIVE = imag(p_loss_Line_TWO_fIVE1)

sq = ((V(3)-V(4))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_6 = conj(Z(6));
p_loss_Line_tHREE_fOUR1 = ((mag_1/z_6)*100);
p_loss_Line_tHREE_fOUR = real(p_loss_Line_tHREE_fOUR1)
Q_loss_Line_tHREE_fOUR = imag(p_loss_Line_tHREE_fOUR1)

sq = ((V(4)-V(5))^2);
[theta, mag_1] = cart2pol(real(sq), imag(sq));
z_7 = conj(Z(7));
p_loss_Line_fOUR_fIVE1 = ((mag_1/z_7)*100);
p_loss_Line_fOUR_fIVE = real(p_loss_Line_fOUR_fIVE1)
Q_loss_Line_fOUR_fIVE = imag(p_loss_Line_fOUR_fIVE1)