x_cal = 0:0.1:1;
RI_cal = [ 1.355855 1.355955 1.355754;
        1.368868 1.368868 1.368167;
        1.386585 1.387086 1.387987;
        1.404603 1.403502 1.404803;
        1.42092 1.420019 1.42172;
        1.44174 1.44164 1.444843;
        1.447446 1.449348 1.447947;
        1.459058 1.45155 1.459258;
        1.475474 1.475975 1.475774;
        1.488988 1.487886 1.488887;
        1.496195 1.497196 1.496495;
    ];

RI_mean = (RI_cal(:,1) + RI_cal(:,2) + RI_cal(:,3))/3;

X = ones(11,2);
X(:,2) = x_cal;
coef = X\RI_mean;

RI_raw = [ 1.380679 1.379178 1.379478 1.386886 1.385985 1.389488;
        1.399298 1.400199 1.399999 1.408007 1.41181 1.412411;
        1.408707 1.404103 1.404904 1.438637 1.44094 1.440539;
        1.411911 1.418517 1.418918 1.448347 1.445945 1.447846;
        1.435834 1.436135 1.433933 1.453852 1.454553 1.455154;
        1.447146 1.450049 1.45105 1.477476 1.478177 1.477076;
    ];

RI_v = (RI_raw(:,1) + RI_raw(:,2) + RI_raw(:,3))/3;

RI_l = (RI_raw(:,4) + RI_raw(:,5) + RI_raw(:,6))/3;

V_v = (RI_v - coef(1))/coef(2);

V_l = (RI_l - coef(1))/coef(2);

rho_A = 729;
rho_B = 879;
M_A = 58.08;
M_B = 78.11;

Y_A = 1 - (V_v.*rho_B./M_B)./(V_v.*rho_B./M_B + (1-V_v).*rho_A./M_A)
X_A = 1 - (V_l.*rho_B/M_B)./(V_l.*rho_B/M_B + (1-V_l).*rho_A/M_A)

K_A = Y_A./X_A;
K_B = (1-Y_A)./(1-X_A)

P=1; %1bar or approximately 1 atm
A_A = 14.3145;
B_A = 2756.22;
C_A = 228.060;

T = [58.2; 59.9; 60.8; 63.0; 65.2; 68.1];

P_As = exp(A_A - B_A./(T+C_A))/100;

A_B = 13.7819;
B_B = 2726.81;
C_B = 217.572;

P_Bs = exp(A_B - B_B./(T+C_B))/100;


gamma_A = K_A.*P./P_As;
gamma_B = K_B.*P./P_Bs;
ln_gamma_A = log(gamma_A)
ln_gamma_B = log(gamma_B)

Y_ln = log(gamma_B./gamma_A);

R = 8.314;

G_E = R.*(273+T).*(X_A.*ln_gamma_A + (1-X_A).*ln_gamma_B)
G_E_RT = G_E./(R.*T);

x = zeros(8,1);
x(2:7) = X_A;
x(1) = 1;
xx = 0:0.01:1;
y_g_A = zeros(7,1);
y_g_A(2:7) = ln_gamma_A;

y_g_B = zeros(7,1);
y_g_B(1:6) = ln_gamma_B;
spline_y_A = spline(x(1:7), y_g_A, xx);
spline_y_B = spline(x(2:8), y_g_B, xx);


% plot(xx, spline_y_A)
% hold on
% plot(xx, spline_y_B)
% plot(x(1:7),y_g_A, 'ko')
% plot(x(2:8),y_g_B, 'ko')
% plot([0 1], [0 0], 'k-')
% legend('log (\gamma_A)', 'log (\gamma_B)')
% ylabel('log (\gamma_A) and log (\gamma_B)')
% xlabel("Mole fraction of acetone (X_A)")
% hold off

Y_ln = spline_y_A - spline_y_B;

% uncomment this code for graph 2
hold on
plot(xx(2:100),Y_ln(2:100))
xlabel("Mole fraction of acetone (X_A)")
ylabel("log (\gamma_1/\gamma_2)")
plot([0 1], [0 0], 'k-')
legend(["spline fit for log (\gamma_1/\gamma_2)"; ""])
hold off

area = trapz(xx(2:100),Y_ln(2:100))

% Van Larr
a = ln_gamma_A.*(1 + ((1-X_A).*ln_gamma_B)./(X_A.*ln_gamma_A)).^2
b = ln_gamma_B.*(1 + (X_A.*ln_gamma_A)./((1-X_A).*ln_gamma_B)).^2