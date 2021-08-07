V=16.5;
I=2.72;
time=[79.2; 77.47; 68.32; 60.03; 53.43; 46.28];
x_ac=[0; 0.2; 0.4; 0.6; 0.8; 1];
cp_w = 4.186;
m = [150*0.997; 133.8; 126.52; 122.41; 119.76; 117.9];
k=(2.7*16.3*time(1) - cp_w*m(1)*2)/2
cp_ac = (2.72*V.*time(6) - k*2)./(2.*m(6))


n=[8.31; 5.14; 3.72; 2.92; 2.39; 2.01];
dT_m = [0; 5.5; 2.9; 0.1; -1.9; 0];

cp_mix = (2.72*V.*time - k*2)./(2.*m)

Q_m = -I*V*time/2 .* dT_m
H_m = Q_m./n;
H_m;

x=0:0.2:1;
xx=0:.01:1;
yy = spline(x,H_m,xx);
yy_p = spline(x,H_m,x);
plot(x,H_m,'o',xx,yy,'-',x,yy_p,'-.', [0; 1], [0; 0], '-k')
ylabel("Value of enthalpy of mixing (H_m) in J/mole")
xlabel("Mole fraction of acetone (x_a)")
%title("H_m VS x_a plot")
legend(["Data points" "Spline fit"],'Location', 'Best')
grid on

H_m_0_x1 = (yy(2)-yy(1))/0.01;
H_m_1_x1 = (yy(101)-yy(100))/0.01;
H_m_5_x1 = (yy(51)-yy(50))/0.01

H_m_5 = yy(51)

H_1p = H_m_5 + (0.5)*H_m_5_x1
H_2p = H_m_5 - (0.5)*H_m_5_x1


H_1inf = H_m_0_x1
H_2inf = -H_m_1_x1