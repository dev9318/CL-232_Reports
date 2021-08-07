d_m = 0.1;
d_I = 0.01;
d_V = 0.01;
d_t = 0.167;
d_T = 0.1;

m = [133.8; 126.52; 122.41; 119.76; 117.9];
V=16.5;
I=2.72;
t=[77.47; 68.32; 60.03; 53.43; 46.28];
T= [5.5; 2.9; 0.1; -1.9; 0];

d_k_rel = d_m./(150*0.997) + d_I./2.7 + d_V/16.3  + d_t./79.2 + abs(d_T./2)

d_cp = d_k_rel + d_m./m + d_I./I + d_V/V + d_t./t + abs(d_T./2)

d_Hm = d_k_rel + d_cp + abs(d_T./T) + d_m./m

H_m = 1.0e3*[-1.8602;
   -1.1952;
   -0.0461;
    0.9532;
    1];
H_m.*d_Hm
Cp_mix = [4.6461;
    3.2906;
    1.8813;
    0.6863;
    0];
Cp_mix .* d_cp 

150*0.997;