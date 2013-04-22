syms c0 c1 c2
syms aWAT a0 a1 a2
syms bbWAT bb0 bb1 bb2
syms bWAT b0 b1 b2
syms h al al1 al2
syms mu0 qf
syms RW0 RW1 RW2
syms KD0 KD1
syms G1 G2 L0 ll


syms s

a = aWAT + a0 * c0 + a1 * c1 + a2 * c2;
bb = bbWAT + bb0 * c0 + bb1 * c1 + bb2 * c2;
b = bbWAT / bWAT + bb0 * c0 / b0 + bb1 * c1 / b1 + bb2 * c2 / b2;
rw = RW0 + RW1 * bb / a + RW2 * (bb / a)^2;
kd = sqrt(a^2 + a * b * (KD0 + KD1 * mu0)) / mu0;
al = 0.25*al1 * exp(- (ll - G1)^2 / (2* G2^2)) + L0 + 0.0005*al2*ll;

rt = rw * (1 - exp(-2 * kd * h)) + al * exp(-2 * kd * h) / qf;
tdeep = (rw - s);
tshal = (rt - s);



%jdeep = diff(tdeep, c0);
%ccode(jdeep)


jshal = diff(tshal, c0);
ccode(jshal)


jshal_al1 = diff(tshal, al1);
ccode(jshal_al1)

jshal_al2 = diff(tshal, al2);
ccode(jshal_al2)
