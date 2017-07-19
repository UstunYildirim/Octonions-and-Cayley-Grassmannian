LIB "primdec.lib";
option(redSB);
ring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;
poly p1 = cn+dg-hk-ij;
poly p2 = -goa+ima+kla-lna+bgp-bhm+bjl+bmn-c+ceo-cfp+chk-chn+cij-ckn-del+dfm-dgh-dgk-eip-eko+fhp-fjo+h+hkn+ijn+k-n;
poly p3 = -goa+ima+kla-lna+bgp-bhm+bjl+bmn+c+ceo-cfp+chk-chn+cij-ckn-del+dfm-dgh-dgk-eip-eko+fhp-fjo-h+hkn+ijn-k+n;
poly p4 = -ka+na+bg-bj-cf+cp+de-dm-ei+fh-go-hp+im+jo+kl-ln;
poly p5 = ka-na+bg+bj-cf-cp-de+dm-ei+fh-go+hp+im-jo+kl-ln;
poly p6 = ga-ia+bh-bn-ce-co-df+dl+ek+fj-gp+hm+ip-jl+ko-mn;
poly p7 = ga+ia-bh+bn-ce+co+df-dl+ek+fj-gp+hm-ip-jl-ko-mn;
ideal I = p1,p2,p3,p4,p5,p6,p7;

I = groebner(I);
matrix H = jacob(I);
ideal J = I,(minor(H,4));

J = groebner(J);
J = radical(J);


matrix M = jacob(J);
ideal K = J,(minor(M,11));
K = groebner (K);
print(J);
print(K);

